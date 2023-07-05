import pandas as pd
import numpy as np
import glob
import logging
import sys
import os

#TODO: add this function to the utils.py?
def explode_immunopepper_coord(mx, coord_col='coord', sep=';'):
    coord_mx = mx[coord_col].str.split(sep, expand=True)  # 7 min

    dummy_fill = -9999

    if 4 not in coord_mx.columns:
        coord_mx = coord_mx.astype(float)
        coord_mx.fillna(dummy_fill, inplace=True)
        coord_mx = coord_mx.astype(int)

    coord_mx['strand'] = None
    coord_mx.loc[coord_mx[1] < coord_mx[2], 'strand'] = '+'
    coord_mx.loc[coord_mx[1] > coord_mx[2], 'strand'] = '-'

    coord_mx['junction_coordinate1'] = None
    coord_mx['junction_coordinate2'] = None

    #The column 4 would appear if we have 3-exon junctions, rather than just 2exon.
    if 4 not in coord_mx.columns:
        coord_mx[4] = dummy_fill
        coord_mx[5] = dummy_fill

    coord_mx = coord_mx.astype(str)  # 7 min

    coord_mx['+first'] = coord_mx[1] + ':' + coord_mx[2]
    coord_mx['+secon'] = coord_mx[3] + ':' + coord_mx[4]
    coord_mx['-first'] = coord_mx[3] + ':' + coord_mx[0]
    coord_mx['-secon'] = coord_mx[5] + ':' + coord_mx[2]

    coord_mx[4] = ['None' if str(dummy_fill) in jx2 else jx2 for jx2 in coord_mx[4]]
    coord_mx[5] = ['None' if str(dummy_fill) in jx2 else jx2 for jx2 in coord_mx[5]]

    coord_mx.loc[(coord_mx['strand'] == '+'), 'junction_coordinate1'] = coord_mx['+first']
    coord_mx.loc[(coord_mx['strand'] == '-'), 'junction_coordinate1'] = coord_mx['-first']
    coord_mx.loc[(coord_mx['strand'] == '+') & (coord_mx[4] != 'None'), 'junction_coordinate2'] = coord_mx['+secon']
    coord_mx.loc[(coord_mx['strand'] == '-') & (coord_mx[4] != 'None'), 'junction_coordinate2'] = coord_mx['-secon']

    return coord_mx

def mode_pepquery(arg):

    if arg.partitioned_tsv and arg.metadata_path and arg.partitioned_coords_tsv:
        logging.info(">>>>> Extracting whole peptides from the filtered kmers file obtained in cancerspecif mode, and "
                     "provided under {} \n".format(arg.partitioned_tsv))
        if not os.path.exists(arg.metadata_path):
            logging.error("Metadata file {} does not exist. Please check --metadata-path".format(arg.metadata_path))
            sys.exit(1)

        if arg.kmer_type == None:
            logging.error("Please specify the type of kmers provided. Please check --kmer-type")
            sys.exit(1)

        meta = pd.read_csv(arg.metadata_path, sep = '\t')
        meta = meta[~meta['peptide'].isna()]
        meta = meta.set_index(np.arange(len(meta)))

        #Load the filtered kmer file
        input_kmers = glob.glob('{}/*part*'.format(arg.partitioned_tsv))
        if len(input_kmers) == 0:
            logging.error("No file partitions in {}. Please check --kmers-path".format(arg.partitioned_tsv))
            sys.exit(1)
        kmers_filt = pd.concat(map(lambda file: pd.read_csv(file , sep ='\t'), input_kmers))

        input_jx_files= glob.glob('{}/*part*'.format(arg.partitioned_coords_tsv))
        if len(input_jx_files) == 0:
            logging.error("No file partitions in {}. Please check --kmers-coord-path".format(arg.partitioned_coords_tsv))
            sys.exit(1)

        input_jx= pd.concat(map(lambda file: pd.read_csv(file , sep ='\t'), input_jx_files))
        input_jx = input_jx[input_jx['kmer'].isin(kmers_filt['kmer'])]

        if arg.kmer_type == 'junctions':
            meta_coord = explode_immunopepper_coord(meta, coord_col='modifiedExonsCoord', sep=';')
            meta = pd.concat([meta, meta_coord[['junction_coordinate1', 'junction_coordinate2']]], axis=1)
            meta = meta[meta['isIsolated'] == 0]
            input_jx_coord= explode_immunopepper_coord(input_jx, coord_col='coord', sep=':')
            kmers_filt = pd.concat([kmers_filt.reset_index(), input_jx_coord[['strand', 'junction_coordinate1', 'junction_coordinate2']].reset_index()], axis = 1)
            meta_filt = meta[(meta['junction_coordinate1'].isin(kmers_filt['junction_coordinate1'])) & (meta['junction_coordinate2'].isin(kmers_filt['junction_coordinate2']))]
            peptides = meta_filt['peptide'].unique()

        elif arg.kmer_type == 'segments':
            peptides = meta[meta['kmer']]

        #Remove the peptides that do not contain any of the kmers
        peptides = [pep for pep in peptides if any(kmer in pep for kmer in kmers_filt['kmer'].unique())]
        #Save the peptides in the output directory
        with open('{}/peptides.fa'.format(arg.output_dir), 'w') as f:
            for i, pep in enumerate(peptides):
                f.write('{}\n'.format(pep))

        logging.info(">>>>> Extracted peptides. Saved to {}/peptides.fa \n".format(arg.output_dir))

    #NOTE: After this point, the peptides are in the input format necessary for the pepQuery prediction




