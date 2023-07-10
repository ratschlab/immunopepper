import pandas as pd
import numpy as np
import glob
import logging
import sys
import os
import shlex
import subprocess

#TODO: add this function to the utils.py?
def explode_immunopepper_coord(mx, coord_col='coord', sep=';'):
    coord_mx = mx[coord_col].str.split(sep, expand=True)

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

    logging.info(">>>>> Running immunopepper in pepquery mode \n")
    args_list = shlex.split(arg.argstring)
    if '-o' not in args_list:
        logging.error('Output directory was not provided in --arg.argstring')
        sys.exit(1)

    if '-i' not in args_list:
        logging.error('Input file was not provided in --arg.argstring')
        sys.exit(1)

    if '-ms' not in args_list:
        logging.error('MS/MS spectra was not provided in --arg.argstring')
        sys.exit(1)

    if '-db' not in args_list:
        logging.error('Reference database was not provided in --arg.argstring')
        sys.exit(1)

    if '-i' in args_list and arg.partitioned_tsv:
        logging.info(">>>>> Extracting whole peptides from the filtered kmers file obtained in cancerspecif mode, and "
                     "provided under {} \n".format(arg.partitioned_tsv))

        if not os.path.exists(arg.metadata_path):
            logging.error("Metadata file {} does not exist or was not provided. Please check --metadata-path".format(arg.metadata_path))
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

        if arg.kmer_type == 'junctions':
            meta_coord = explode_immunopepper_coord(meta, coord_col='modifiedExonsCoord', sep=';')
            meta = pd.concat([meta, meta_coord[['junction_coordinate1', 'junction_coordinate2']]], axis=1)
            meta = meta[meta['isIsolated'] == 0]
            kmers_coord= explode_immunopepper_coord(kmers_filt, coord_col='coord', sep=':')
            kmers_filt = pd.concat([kmers_filt.reset_index(), kmers_coord[['strand', 'junction_coordinate1', 'junction_coordinate2']].reset_index()], axis = 1)
            meta_filt = meta[(meta['junction_coordinate1'].isin(kmers_filt['junction_coordinate1'])) & (meta['junction_coordinate2'].isin(kmers_filt['junction_coordinate2']))]
            peptides = meta_filt['peptide'].unique()

        elif arg.kmer_type == 'segments':
            peptides = meta['peptide'].unique()

        #Remove the peptides that do not contain any of the kmers
        peptides = [pep for pep in peptides if any(kmer in pep for kmer in kmers_filt['kmer'].unique())]
        #Save the peptides in the output directory

        input_peptides_file_idx = [i + 1 for i, j in enumerate(args_list) if j == '-i']
        input_peptides_file = args_list[input_peptides_file_idx[0]]

        with open(input_peptides_file, 'w') as f:
            for i, pep in enumerate(peptides):
                f.write('{}\n'.format(pep))

        logging.info(">>>>> Extracted peptides. Saved to {}/peptides.fa \n".format(arg.output_dir))

    logging.info(">>>>> Launching PepQuery with command {} \n".format(arg.argstring))
    command = 'java '+ '-jar '+ arg.pepquery_software_path + ' '+ arg.argstring
    try:
        subprocess.run(command, check= True, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(">>>>> PepQuery failed with error: {} \n".format(e))
        sys.exit(1)
    logging.info(">>>>> PepQuery finished successfully \n")

    #Now we preprocess the output files to add the filtering information

    output_peptides_file_idx = [i + 1 for i, j in enumerate(args_list) if j == '-o']
    output_peptides_file = args_list[output_peptides_file_idx[0]]

    if f'{output_peptides_file}/psm_rank.txt' not in glob.glob(f'{output_peptides_file}/*'):
        logging.error(">>>>> PepQuery failed to generate the psm_rank.txt file. Please check the output directory")
        sys.exit(1)
    else:
        psm_rank = pd.read_csv(f'{output_peptides_file}/psm_rank.txt', sep = '\t')
        ip_out = pd.DataFrame(columns=['peptide', 'modification', 'spectrum', 'score', 'confident', 'pvalue', 'peptides in the reference database that match better the matched MS/MS spectrum', 'random shuffled peptides that obtained better score than the peptide under study', 'total number of random shuffled peptides studied', 'modified reference proteins that match better the spectra', 'filtering summary'])

        for i, row in psm_rank.iterrows():
            if row['pvalue'] == 100 and row['n_db'] > 0:
                new_row = {'peptide': row['peptide'], 'modification': row['modification'],
                           'spectrum': row['spectrum_title'], 'score': row['score'],
                           'confident': row['confident'], 'pvalue': 'NaN',
                           'peptides in the reference database that match better the matched MS/MS spectrum': row['n_db'],
                           'random shuffled peptides that obtained better score than the peptide under study' : 'NaN',
                           'total number of random shuffled peptides studied' : 'NaN',
                           'modified reference proteins that match better the spectra': 'NaN',
                           'filtering summary': 'Failed at competitive filtering based on reference sequences (step 3).'}
                ip_out.loc[i] = new_row
            elif row['pvalue'] < 100 and row['n_db'] == 0:
                if row['pvalue'] > 0.01:
                    new_row = {'peptide': row['peptide'], 'modification': row['modification'],
                               'spectrum': row['spectrum_title'], 'score': row['score'],
                               'confident': row['confident'], 'pvalue': row['pvalue'],
                               'peptides in the reference database that match better the matched MS/MS spectrum': row['n_db'],
                               'random shuffled peptides that obtained better score than the peptide under study' : row['n_random'],
                               'total number of random shuffled peptides studied' : row['total_random'],
                               'modified reference proteins that match better the spectra': 'NaN',
                               'filtering summary': 'Failed at the statistical evaluation based on random shuffling (step 4). The pvalue is >0.01.'}
                    ip_out.loc[i] = new_row
                elif row['pvalue'] <= 0.01 and row['n_ptm'] != 0:
                    new_row = {'peptide': row['peptide'], 'modification': row['modification'],
                               'spectrum': row['spectrum_title'], 'score': row['score'],
                               'confident': row['confident'], 'pvalue': row['pvalue'],
                               'peptides in the reference database that match better the matched MS/MS spectrum': row['n_db'],
                               'random shuffled peptides that obtained better score than the peptide under study' : row['n_random'],
                               'total number of random shuffled peptides studied': row['total_random'],
                               'modified reference proteins that match better the spectra': row['n_ptm'],
                               'filtering summary': 'Failed at the competitive filtering based on reference proteins with post translational modifications (step 5).'}
                    ip_out.loc[i] = new_row
                elif row['pvalue'] <= 0.01 and row['n_ptm'] == 0 and row['confident'] == 'Yes':
                    new_row = {'peptide': row['peptide'], 'modification': row['modification'],
                               'spectrum': row['spectrum_title'], 'score': row['score'],
                               'confident': row['confident'], 'pvalue': row['pvalue'],
                               'peptides in the reference database that match better the matched MS/MS spectrum': row['n_db'],
                               'random shuffled peptides that obtained better score than the peptide under study' : row['n_random'],
                               'total number of random shuffled peptides studied': row['total_random'],
                               'modified reference proteins that match better the spectra': row['n_ptm'],
                               'filtering summary': 'The peptide passed all the filters and the identified spectra is considered confident'}
                    ip_out.loc[i] = new_row

        confident_categories = ["Yes", "No"]
        ip_out["confident"] = pd.Categorical(ip_out["confident"], categories=confident_categories)
        ip_out.sort_values(by=["confident", "score"], ascending=[True, False], inplace=True)
        ip_out.to_csv(f'{arg.output_dir}/peptides_validated.tsv.gz', sep='\t', index=False, compression='gzip')
        logging.info(">>>>> Processed output file saved to {}/peptides_validated.txt \n".format(arg.output_dir))
        logging.info(">>>>> Finished running immunopepper in pepquery mode  \n")



