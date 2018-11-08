import sys
import os
import scipy as sp
import re


if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <dataset>\n' % sys.argv[0])
    sys.exit(1)
dataset = sys.argv[1]

basedir = '/cluster/work/grlab/projects/coukos_immuno/results_%s/output_peptide/ref' % dataset
filt_tag = '.filt'
metafile = os.path.join(basedir, 'ref_metadata%s_REF.tsv.gz' % filt_tag)
metafile_gtex = os.path.join(basedir, 'ref_metadata%s_GTEX.tsv.gz' % filt_tag)
fasta_file = os.path.join(basedir, 'ref_peptides%s.fa' % filt_tag)

### set filter criteria
if dataset == '20171031':
    celllines = ["OD5P", "EOC24", "MM12", "EOC11", "MM33", "OMM475", "ONVC"]
elif dataset == '20180802':
    celllines = ['C3N-02289', 'C3N-02671', 'EOC04', 'GY150806', 'JM-2', 'JM-4', 'MM170111', 'MM909362', 'MM909432', 'Me275', 'Me290', 'RCC14', 'RCC16', 'T1015A', 'T1185B']
else:
    celllines = []
filter_levels = ['1', '2', '3']

### load metadata
metadata = sp.loadtxt(metafile, dtype='str' , delimiter='\t')
metaheader = metadata[0, :]
metadata = metadata[1:, :]

### load metadata GTEX
metadata_gtex = sp.loadtxt(metafile_gtex, dtype='str' , delimiter='\t')
metaheader_gtex = metadata_gtex[0, :]
metadata_gtex = metadata_gtex[1:, :]

for cl in celllines:

    print('Processing %s' % cl)
    
    ### specific cell line expression 
    cl_idx = sp.where([_.startswith('junc_exp:%s' % cl) for _ in metaheader])[0]
    k_idx1 = (metadata[:, cl_idx].astype('float').max(axis=1) > 0)

    ### novel peptide
    idx = sp.where(metaheader == 'peptide_annotated')[0][0]
    k_idx2 = (metadata[:, idx].astype('int') > 0)
    ### novel junction
    idx = sp.where(metaheader == 'junction_annotated')[0][0]
    k_idx3 = (metadata[:, idx].astype('int') > 0)
    ### non-gtex junction
    #idx = sp.where(metaheader == 'is_in_junction_list')[0][0]
    #k_idx4 = (metadata[:, idx].astype('int') == 0)
    gt_idx = sp.where([_.startswith('junc_exp:gtex') for _ in metaheader_gtex])[0]
    k_idx4 = (metadata_gtex[:, gt_idx].astype('float').max(axis=1) < 3)

    idx = sp.where(metaheader == 'output_id')[0][0]

    for fl in filter_levels:
        print('... filter level %s' % fl)
        if fl == '1':    
            kk_idx = k_idx1
        elif fl == '2':    
            kk_idx = (k_idx1 & k_idx2 & k_idx3)
        elif fl == '3':
            kk_idx = (k_idx1 & k_idx2 & k_idx3 & k_idx4)

        sp.savetxt(re.sub(r'.tsv', '.%s.filt_L%s.tsv' % (cl, fl), metafile), sp.r_[metaheader[sp.newaxis, :], metadata[kk_idx]], fmt='%s', delimiter='\t')
        keep_ids = metadata[kk_idx, idx]

        fasta_out = open(re.sub(r'.fa$', '', fasta_file) + '.%s.filt_L%s.fa' % (cl, fl), 'w')
        i = 0
        is_header = True
        for line in open(fasta_file, 'r'):
            if is_header:
                name = line.strip()[1:]
            else:
                if i < keep_ids.shape[0] and name == keep_ids[i]:
                    fasta_out.write('>' + name + '\n' + line)
                    i += 1
            is_header = not is_header
        fasta_out.close()

