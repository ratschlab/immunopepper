import numpy as np
import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects import pandas2ri, Formula, r
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
pandas2ri.activate()
deseq = importr('DESeq2')
BiocParallel = importr('BiocParallel')
BiocGenerics = importr("BiocGenerics")


def DESeq2(count_matrix, design_matrix, normalize, cores=1):
    # gene_column = ''
    to_dataframe = ro.r('function(x) data.frame(x)')

    count_matrix = count_matrix.astype(float)
    count_matrix[count_matrix.isnull()] = 0.0
    count_matrix = round(count_matrix)

    count_matrix = pandas2ri.py2rpy(count_matrix)
    design_matrix = pandas2ri.py2rpy(design_matrix)

    #count_matrix.rownames = count_matrix.get_attrib('kmer')
    design_formula = Formula(' ~ 1')


    dds0 = deseq.DESeqDataSetFromMatrix(countData=count_matrix,
                                            colData=design_matrix,
                                            design=design_formula)
    dds0 = BiocGenerics.estimateSizeFactors(dds0, type="poscounts")
    order_size_factor = list(dds0.do_slot('colData').do_slot('rownames'))
    dds0.do_slot('colData').do_slot('listData')[1] = ro.vectors.FloatVector(list(normalize.loc[order_size_factor, 'libsize_75percent'])) # Enforce size factors



    dds = deseq.DESeq(dds0, parallel=True, BPPARAM=BiocParallel.MulticoreParam(cores),
                      sfType="poscounts", # Will run 1. estimation of size factors: estimateSizeFactors # parameter "poscounts"
                      fitType="parametric" # 2. estimation of dispersion: estimateDispersions # parameter "parametric"
                      )

    #normalized_count_matrix = deseq.counts(self, normalized=Tru os.path.join(cancer_dir, 'expression_counts.libsize.tsv'),e)
    deseq_result = deseq.results(dds)
    fit_res = to_dataframe(deseq_result)
    disp = to_dataframe(deseq.dispersions(dds)).rename({'x': 'dispersion'}, axis = 1)
    disp.index = fit_res.index
    fit_res = pd.concat([fit_res['baseMean'], disp], axis=1)


    return fit_res

def mode_cancerspecif(arg):
    normal_matrix = pd.read_parquet(arg.normal_matrix, engine = 'pyarrow')
    normal_matrix  = normal_matrix.drop(['is_cross_junction'], axis=1)
    normal_matrix = normal_matrix.set_index('kmer')
    design_matrix = pd.DataFrame([1]* normal_matrix.shape[1], columns = ["design"])
    design_matrix['sample'] = normal_matrix.columns
    design_matrix = design_matrix.set_index('sample')
    libsize = pd.read_csv(arg.normal_libsize, sep = '\t')
    libsize['libsize_75percent'] = libsize['libsize_75percent'] / np.median(libsize['libsize_75percent'])
    libsize = libsize.set_index('sample')
    DESeq2(normal_matrix, design_matrix, normalize=libsize, cores=1)



#with localconverter(ro.default_converter + pandas2ri.converter):
  #r_from_pd_df = ro.conversion.py2rpy(normal_matrix)

