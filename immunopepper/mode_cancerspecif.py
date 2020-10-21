
import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
pandas2ri.activate()
deseq = importr('DESeq2')
BiocParallel = importr('BiocParallel')

def DESeq2(count_matrix, design_matrix, id_column='kmer', cores=1):
    # gene_column = ''
    to_dataframe = ro.r('function(x) data.frame(x)')

    count_matrix = count_matrix.astype(float)
    count_matrix[count_matrix.isnull()] = 0.0
    count_matrix = round(count_matrix)

    count_matrix = pandas2ri.py2rpy(count_matrix)
    design_matrix = pandas2ri.py2rpy(design_matrix)
    design_formula = Formula(' ~ 1')

    dds = deseq.DESeqDataSetFromMatrix(countData=count_matrix,
                                            colData=design_matrix,
                                            design=design_formula)

    dds = deseq.DESeq(dds, parallel=True, BPPARAM=BiocParallel.MulticoreParam(cores),
                      sfType="poscounts", # Will run 1. estimation of size factors: estimateSizeFactors # parameter "poscounts"
                      fitType="parametric" # 2. estimation of dispersion: estimateDispersions # parameter "parametric"
                      )

    #normalized_count_matrix = deseq.counts(self, normalized=True)
    deseq_result = deseq.results(dds)
    deseq_result = to_dataframe(deseq_result) #already pandas dataframe
    #deseq_result = pandas2ri.rpy2py(deseq_result)  ## back to pandas dataframe

    return deseq_result

def mode_cancerspecif(arg):
    normal_matrix = pd.read_parquet(arg.normal_matrix, engine = 'pyarrow')
    normal_matrix  = normal_matrix.drop(['is_cross_junction'], axis=1)
    normal_matrix = normal_matrix.set_index(['kmer'])
    design_matrix = pd.DataFrame(['normal']* normal_matrix.shape[1], columns = ["design"])
    design_matrix.index  = normal_matrix.columns

    DESeq2(normal_matrix, design_matrix, id_column='kmer', cores=1)



#with localconverter(ro.default_converter + pandas2ri.converter):
  #r_from_pd_df = ro.conversion.py2rpy(normal_matrix)

