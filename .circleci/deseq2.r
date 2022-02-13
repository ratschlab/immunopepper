# Installs the R deseq2 package which is needed by mhctools
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
