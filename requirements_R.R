#!/usr/bin/env Rscript
# Install required packages

# requirements
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# annotation
BiocManager::install("biomaRt")
BiocManager::install("biomartr")
remotes::install_github("GangLiLab/genekitr")
BiocManager::install("GO.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("AnnotationHub")


# Differential expression
BiocManager::install("DESeq2")
BiocManager::install("apeglm")

# Gene set enrichment
BiocManager::install("clusterProfiler")
BiocManager::install("GSVA")


# batch effect correction
devtools::install_github("zhangyuqing/sva-devel")
BiocManager::install("BatchQC")


# data handling
install.packages("tidyverse")
install.packages("stringr")
install.packages("dplyr")
install.packages("readxl")
install.packages(data.table)

# plotting
install.packages("ggplot2")
install.packages("scales")
install.packages("ggpubr")
install.packages("reshape2")
install.packages("gridExtra")
install.packages("pheatmap")
install.packages("igraph")
install.packages(c("gplots", "amap"))
devtools::install_github("kevinblighe/EnhancedVolcano")
BiocManager::install("apeglm")
install.packages("dendextend")
devtools::install_github("gaospecial/ggVennDiagram")
install.packages("Rgraphviz")
install.packages("ggraph")
BiocManager::install("KEGGgraph")
BiocManager::install("pathview")
BiocManager::install("rrvgo")


# WGCNA
remotes::install_github("kevinblighe/CorLevelPlot")
install.packages("WGCNA", dependencies = TRUE)
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("minet")
install.packages("venneuler")
remotes::install_github("jtlovell/limmaDE2") 
