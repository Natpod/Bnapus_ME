#!/usr/bin/env Rscript

# Normalize batch adjusted data
# Author: Natalia García Sánchez

# Variance Stabilizing Transformation with DESeq2

data_to_norm <- read.csv("/home/famgarcia/Escritorio/genes_all_experiments_batch_correction.csv", sep=",")


# create dds to perform a Variance Stabilizing Transformation Normalization from the counts
dds_tonorm <- DESeqDataSetFromMatrix(countData = data_to_norm,
                                     colData = colData,
                                     design = ~ 1) # not spcifying model
dds_tonorm<-vst(dds_tonorm)
write.table(assay(dds_tonorm),"/home/famgarcia/Escritorio/VST_norm_countscsv")


# filter by gene id
data_to_norm <- assay(dds_tonorm)
data_to_norm$ensembl_id = rownames(data_to_norm)

# 2. Filter by having Entrez ID

filter <- read.csv("/home/famgarcia/Descargas/mart_export.csv", sep=",")
colnames(filter)<-c("ensembl_id","gid")
filter = filter[!is.na(filter$gid),]

# Inner join

data_to_norm = merge(data_to_norm, filter, by = "ensembl_id")

# prepare data
data_to_norm = data_to_norm[!duplicated(data_to_norm$gid),]
data_to_norm$ensembl_id = NULL
rownames(data_to_norm) = data_to_norm$gid
data_to_norm$gid = NULL

write.table(data_to_norm,"/home/famgarcia/Escritorio/VST_norm_countscsv")
