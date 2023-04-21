#!/usr/bin/env Rscript

# Bna_ensembl_ID_gene_mappings.R

# Author : Natalia Garcia Sanchez
# Date 21/04/2023
# Description : Get external attribute (e.g. gene name ENA - Ensembl mapping) through BioMART API programmatic access
# ---
args = commandArgs(trailingOnly = FALSE)

print("~~~~~~~~~~~~ Bna_ensembl_gene_mappings.R ~~~~~~~~~~~~")
writeLines(print("\nGet Brassica napus Ensembl Gene ID mappings to external attribute (e.g. gene name ENA - Ensembl mapping) through BioMART API programmatic access\n"))
print("To see a list of all attributes available in Ensembl Biomart, go to https://plants.ensembl.org/biomart/martview/")
print ("or perform the Bioconductor-biomaRt function `listAttributes()` on `bnapus` mart dataset object")
print("-------------------")
print("Script usage : Rscript Bna_ensembl_gene_mappings.R <table_wth_ensembl_ID.csv> <BIOMART_FEATURE_ATTRIBUTE> <output_table_wth_mappings.csv>")
print("---")

if(length(args)<3){stop("You have to provide three arguments in order : -INPUT -ATTRIBUTE -OUTPUT", call.=FALSE)}
fname <- toString(args[1])
attribute <- toString(args[2])
fname_out<- toString(args[3])

##################################################
# Installing and loading required packages
##################################################

print("Installing and loading required packages through Biomanager(biomaRt, biomartr, tidyverse)")
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("biomartr")
install.packages("tidyverse")}

# Load libraries
library("tidyverse")
library("biomaRt")
library("biomartr")

###########################################################################
# Create an object with the directory containing HTseq counts:
###########################################################################


#"/home/famgarcia/TFM/data/Testillano col/Report RNAseq MV vs PE/DE_ANALYSIS_SARA/ALL/Final_condition_PErep_vs_VMrep_total.csv"
MV_vsPE_deseq_generes <- read.table(fname)
index_g_id = which(colnames(MV_vsPE_deseq_generes) == "ensembl_gene_id")
head(MV_vsPE_deseq_generes[,which(colnames(MV_vsPE_deseq_generes) == "ensembl_gene_id")])
rownames(MV_vsPE_deseq_generes) <- MV_vsPE_deseq_generes[,index_g_id]

#########################
# Query the mart/dataset
#########################

# see a list of "marts" available at host "plants.ensembl.org"
listMarts <-listMarts(host="https://plants.ensembl.org")

# create an object 
mart <- useMart(biomart="plants_mart", host="https://plants.ensembl.org")

# see a list of data sets within the mart

listdatasets <- listDatasets(mart)
(listDatasets(mart))[grep("bnapus",listDatasets(mart)[,1]),]

# create an object for the plants_mart-Ensembl 

bnapus_mart <- useMart(biomart = "plants_mart", host = "https://plants.ensembl.org", dataset = "bnapus_eg_gene"	)

# Uncomment to see a list of all "attributes" available
# attributes <- listAttributes(mart = bnapus_mart, page="feature_page")
# attributes

# get gene names and transcript lengths when they exist
ann <- getBM(filter="ensembl_gene_id",value=rownames(MV_vsPE_deseq_generes),attributes=c("ensembl_gene_id",attribute), mart=bnapus_mart, useCache = FALSE)

MV_vsPE_deseq_generes <-merge(MV_vsPE_deseq_generes,ann,by='ensembl_gene_id',all.x = TRUE, all.y = TRUE)
#"/home/famgarcia/TFM/data/Testillano col/Report RNAseq MV vs PE/DE_ANALYSIS_SARA/ALL/Final_condition_PErep_vs_VMrep_total_gene_name_annotated.csv"
write.csv(MV_vsPE_deseq_generes, fname_out)
