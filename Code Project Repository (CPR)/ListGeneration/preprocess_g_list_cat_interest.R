#!/bin/usr/env Rscript

# Author: Natalia García Sánchez
# Description : index crafted gene category association file for terms of interest, get gene-category count.
# Objective numb 1 : get counts of categories of interest
# Objective numb 2 : filter lists of DEGs by gene-category associations of interest.
# Objective numb 3 : generate GSEA results for categories of interest


##### Load packages

library("readxl")
library("data.table")
library("tidyverse")
library("stringr")
library("clusterProfiler")
library("genekitr")



######################################################
# Load bnapus object
######################################################

library("AnnotationDbi")
library("AnnotationHub")
library("readxl")

# Making annotation package
ah <- AnnotationHub()
#Ensembl 56
query(ah, c( "OrgDb", "Brassica napus"))
org.Brassicanapus.eg.db <- ah[["AH107389"]]


# custom function to unlist genes
unlist_deg_genes <- function(x) {unlist(strsplit(x, "/"))}

###################### LOAD ANNOTATIONS #############################################

ensembl_c_annotations = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\ensembl_centric_annotations.tsv", sep="\t")
ncbi_c_annotations = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\ncbi_centric_annotations.tsv", sep="\t")
PANTHER_annot = read_excel("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\PANTHER_proteins.xlsx")
#TPM_values = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\TPM_values_bcorrection_DevStage.csv")
TPM_values = read.table("D:\\VST_norm_counts.csv")
TPM_values$ensembl_id =  rownames(TPM_values)

VST_values = read.table("C:\\Users\\naata\\MASTER\\Enrichment\\VST_norm_counts_no_bce.csv")
VST_values$ensembl_id = rownames(VST_values)
map_ensembl_name = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\bna_ID_to_name.csv")
VST_values = merge(VST_values, map_ensembl_name, by="ensembl_id")
VST_values$ensembl_id<-NULL


# # external_gene_name (ENA) to GID
# mappings = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\bna_EID.csv")
# map_ensembl_name = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\bna_ID_to_name.csv")
# TPM_values = merge(TPM_values, map_ensembl_name, by="ensembl_id")
# TPM_values$ensembl_id <- NULL
# write.table(TPM_values, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\VST_external_name_all.csv")


logFC_HDAC_Stress<-read.table("C:\\Users\\naata\\MASTER\\Enrichment\\Condition_T_vs_VM_result_padj_0.05_DevStage.csv", sep= ",", header=TRUE)
logFC_HDAC<-read.table("C:\\Users\\naata\\MASTER\\Enrichment\\Final_Condition_T_vs_C_LFC1_padj05.tsv", sep= "\t", header=TRUE)
logFC_Stress<-read_csv("C:\\Users\\naata\\MASTER\\Enrichment\\Condition_PEall_vs_VM_result_padj_0.05_DevStage.csv")



###################### LOAD TERMS OF INTEREST #############################################

hormone_terms = read_excel("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\final_list_GO_lookup.xlsx", sheet="hormonal_screening")
epi_terms = read_excel("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\final_list_GO_lookup.xlsx", sheet="epigenetic_screening")

##################### IMPORT AND PREPROCESS RMD FILES WITH GENE-CATEGORY ASSOCIATION TUPLES #######################

# Biological Process

ml_BP_GO <- readRDS("C:\\Users\\naata\\Documents\\MASTER\\BioSynth\\synthetic_env\\GSEA\\full_bp_genesets_noorth.rds")

# Handling pivoted multilevel lists and converting into table variables

# GO_BP
listindex <- seq(ml_BP_GO)
pathwayindex <- names(ml_BP_GO)

# Pass to df format
tmp_list <- lapply(listindex, function(i){
  gl = list(ml_BP_GO[[i]])[[1]]
  cbind(rep(names(ml_BP_GO)[i],length(gl)),gl)})

pathways_BP_gs<-do.call(rbind,  tmp_list)
pathways_BP_gs<-as.data.frame(pathways_BP_gs)
colnames(pathways_BP_gs)=c("bp", "gene")


# Load pathway terms

pathways_BP_terms <-read.csv('C:\\Users\\naata\\Documents\\MASTER\\BioSynth\\synthetic_env\\GSEA\\full_bp_terms.csv')


##################################################################################################################

# get gene-category associations for categories of interest

pathways_BP_gs_horm = pathways_BP_gs[pathways_BP_gs$bp %in% hormone_terms$ID, ]
pathways_BP_gs_epi = pathways_BP_gs[pathways_BP_gs$bp %in% epi_terms$ID, ]

################################## OBJ-1 .Annotate pathway count in categories of interest

hormone_terms_count = as.data.frame(table(pathways_BP_gs_horm$bp))
hormone_terms_count = merge(pathways_BP_terms, hormone_terms_count, by.x="bp", by.y="Var1")
hormone_terms_count = hormone_terms_count[order(hormone_terms_count$Freq),]
colnames(hormone_terms_count) <-c("ID", "name", "Count")
hormone_terms_count = merge(hormone_terms, hormone_terms_count, by="ID")[,c(1,2,3,5)]
#write.table(hormone_terms_count, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\hormone_finalterms_report.tsv", sep="\t", quote=F, row.names=F)

epigenetic_terms_count = as.data.frame(table(pathways_BP_gs_epi$bp))
epigenetic_terms_count = merge(pathways_BP_terms, epigenetic_terms_count, by.x="bp", by.y="Var1")
epigenetic_terms_count = epigenetic_terms_count[order(epigenetic_terms_count$Freq),]
colnames(epigenetic_terms_count) <-c("ID", "name", "Count")
epigenetic_terms_count = merge(epi_terms, epigenetic_terms_count, by="ID")[,c(1,2,3,5)]
#write.table(epigenetic_terms_count, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\epigenetic_finalterms_report.tsv", sep="\t", quote=F, row.names=F)


################################## OBJ-2 .GSEA, category specific

##################### LOAD and preprocess ranking files

ranks_HDAC_Stress <- read.csv("C:\\Users\\naata\\MASTER\\Enrichment\\Rank_T_vs_VM.rnk",
                              header=FALSE, colClasses = c("character", "numeric"))


ranks_HDAC <- read.csv("C:\\Users\\naata\\Documents\\MASTER\\BioSynth\\synthetic_env\\GSEA\\bna_HDAC_deseq2stat.rnk",
                       header=FALSE, colClasses = c("character", "numeric"))


ranks_Stress <- read.csv("C:\\Users\\naata\\MASTER\\Enrichment\\Rank_adj_VM_vs_PEall.csv",
                         header=FALSE, colClasses = c("character", "numeric"))


# Remove duplicates (transcript isoform with more differential expression will be conserved from the same gene)
ranks_HDAC_Stress <- ranks_HDAC_Stress[!duplicated(ranks_HDAC_Stress[,1]),]
ranks_HDAC <- ranks_HDAC[!duplicated(ranks_HDAC[,1]),]
ranks_Stress <- ranks_Stress[!duplicated(ranks_Stress[,1]),]


# Labelling columns (gene_name "ID" and statistic "t")
colnames(ranks_HDAC_Stress) <-c("ID","t")
colnames(ranks_HDAC) <-c("ID","t")
colnames(ranks_Stress) <-c("ID","t")


# Sorting and creating ordered gene ranks column
ranks_HDAC_Stress <- ranks_HDAC_Stress %>% mutate(rank = rank(t,  ties.method = "random"))
ranks_HDAC_Stress <- ranks_HDAC_Stress[order(ranks_HDAC_Stress$rank, decreasing=TRUE), ]

ranks_HDAC <- ranks_HDAC %>% mutate(rank = rank(t,  ties.method = "random"))
ranks_HDAC <- ranks_HDAC[order(ranks_HDAC$rank, decreasing=TRUE), ]

ranks_Stress <- ranks_Stress %>% mutate(rank = rank(t,  ties.method = "random"))
ranks_Stress <- ranks_Stress[order(ranks_Stress$rank, decreasing=TRUE), ]


# Converting to set of lists
ranks_HDAC_Stress <- setNames(ranks_HDAC_Stress$t, ranks_HDAC_Stress$ID)
ranks_HDAC <- setNames(ranks_HDAC$t, ranks_HDAC$ID)
ranks_Stress <- setNames(ranks_Stress$t, ranks_Stress$ID)


#################### PERFORM GSEA


############# EPIGENETIC TERMS

gsea_HDAC_Stress_EPI_O <- clusterProfiler::GSEA(geneList = ranks_HDAC_Stress,nPerm =10000, minGSSize=5,maxGSSize=2000, eps=0.0, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs_epi, TERM2NAME = pathways_BP_terms, seed=TRUE, by="fgsea")
gsea_HDAC_Stress_EPI<-gsea_HDAC_Stress_EPI_O[,1:11]


gsea_HDAC_EPI_O <- clusterProfiler::GSEA(geneList = ranks_HDAC,nPerm =10000, minGSSize=5,maxGSSize=2000, eps=0.0, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs_epi, TERM2NAME = pathways_BP_terms, seed=TRUE, by="fgsea")
gsea_HDAC_EPI<-gsea_HDAC_EPI_O[,1:11]

gsea_Stress_EPI_O <- clusterProfiler::GSEA(geneList = ranks_Stress,nPerm =10000, minGSSize=5,maxGSSize=2000, eps=0.0, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs_epi, TERM2NAME = pathways_BP_terms, seed=TRUE, by="fgsea")
gsea_Stress_EPI<-gsea_Stress_EPI_O[,1:11]

############# HORMONE TERMS
gsea_HDAC_Stress_H_O <- clusterProfiler::GSEA(geneList = ranks_HDAC_Stress,nPerm =10000, minGSSize=5,maxGSSize=2000, eps=0.0, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs_horm, TERM2NAME = pathways_BP_terms, seed=TRUE, by="fgsea")
gsea_HDAC_Stress_H<-gsea_HDAC_Stress_H_O[,1:11]

gsea_HDAC_H_O <- clusterProfiler::GSEA(geneList = ranks_HDAC,nPerm =10000, minGSSize=5,maxGSSize=2000, eps=0.0, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs_horm, TERM2NAME = pathways_BP_terms, seed=TRUE, by="fgsea")
gsea_HDAC_H<-gsea_HDAC_H_O[,1:11]

gsea_Stress_H_O <- clusterProfiler::GSEA(geneList = ranks_Stress,nPerm =10000, minGSSize=5,maxGSSize=2000, eps=0.0, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs_horm, TERM2NAME = pathways_BP_terms, seed=TRUE, by="fgsea")
gsea_Stress_H<-gsea_Stress_H_O[,1:11]

write.table(gsea_HDAC_Stress_EPI, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\GSEA_episcat_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)
write.table(gsea_HDAC_EPI, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\GSEA_episcat_PE_SAHA.tsv", sep="\t", row.names=F, quote=F)
write.table(gsea_Stress_EPI, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\GSEA_episcat_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)
write.table(gsea_HDAC_Stress_H, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\GSEA_hormonescat_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)
write.table(gsea_HDAC_H, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\GSEA_hormonescat_PE_SAHA.tsv", sep="\t", row.names=F, quote=F)
write.table(gsea_Stress_H, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\GSEA_hormonescat_HeatStress.tsv", sep="\t", row.names=F, quote=F)




################################# DEG ANNOTATION
epigenetic_terms_count <- epigenetic_terms_count[order(epigenetic_terms_count$Count),]
hormone_terms_count <- hormone_terms_count[order(hormone_terms_count$Count),]
  
df_terms_ep = data.frame(ID=c("_"), gene=c("_"))
for (i in epigenetic_terms_count$ID){
  # index in more specific terms, get genes in order and annotate GO terms in order
  genes = pathways_BP_gs_epi[pathways_BP_gs_epi$bp == i,2]
  genes = genes[!(genes %in% intersect(genes, df_terms_ep[,2]))]
  tmp_df = cbind(rep(i, length(genes)), genes)
  colnames(tmp_df) <- c("ID", "gene")
  df_terms_ep <- rbind(df_terms_ep, tmp_df) 
  
}

df_terms_h = data.frame(ID=c("_"), gene=c("_"))
for (i in hormone_terms_count$ID){
  # index in more specific terms, get genes in order and annotate GO terms in order
  genes = pathways_BP_gs_horm[pathways_BP_gs_horm$bp == i,2]
  genes = genes[!(genes %in% intersect(genes, df_terms_h[,2]))]
  tmp_df = cbind(rep(i, length(genes)), genes)
  colnames(tmp_df) <- c("ID", "gene")
  df_terms_h <- rbind(df_terms_h, tmp_df) 
  
}
############### HDACi STRESS

genes_epi_HDAC_Stress = unique(unlist(lapply(gsea_HDAC_Stress_EPI$core_enrichment, unlist_deg_genes)))
genes_epi_HDAC_Stress <- data.frame(external_gene_name = genes_epi_HDAC_Stress)
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, PANTHER_annot, by="external_gene_name")
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress,df_terms_ep,by.x="external_gene_name", by.y="gene")
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, epi_terms[,1:3], by = "ID")
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, VST_values, by="external_gene_name")
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, logFC_HDAC_Stress[,c(4,9)], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)

genes_epi_HDAC_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", round(genes_epi_HDAC_Stress$`log2FoldChange.x`,2), sep="")
genes_epi_HDAC_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", round(genes_epi_HDAC_Stress$`log2FoldChange.y`,2), sep="")
genes_epi_HDAC_Stress <- genes_epi_HDAC_Stress[!duplicated(genes_epi_HDAC_Stress$external_gene_name),]
genes_epi_HDAC_Stress <- merge(genes_epi_HDAC_Stress, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC_Stress <- merge(genes_epi_HDAC_Stress, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)

genes_epi_HDAC_Stress<-genes_epi_HDAC_Stress[,c(1,26,27,9,10,11,12,24,25,13:23)]


write.table(genes_epi_HDAC_Stress, "C:\\Users\\naata\\Downloads\\Leading_edges_epigenetic_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)

#write.table(genes_epi_HDAC_Stress, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_epigenetic_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)



genes_hormone_HDAC_Stress = unique(unlist(lapply(gsea_HDAC_Stress_H_O$core_enrichment, unlist_deg_genes)))
genes_hormone_HDAC_Stress <- data.frame(external_gene_name = genes_hormone_HDAC_Stress)
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, PANTHER_annot, by="external_gene_name")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress,df_terms_h,by.x="external_gene_name", by.y="gene")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, hormone_terms[,1:3], by = "ID")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, VST_values, by="external_gene_name")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, logFC_HDAC_Stress[,c(4,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", round(genes_hormone_HDAC_Stress$`log2FoldChange.x`,2), sep="")
genes_hormone_HDAC_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", round(genes_hormone_HDAC_Stress$`log2FoldChange.y`,2), sep="")
genes_hormone_HDAC_Stress <- genes_hormone_HDAC_Stress[!duplicated(genes_hormone_HDAC_Stress$external_gene_name),]
genes_hormone_HDAC_Stress <- merge(genes_hormone_HDAC_Stress, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC_Stress <- merge(genes_hormone_HDAC_Stress, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)

genes_hormone_HDAC_Stress<-genes_hormone_HDAC_Stress[,c(1,26,27,9,10,11,12,24,25,13:23)]
write.table(genes_hormone_HDAC_Stress, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_hormone_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)



############## STRESS

genes_epi_Stress = unique(unlist(lapply(gsea_Stress_EPI$core_enrichment, unlist_deg_genes)))
genes_epi_Stress <- data.frame(external_gene_name = genes_epi_Stress)
genes_epi_Stress = merge(genes_epi_Stress, PANTHER_annot, by="external_gene_name")
genes_epi_Stress = merge(genes_epi_Stress,df_terms_ep,by.x="external_gene_name", by.y="gene")
genes_epi_Stress = merge(genes_epi_Stress, epi_terms[,1:3], by = "ID")
genes_epi_Stress = merge(genes_epi_Stress, VST_values, by="external_gene_name")
genes_epi_Stress = merge(genes_epi_Stress, logFC_HDAC_Stress[,c(4,9)], by="external_gene_name", all.x=TRUE)
genes_epi_Stress = merge(genes_epi_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)

genes_epi_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", round(genes_epi_Stress$`log2FoldChange.x`,2), sep="")
genes_epi_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", round(genes_epi_Stress$`log2FoldChange.y`,2) , sep="")
genes_epi_Stress <- genes_epi_Stress[!duplicated(genes_epi_Stress$external_gene_name),]
genes_epi_Stress <- merge(genes_epi_Stress, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_epi_Stress <- merge(genes_epi_Stress, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)

genes_epi_Stress<-genes_epi_Stress[,c(1,26,27,9,10,11,12,24,25,13:23)]
write.table(genes_epi_Stress, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_epigenetic_HeatStress.tsv", sep="\t", row.names=F, quote=F)



genes_hormone_Stress = unique(unlist(lapply(gsea_Stress_H_O$core_enrichment, unlist_deg_genes)))
genes_hormone_Stress = unique(unlist(lapply(gsea_Stress_H_O$core_enrichment, unlist_deg_genes)))
genes_hormone_Stress <- data.frame(external_gene_name = genes_hormone_Stress)
genes_hormone_Stress = merge(genes_hormone_Stress, PANTHER_annot, by="external_gene_name")
genes_hormone_Stress = merge(genes_hormone_Stress,df_terms_h,by.x="external_gene_name", by.y="gene")
genes_hormone_Stress = merge(genes_hormone_Stress, hormone_terms[,1:3], by = "ID")
genes_hormone_Stress = merge(genes_hormone_Stress, VST_values, by="external_gene_name")
genes_hormone_Stress = merge(genes_hormone_Stress, logFC_HDAC_Stress[,c(4,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_Stress = merge(genes_hormone_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", round(genes_hormone_Stress$`log2FoldChange.x`,2), sep="")
genes_hormone_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", round(genes_hormone_Stress$`log2FoldChange.y`,2), sep="")
genes_hormone_Stress <- genes_hormone_Stress[!duplicated(genes_hormone_Stress$external_gene_name),]
genes_hormone_Stress <- merge(genes_hormone_Stress, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_hormone_Stress <- merge(genes_hormone_Stress, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)

genes_hormone_Stress<-genes_hormone_Stress[,c(1,26,27,9,10,11,12,24,25,13:23)]
write.table(genes_hormone_HDAC_Stress, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_hormone_HeatStress.tsv", sep="\t", row.names=F, quote=F)


############## HDACi proembryos

genes_epi_HDAC = unique(unlist(lapply(gsea_HDAC_EPI$core_enrichment, unlist_deg_genes)))

genes_epi_HDAC <- data.frame(external_gene_name = genes_epi_HDAC)
genes_epi_HDAC = merge(genes_epi_HDAC, PANTHER_annot, by="external_gene_name")
genes_epi_HDAC = merge(genes_epi_HDAC,df_terms_ep,by.x="external_gene_name", by.y="gene")
genes_epi_HDAC = merge(genes_epi_HDAC, epi_terms[,1:3], by = "ID")
genes_epi_HDAC = merge(genes_epi_HDAC, VST_values, by="external_gene_name")
genes_epi_HDAC = merge(genes_epi_HDAC, logFC_HDAC[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC$log2FoldChange = paste("LFC PE-SAHA: ", round(genes_epi_HDAC$log2FoldChange,2), sep="")
genes_epi_HDAC <- genes_epi_HDAC[!duplicated(genes_epi_HDAC$external_gene_name),]
genes_epi_HDAC <- merge(genes_epi_HDAC, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC <- merge(genes_epi_HDAC, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC <- genes_epi_HDAC[,c(1,25,26,9,10,11,12, 24,13:23)]
write.table(genes_epi_HDAC, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_epi_PE_SAHA.tsv", sep="\t", row.names=F, quote=F)



genes_hormone_HDAC = unique(unlist(lapply(gsea_HDAC_H_O$core_enrichment, unlist_deg_genes)))

genes_hormone_HDAC <- data.frame(external_gene_name = genes_hormone_HDAC)
genes_hormone_HDAC = merge(genes_hormone_HDAC, PANTHER_annot, by="external_gene_name")
genes_hormone_HDAC = merge(genes_hormone_HDAC,df_terms_h,by.x="external_gene_name", by.y="gene")
genes_hormone_HDAC = merge(genes_hormone_HDAC, hormone_terms[,1:3], by = "ID")
genes_hormone_HDAC = merge(genes_hormone_HDAC, VST_values, by="external_gene_name")
genes_hormone_HDAC = merge(genes_hormone_HDAC, logFC_HDAC[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC$log2FoldChange = paste("LFC PE-SAHA: ", round(genes_hormone_HDAC$log2FoldChange,2), sep="")
genes_hormone_HDAC <- genes_hormone_HDAC[!duplicated(genes_epi_HDAC$external_gene_name),]
genes_hormone_HDAC <- merge(genes_hormone_HDAC, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC <- merge(genes_hormone_HDAC, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC <- genes_hormone_HDAC[,c(1,25,26,9,10,11,12, 24,13:23)]
write.table(genes_hormone_HDAC, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_hormone_PE_SAHA.tsv", sep="\t", row.names=F, quote=F)


##################################### DEG intersection annotation

attr_order = c("external_gene_name","ENTREZID", "ALIAS", "Description_Ensembl", "Description_NCBI", "Uniprot_ID", "PANTHER_Family_name", "PANTHER_SF_name", "GO", "GOTERM", "ID", "term", "context","lfc_HeatStress_SAHA", "lfc_HeatStress", "lfc_PE_SAHA", "Crep1", "Crep2", "Crep3", "Trep1", "Trep2", "Trep3", "VMrep1", "VMrep2", "VMrep3", "PErep1", "PErep2")

#################################### EPIGENETIC TERMS

############# Outersection - unique DEGS in SAHA Embryogenesis


DEGs_relevant = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\outersec_Embryogenesis_in_HeatStress_SAHA.tsv", sep="\t")

DEGs_relevant_f_epi = merge(DEGs_relevant, data.frame(external_gene_name = pathways_BP_gs_epi[,2]), by="external_gene_name")
DEGS_epi_SAHA = DEGs_relevant_f_epi
DEGS_epi_SAHA = merge(DEGS_epi_SAHA,df_terms_ep,by.x="external_gene_name", by.y="gene",all.x=T)
DEGS_epi_SAHA = merge(DEGS_epi_SAHA, epi_terms, by = "ID", all.x=T)
DEGS_epi_SAHA <- DEGS_epi_SAHA[!duplicated(DEGS_epi_SAHA$external_gene_name),]

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(DEGS_epi_SAHA) %in% i)) }

DEGS_epi_SAHA<-DEGS_epi_SAHA[,orderlist]
write.table(DEGS_epi_SAHA, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_HeatStressSAHAonly_interest_epigenetics_general.tsv", sep="\t", row.names=F, quote=F)


table(DEGS_epi_SAHA[,colnames(DEGS_epi_SAHA)=="context"])
DEGS_epi_SAHA[grepl("methylation", DEGS_epi_SAHA$context),]


##############  Intersection all embryogenesis DEGs genes

DEGs_relevant = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\intersec_Embryogenesis_all.tsv", sep="\t")

DEGs_relevant_f_epi = merge(DEGs_relevant, data.frame(external_gene_name = pathways_BP_gs_epi[,2]), by="external_gene_name")
DEGS_epi_common = DEGs_relevant_f_epi
DEGS_epi_common = merge(DEGS_epi_common,df_terms_ep,by.x="external_gene_name", by.y="gene",all.x=T)
DEGS_epi_common = merge(DEGS_epi_common, epi_terms, by = "ID", all.x=T)
DEGS_epi_common <- DEGS_epi_common[!duplicated(DEGS_epi_common$external_gene_name),]

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(DEGS_epi_common) %in% i)) }

DEGS_epi_common<-DEGS_epi_common[,orderlist]

write.table(DEGS_epi_common, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_Embryogenesisintersection_epigenetics_general.tsv", sep="\t", row.names=F, quote=F)
table(DEGS_epi_common[,colnames(DEGS_epi_common)=="term"])


############## Outersection - unique DEGs in Conventional Embryogenesis

DEGs_relevant = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\outersec_Embryogenesis_in_HeatStress.tsv", sep="\t")

DEGs_relevant_f_epi = merge(DEGs_relevant, data.frame(external_gene_name = pathways_BP_gs_epi[,2]), by="external_gene_name")
DEGS_epi_HS = DEGs_relevant_f_epi
DEGS_epi_HS = merge(DEGS_epi_HS,df_terms_ep,by.x="external_gene_name", by.y="gene",all.x=T)
DEGS_epi_HS = merge(DEGS_epi_HS, epi_terms, by = "ID", all.x=T)
DEGS_epi_HS <- DEGS_epi_HS[!duplicated(DEGS_epi_HS$external_gene_name),]

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(DEGS_epi_HS) %in% i)) }

DEGS_epi_HS<-DEGS_epi_HS[,orderlist]
write.table(DEGS_epi_HS, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_HeatStressonly_interest_epigenetics_general.tsv", sep="\t", row.names=F, quote=F)
table(DEGS_epi_HS[,colnames(DEGS_epi_HS)=="term"])
DEGS_epi_HS[grepl("acetylation", DEGS_epi_HS$term),]



# Saving term frequency statistics in file
summary_epi = merge( as.data.frame(table(DEGS_epi_HS[,colnames(DEGS_epi_HS)=="term"])), as.data.frame(table(DEGS_epi_common[,colnames(DEGS_epi_common)=="term"])), by="Var1", all.x = T, all.y=T)
summary_epi = merge( summary_epi, as.data.frame(table(DEGS_epi_SAHA[,colnames(DEGS_epi_SAHA)=="term"])), by="Var1", all.x = T, all.y=T)
colnames(summary_epi)<-c("GOterm", "HeatStress_freq", "Intersection_freq", "HeatStress_SAHA")
summary_epi[which(is.na(summary_epi[,2])),2]<-0
summary_epi[which(is.na(summary_epi[,3])),3]<-0
summary_epi[which(is.na(summary_epi[,4])),4]<-0
summary_epi$Total <- rowSums(summary_epi[,c(2,3,4)])
summary_epi<- summary_epi[order(summary_epi$Total, decreasing=T),]
write.table(summary_epi, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\DEG_subsets_epigenetics_statistics.tsv", sep="\t", row.names=F, quote=F)


####################### HORMONE TERM FILTERING


############# Outersection - unique DEGS in SAHA Embryogenesis


DEGs_relevant = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\outersec_Embryogenesis_in_HeatStress_SAHA.tsv", sep="\t")

DEGs_relevant_f_hormone = merge(DEGs_relevant, data.frame(external_gene_name = pathways_BP_gs_horm[,2]), by="external_gene_name")
DEGS_hormone_SAHA = DEGs_relevant_f_hormone
DEGS_hormone_SAHA = merge(DEGS_hormone_SAHA,df_terms_h,by.x="external_gene_name", by.y="gene",all.x=T)
DEGS_hormone_SAHA = merge(DEGS_hormone_SAHA, hormone_terms, by = "ID", all.x=T)
DEGS_hormone_SAHA <- DEGS_hormone_SAHA[!duplicated(DEGS_hormone_SAHA$external_gene_name),]

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(DEGS_hormone_SAHA) %in% i)) }

DEGS_hormone_SAHA<-DEGS_hormone_SAHA[,orderlist]

write.table(DEGS_hormone_SAHA, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_HeatStressSAHAonly_interest_hormone_general.tsv", sep="\t", row.names=F, quote=F)


table(DEGS_hormone_SAHA[,colnames(DEGS_hormone_SAHA)=="context"])
table(DEGS_hormone_SAHA[,colnames(DEGS_hormone_SAHA)=="term"])
DEGS_hormone_SAHA[grepl("methylation", DEGS_hormone_SAHA$context),]


##############  Intersection all embryogenesis DEGs genes

DEGs_relevant = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\intersec_Embryogenesis_all.tsv", sep="\t")

DEGs_relevant_f_hormone = merge(DEGs_relevant, data.frame(external_gene_name = pathways_BP_gs_horm[,2]), by="external_gene_name")
DEGS_hormone_common = DEGs_relevant_f_hormone
DEGS_hormone_common = merge(DEGS_hormone_common,df_terms_h,by.x="external_gene_name", by.y="gene",all.x=T)
DEGS_hormone_common = merge(DEGS_hormone_common, hormone_terms, by = "ID", all.x=T)
DEGS_hormone_common <- DEGS_hormone_common[!duplicated(DEGS_hormone_common$external_gene_name),]

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(DEGS_hormone_common) %in% i)) }

DEGS_hormone_common<-DEGS_hormone_common[,orderlist]

write.table(DEGS_hormone_common, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_Embryogenesisintersection_hormone_general.tsv", sep="\t", row.names=F, quote=F)
table(DEGS_hormone_common[,colnames(DEGS_hormone_common)=="term"])


############## Outersection - unique DEGs in Conventional Embryogenesis

DEGs_relevant = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\outersec_Embryogenesis_in_HeatStress.tsv", sep="\t")

DEGs_relevant_f_hormone = merge(DEGs_relevant, data.frame(external_gene_name = pathways_BP_gs_horm[,2]), by="external_gene_name")
DEGS_hormone_HS = DEGs_relevant_f_hormone
DEGS_hormone_HS = merge(DEGS_hormone_HS,df_terms_h,by.x="external_gene_name", by.y="gene",all.x=T)
DEGS_hormone_HS = merge(DEGS_hormone_HS, hormone_terms, by = "ID", all.x=T)
DEGS_hormone_HS <- DEGS_hormone_HS[!duplicated(DEGS_hormone_HS$external_gene_name),]

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(DEGS_hormone_HS) %in% i)) }

DEGS_hormone_HS<-DEGS_hormone_HS[,orderlist]
write.table(DEGS_hormone_HS, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_HeatStressonly_interest_hormone_general.tsv", sep="\t", row.names=F, quote=F)

table(DEGS_hormone_HS[,colnames(DEGS_hormone_HS)=="term"])
DEGS_hormone_HS[grepl("acetylation", DEGS_hormone_HS$term),]



# Saving term frequency statistics in file
summary_hormone = merge( as.data.frame(table(DEGS_hormone_HS[,colnames(DEGS_hormone_HS)=="term"])), as.data.frame(table(DEGS_hormone_common[,colnames(DEGS_hormone_common)=="term"])), by="Var1", all.x = T, all.y=T)
summary_hormone = merge( summary_hormone, as.data.frame(table(DEGS_hormone_SAHA[,colnames(DEGS_hormone_SAHA)=="term"])), by="Var1", all.x = T, all.y=T)
colnames(summary_hormone)<-c("GOterm", "HeatStress_freq", "Intersection_freq", "HeatStress_SAHA")
summary_hormone[which(is.na(summary_hormone[,2])),2]<-0
summary_hormone[which(is.na(summary_hormone[,3])),3]<-0
summary_hormone[which(is.na(summary_hormone[,4])),4]<-0
summary_hormone$Total <- rowSums(summary_hormone[,c(2,3,4)])
summary_hormone<- summary_hormone[order(summary_hormone$Total, decreasing=T),]
write.table(summary_hormone, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\DEG_subsets_hormone_statistics.tsv", sep="\t", row.names=F, quote=F)


DEGS_hormone_SAHA[DEGS_hormone_SAHA$term %in% summary_hormone[(summary_hormone$HeatStress_SAHA>0)&((summary_hormone$HeatStress_freq==0)&(summary_hormone$Intersection_freq == 0)),]$GOterm , ]
DEGS_epi_SAHA[DEGS_epi_SAHA$term %in% summary_epi[(summary_epi$HeatStress_SAHA>0)&((summary_epi$HeatStress_freq==0)&(summary_epi$Intersection_freq == 0)),]$GOterm , ]

DEGS_hormone_SAHA[DEGS_hormone_SAHA$term %in% summary_hormone[(summary_hormone$HeatStress_SAHA>0)&((summary_hormone$HeatStress_freq==0)),1] , ]
DEGS_epi_SAHA[DEGS_epi_SAHA$term %in% summary_epi[(summary_epi$HeatStress_SAHA>0)&((summary_epi$HeatStress_freq==0)),1]  , ]



summary_hormone[(summary_hormone$HeatStress_freq>0)&((summary_hormone$HeatStress_SAHA==0)&(summary_hormone$Intersection_freq == 0)),]
summary_epi[(summary_epi$HeatStress_freq>0)&((summary_epi$HeatStress_SAHA==0)&(summary_epi$Intersection_freq == 0)),]


########## OVERRREPRESENTATION ANALYSIS

bckgr_geneset_HDAC = as.data.frame( names(ranks_HDAC_Stress[(ranks_HDAC_Stress > 0)]) ) 
colnames(bckgr_geneset_HDAC) = "gene"
bckgr_geneset_HDAC <- merge(background_UP_gl_HDAC, pathways_BP_gs, by="gene")
colnames(bckgr_geneset_HDAC) <- c("ID","bp")
bckgr_geneset_HDAC<-bckgr_geneset_HDAC[,c(2,1)]
colnames(DEGs_relevant) <- "gene"

ORA_sDEGs_UP <- clusterProfiler::enricher(gene = DEGs_relevant[21:length(DEGs_relevant[,1]),], minGSSize=2, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = bckgr_geneset_HDAC, TERM2NAME = pathways_BP_terms)
ORA_sDEGs_UP@organism="bnapus"
ORA_sDEGs_UP@ontology="BP"


sDEG_ORA <- importCP(ORA_sDEGs_UP, type = "other")
write_csv(sDEG_ORA,)
sDEG_cinterest = rbind(sDEG_ORA[sDEG_ORA$ID %in% epigenetic_terms_count$ID,], sDEG_ORA[sDEG_ORA$ID %in% hormone_terms_count$ID,])
# Sort by count to get more specific genes first
sDEG_cinterest = sDEG_cinterest[order(sDEG_cinterest$Count),]
# Get genes in list in deg_all, annotate GO terms that were selected in a separate df after processing
list_term_gene = list()
deg_all <- c()
count=0
for ( i in 1:dim(sDEG_cinterest)[1] ){
  newlist <- unlist_deg_genes(sDEG_cinterest[i,9])
  deg_all <- c(deg_all, newlist)
  deg_all <- unique(deg_all)
}
deg_all<-data.frame(external_gene_name=deg_all)
d_all_p = merge(deg_all, PANTHER_annot, by="external_gene_name")[,c(1,8,9)]
d_all_ncbi = merge(deg_all, ncbi_c_annotations, by="external_gene_name")[,c(1,3,4)]
d_all_ncbi = d_all_ncbi[!duplicated(d_all_ncbi$external_gene_name),]
d_all_ensembl = merge(deg_all, ensembl_c_annotations, by="external_gene_name")

all_annot = merge(deg_all, logFC_HDAC_Stress[,c(4,10)], by="external_gene_name")
all_annot = merge(all_annot, logFC_Stress[,c(3,9)], by="external_gene_name")
all_annot = merge(all_annot, d_all_p, by="external_gene_name", all.x=T, all.=T)
all_annot = merge(all_annot, d_all_ncbi, by="external_gene_name", all.x=T, all.=T)
all_annot = merge(all_annot, d_all_ensembl, by="external_gene_name", all.x=T, all.=T)
all_annot = merge(all_annot, TPM_values, by="external_gene_name")
all_annot$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", all_annot$`log2FoldChange.x`, sep="")
all_annot$`log2FoldChange.y` = paste("LFC HeatStress: ", all_annot$`log2FoldChange.y`, sep="")
write.table(all_annot, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_subset_interest_in_episcat.csv", row.names=F, quote=F)






############## WGCNA DEG

# INTERCONVERSION - GID TO ENSEMBL ID

# Behaviour A: negatively correlated in PE, positively correlated in VM: blue 

bluemodule = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\SignNet_genes_modules_interest-blue.txt", header=F)
yellowmodule = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\SignNet_genes_modules_interest-yellow.txt", header=F)
behaviourA = data.frame(entrezgene_id = unique(c(bluemodule[,1], yellowmodule[,1])))


g_bA <- merge(behaviourA, ncbi_c_annotations, by="entrezgene_id", all.x=T)
g_bA <- merge(g_bA, ensembl_c_annotations, by="external_gene_name", all.x=T)
g_bA <- merge(g_bA, PANTHER_annot[,c(1,8,9)], by="external_gene_name", all.x=T)
g_bA <- merge(g_bA, df_terms_ep, by.x="external_gene_name", by.y="gene", all.x=T)
g_bA <- merge(g_bA, df_terms_h, by.x="external_gene_name", by.y="gene", all.x=T)
g_bA <- merge(g_bA, VST_values, by.x="external_gene_name", all.x=T)
g_bA = merge(g_bA, logFC_HDAC_Stress[,c(4,9)], by="external_gene_name", all.x=T)
g_bA = merge(g_bA, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=T)
g_bA$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", round(g_bA$`log2FoldChange.x`,2), sep="")
g_bA$`log2FoldChange.y` = paste("LFC HeatStress: ", round(g_bA$`log2FoldChange.y`,2), sep="")

g_bA_GO = AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(g_bA$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
g_bA_GO = g_bA_GO[g_bA_GO$ONTOLOGY=="BP",]
g_bA_GO= g_bA_GO[!(g_bA_GO$GO == "NA"),]

xx <- as.list(g_bA_GO$GO)

# Get the TERMS for the first element of xx

list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
}
g_bA_GO$GOTERM = list_terms
g_bA_GO= g_bA_GO[!is.na(g_bA_GO$GO),]
g_bA = merge( g_bA, g_bA_GO, by.y="ENTREZID", by.x="entrezgene_id", all.x=TRUE)

g_bA$ID.x = unlist(lapply(g_bA$ID.x, function(x){if (is.na(x)){""}else{ pathways_BP_terms[pathways_BP_terms$bp %in% x,2]}}))
g_bA$ID.y = unlist(lapply(g_bA$ID.y, function(x){if (is.na(x)){""}else{ pathways_BP_terms[pathways_BP_terms$bp %in% x,2]}}))

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(g_bA) %in% i)) }
g_bA <- g_bA[,orderlist]

g_bA_summary = data.frame(table(g_bA$GOTERM))
g_bA_summary = g_bA_summary[order(g_bA_summary$Freq, decreasing=T),]


write.table(g_bA, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_UP_PE_DOWN_VM.tsv", sep="\t", row.names=F, quote=F)



behaviourA = merge(mappings, behaviourA, by="entrezgene_id")
behaviourA$entrezgene_id <- NULL

colnames(behaviourA)<-"gene"
behaviourA = behaviourA[behaviourA$gene %in% pathways_BP_gs$gene,]

ORA_bA <- clusterProfiler::enricher(gene = behaviourA, minGSSize=2, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs, TERM2NAME = pathways_BP_terms)

ORA_bA_selected <- rbind(ORA_bA[ORA_bA[,1] %in% hormone_terms$ID,], ORA_bA[ORA_bA[,1] %in% epigenetic_terms_count$ID,])



# Behaviour B: positively correlated in PE, negatively correlated in VM: blue 

turquoisemodule = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\SignNet_genes_modules_interest-turquoise.txt", header=F)
brownmodule = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\SignNet_genes_modules_interest-brown.txt", header=F)
behaviourB = data.frame(entrezgene_id = unique(c(turquoisemodule[,1], brownmodule[,1])))


attr_order = c("external_gene_name","entrezgene_id", "ALIAS", "Description_Ensembl", "Description_NCBI", "Uniprot_Acc", "PANTHER_Family_name", "PANTHER_SF_name",  "ID.x", "ID.y", "GO", "GOTERM", "log2FoldChange.x", "log2FoldChange.y", "Crep1", "Crep2", "Crep3", "Trep1", "Trep2", "Trep3", "VMrep1", "VMrep2", "VMrep3", "PErep1", "PErep2")


g_bB <- merge(behaviourB, ncbi_c_annotations, by="entrezgene_id", all.x=T)
g_bB <- merge(g_bB, ensembl_c_annotations, by="external_gene_name", all.x=T)
g_bB <- merge(g_bB, PANTHER_annot[,c(1,8,9)], by="external_gene_name", all.x=T)
g_bB <- merge(g_bB, df_terms_ep, by.x="external_gene_name", by.y="gene", all.x=T)
g_bB <- merge(g_bB, df_terms_h, by.x="external_gene_name", by.y="gene", all.x=T)

g_bB_GO = AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(g_bB$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
g_bB_GO = g_bB_GO[g_bB_GO$ONTOLOGY=="BP",]
g_bB_GO= g_bB_GO[!(g_bB_GO$GO == "NA"),]

xx <- as.list(g_bB_GO$GO)

# Get the TERMS for the first element of xx

list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
}
g_bB_GO$GOTERM = list_terms
g_bB_GO= g_bB_GO[!is.na(g_bB_GO$GO),]
g_bB = merge( g_bB, g_bB_GO, by.y="ENTREZID", by.x="entrezgene_id", all.x=TRUE)



g_bB <- merge(g_bB, VST_values, by.x="external_gene_name", all.x=T)
g_bB = merge(g_bB, logFC_HDAC_Stress[,c(4,9)], by="external_gene_name", all.x=T)
g_bB = merge(g_bB, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=T)

gB_summary = data.frame(g_bB$GOTERM)
gB_summary[order(gB_summary$Freq, decreasing = T),]

g_bB$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", round(g_bB$`log2FoldChange.x`,2), sep="")
g_bB$`log2FoldChange.y` = paste("LFC HeatStress: ", round(g_bB$`log2FoldChange.y`,2), sep="")

g_bB$ID.x = unlist(lapply(g_bB$ID.x, function(x){if (is.na(x)){""}else{ pathways_BP_terms[pathways_BP_terms$bp %in% x,2]}}))
g_bB$ID.y = unlist(lapply(g_bB$ID.y, function(x){if (is.na(x)){""}else{ pathways_BP_terms[pathways_BP_terms$bp %in% x,2]}}))

orderlist = c()
for (i in attr_order){orderlist <- c(orderlist, which(colnames(g_bB) %in% i)) }
g_bB <- g_bB[,orderlist]
g_bB<-g_bB[!duplicated(g_bB$external_gene_name),]

write.table(g_bB, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_DOWN_PE_UP_VM.tsv", sep="\t", row.names=F, quote=F)




behaviourB = g_bB[g_bB$external_gene_name %in% pathways_BP_gs$gene,]

ORA_bB <- clusterProfiler::enricher(gene = behaviourB, minGSSize=2, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs, TERM2NAME = pathways_BP_terms)


# SAVE GENES FROM OVERREPRESENTED CATEFORIES

ORA_bB_selected <- rbind(ORA_bB[ORA_bB[,1] %in% hormone_terms$ID,], ORA_bB[ORA_bB[,1] %in% epigenetic_terms_count$ID,])
write.table(ORA_bB[,], "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_UP_PE_DOWN_VM_GO_analysis.tsv", sep="\t", row.names=F, quote=F)
write.table(ORA_bA[,], "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_DOWN_PE_UP_VM_GO_analysis.tsv", sep="\t", row.names=F, quote=F)


summary_WGCNA = merge(g_bA_summary, gB_summary, by="Var1", all.x=T, all.y=T)
colnames(summary_WGCNA) <- c("GOterm","Coexpression_Proembryo", "Coexpression_VM")
summary_WGCNA[which(is.na(summary_WGCNA[,2])),2]<-0
summary_WGCNA[which(is.na(summary_WGCNA[,3])),3]<-0
summary_WGCNA$Total <- rowSums(summary_WGCNA[,c(2,3)])

summary_WGCNA= summary_WGCNA[order(summary_WGCNA$Total, decreasing=T),]
write.table(summary_WGCNA, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\WGCNA_behaviours_summary.tsv", sep="\t", row.names=F, quote=F)

