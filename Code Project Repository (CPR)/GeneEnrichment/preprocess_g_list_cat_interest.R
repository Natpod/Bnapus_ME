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
# custom function to unlist genes
unlist_deg_genes <- function(x) {unlist(strsplit(x, "/"))}

###################### LOAD ANNOTATIONS #############################################

ensembl_c_annotations = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\ensembl_centric_annotations.tsv", sep="\t")
ncbi_c_annotations = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\ncbi_centric_annotations.tsv", sep="\t")
PANTHER_annot = read_excel("C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\PANTHER_proteins.xlsx")
TPM_values = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\TPM_values_bcorrection_DevStage.csv")

# external_gene_name (ENA) to GID
mappings = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\bna_EID.csv")


logFC_HDAC_Stress<-read.table("D:\\other\\Condition_T_vs_VM_results_total_bce_DevStage.csv", sep= ",", header=TRUE)
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
write.table(hormone_terms_count, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\hormone_finalterms_report.tsv", sep="\t", quote=F, row.names=F)

epigenetic_terms_count = as.data.frame(table(pathways_BP_gs_epi$bp))
epigenetic_terms_count = merge(pathways_BP_terms, epigenetic_terms_count, by.x="bp", by.y="Var1")
epigenetic_terms_count = epigenetic_terms_count[order(epigenetic_terms_count$Freq),]
colnames(epigenetic_terms_count) <-c("ID", "name", "Count")
epigenetic_terms_count = merge(epi_terms, epigenetic_terms_count, by="ID")[,c(1,2,3,5)]
write.table(epigenetic_terms_count, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\epigenetic_finalterms_report.tsv", sep="\t", quote=F, row.names=F)


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
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, TPM_values, by="external_gene_name")
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, logFC_HDAC_Stress[,c(4,10)], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC_Stress = merge(genes_epi_HDAC_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)

genes_epi_HDAC_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", genes_epi_HDAC_Stress$`log2FoldChange.x`, sep="")
genes_epi_HDAC_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", genes_epi_HDAC_Stress$`log2FoldChange.y`, sep="")
genes_epi_HDAC_Stress <- genes_epi_HDAC_Stress[!duplicated(genes_epi_HDAC_Stress$external_gene_name),]
genes_epi_HDAC_Stress <- merge(genes_epi_HDAC_Stress, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC_Stress <- merge(genes_epi_HDAC_Stress, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)

genes_epi_HDAC_Stress<-genes_epi_HDAC_Stress[,c(1,26,27,9,10,11,12,24,25,13:23)]
write.table(genes_epi_HDAC_Stress, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_epigenetic_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)



genes_hormone_HDAC_Stress = unique(unlist(lapply(gsea_HDAC_Stress_H_O$core_enrichment, unlist_deg_genes)))
genes_hormone_HDAC_Stress <- data.frame(external_gene_name = genes_hormone_HDAC_Stress)
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, PANTHER_annot, by="external_gene_name")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress,df_terms_h,by.x="external_gene_name", by.y="gene")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, hormone_terms[,1:3], by = "ID")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, TPM_values, by="external_gene_name")
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, logFC_HDAC_Stress[,c(4,10)], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC_Stress = merge(genes_hormone_HDAC_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", genes_hormone_HDAC_Stress$`log2FoldChange.x`, sep="")
genes_hormone_HDAC_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", genes_hormone_HDAC_Stress$`log2FoldChange.y`, sep="")
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
genes_epi_Stress = merge(genes_epi_Stress, TPM_values, by="external_gene_name")
genes_epi_Stress = merge(genes_epi_Stress, logFC_HDAC_Stress[,c(4,10)], by="external_gene_name", all.x=TRUE)
genes_epi_Stress = merge(genes_epi_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)

genes_epi_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", genes_epi_Stress$`log2FoldChange.x`, sep="")
genes_epi_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", genes_epi_Stress$`log2FoldChange.y`, sep="")
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
genes_hormone_Stress = merge(genes_hormone_Stress, TPM_values, by="external_gene_name")
genes_hormone_Stress = merge(genes_hormone_Stress, logFC_HDAC_Stress[,c(4,10)], by="external_gene_name", all.x=TRUE)
genes_hormone_Stress = merge(genes_hormone_Stress, logFC_Stress[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_Stress$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", genes_hormone_Stress$`log2FoldChange.x`, sep="")
genes_hormone_Stress$`log2FoldChange.y` = paste("LFC HeatStress: ", genes_hormone_Stress$`log2FoldChange.y`, sep="")
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
genes_epi_HDAC = merge(genes_epi_HDAC, TPM_values, by="external_gene_name")
genes_epi_HDAC = merge(genes_epi_HDAC, logFC_HDAC[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_epi_HDAC$log2FoldChange = paste("LFC PE-SAHA: ", genes_epi_HDAC$log2FoldChange, sep="")
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
genes_hormone_HDAC = merge(genes_hormone_HDAC, TPM_values, by="external_gene_name")
genes_hormone_HDAC = merge(genes_hormone_HDAC, logFC_HDAC[,c(3,9)], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC$log2FoldChange = paste("LFC PE-SAHA: ", genes_hormone_HDAC$log2FoldChange, sep="")
genes_hormone_HDAC <- genes_hormone_HDAC[!duplicated(genes_epi_HDAC$external_gene_name),]
genes_hormone_HDAC <- merge(genes_hormone_HDAC, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC <- merge(genes_hormone_HDAC, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)
genes_hormone_HDAC <- genes_hormone_HDAC[,c(1,25,26,9,10,11,12, 24,13:23)]
write.table(genes_hormone_HDAC, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\Leading_edges_hormone_PE_SAHA.tsv", sep="\t", row.names=F, quote=F)




############# INTERSECT DEGS

DEGs_relevant = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEG_subset_interest_HDACi_ME.csv")

DEGs_relevant_f_epi = merge(DEGs_relevant, pathways_BP_gs_epi, by.x="external_gene_name", by.y="gene")

epigenetic_terms_count[epigenetic_terms_count$ID %in% unique(DEGs_relevant_f_epi$bp),]
# get counts present in subset of interest - epigenetics
epi_DEGrev_ctable = merge(epigenetic_terms_count,as.data.frame(table(DEGs_relevant_f_epi$bp)), by.x="ID", by.y="Var1")
epi_DEGrev_ctable = epi_DEGrev_ctable[order(epi_DEGrev_ctable$Freq),]
epi_DEGrev_ctable[epi_DEGrev_ctable$context == ""]
DEGs_epi_pairs = do.call(rbind, 
                             lapply( unique(epi_DEGrev_ctable$context), 
                                     function(x){
                                       list_index = epi_DEGrev_ctable[epi_DEGrev_ctable$context == x,1]
                                       gl = data.frame(external_gene_name = unique(pathways_BP_gs_epi[pathways_BP_gs_epi$bp %in% list_index, 2]))
                                       gl$context = rep( x, dim(gl)[1])
                                       gl
                                     }
                                     
                             )
)
df_terms = data.frame(ID=c("_"), gene=c("_"))
for (i in epi_DEGrev_ctable$ID){
  # index in more specific terms, get genes in order and annotate GO terms in order
  genes = pathways_BP_gs_epi[pathways_BP_gs_epi$bp == i,2]
  genes = genes[!(genes %in% intersect(genes, df_terms[,2]))]
  tmp_df = cbind(rep(i, length(genes)), genes)
  colnames(tmp_df) <- c("ID", "gene")
  df_terms <- rbind(df_terms, tmp_df) 
  
}

DEGs_epi_pairs_f_ncbi = merge(DEGs_epi_pairs, PANTHER_annot, by="external_gene_name")
DEGs_epi = merge(DEGs_epi,df_terms,by.x="external_gene_name", by.y="gene")
DEGs_epi = merge(DEGs_epi, epi_terms[,1:2], by = "ID")
DEGs_epi = merge(DEGs_epi, TPM_values, by="external_gene_name")
DEGs_epi = merge(DEGs_epi, logFC_HDAC_Stress[,c(4,10)], by="external_gene_name")
DEGs_epi = merge(DEGs_epi, logFC_Stress[,c(3,9)], by="external_gene_name")

DEGs_epi$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", DEGs_epi$`log2FoldChange.x`, sep="")
DEGs_epi$`log2FoldChange.y` = paste("LFC HeatStress: ", DEGs_epi$`log2FoldChange.y`, sep="")

DEGs_epi<-DEGs_epi[,c(1,3,10,11,12,24,25,13:23)]
DEGs_epi<-DEGs_epi[,c(1,3,4,5,17,18,6:16)]
DEGs_epi <- DEGs_epi[!duplicated(DEGs_epi$external_gene_name),]
DEGs_epi <- merge(DEGs_epi, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
DEGs_epi <- merge(DEGs_epi, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)
DEGs_epi<-DEGs_epi[,c(1,2,19,20,3,4,5,6,7,8:18)]
write.table(DEGs_epi, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_subset_interest_epigenetics_general.tsv", sep="\t", row.names=F, quote=F)

##############

####################### HORMONE TERM FILTERING


DEGs_relevant_f_hormone = merge(DEGs_relevant, pathways_BP_gs_horm, by.x="external_gene_name", by.y="gene")

hormone_terms_count[hormone_terms_count$ID %in% unique(DEGs_relevant_f_hormone$bp),]

# get counts present in subset of interest - hormones
hormone_DEGrev_ctable = merge(hormone_terms_count,as.data.frame(table(DEGs_relevant_f_hormone$bp)), by.x="ID", by.y="Var1")
hormone_DEGrev_ctable = hormone_DEGrev_ctable[order(hormone_DEGrev_ctable$Freq),]
DEGs_hormone_pairs = do.call(rbind, 
        lapply( unique(hormone_DEGrev_ctable$context), 
                function(x){
                list_index = hormone_DEGrev_ctable[hormone_DEGrev_ctable$context == x,1]
                gl = data.frame(external_gene_name = unique(pathways_BP_gs_horm[pathways_BP_gs_horm$bp %in% list_index, 2]))
                gl$context = rep( x, dim(gl)[1])
                gl
                }
        
                )
              )
# get list of specific terms
df_terms = data.frame(ID=c("_"), gene=c("_"))
for (i in hormone_terms_count$ID){
  # index in more specific terms, get genes in order and annotate GO terms in order
  genes = pathways_BP_gs_horm[pathways_BP_gs_horm$bp == i,2]
  genes = genes[!(genes %in% intersect(genes, df_terms[,2]))]
  tmp_df = cbind(rep(i, length(genes)), genes)
  colnames(tmp_df) <- c("ID", "gene")
  df_terms <- rbind(df_terms, tmp_df) 
  
}



DEGs_hormone_pairs_f_ncbi = merge(DEGs_hormone_pairs, PANTHER_annot, by="external_gene_name")
DEGs_auxin = DEGs_hormone_pairs_f_ncbi
DEGs_auxin = merge(DEGs_auxin,df_terms,by.x="external_gene_name", by.y="gene")
DEGs_auxin = merge(DEGs_auxin, hormone_terms[,1:2], by = "ID")
DEGs_auxin = merge(DEGs_auxin, TPM_values, by="external_gene_name")
DEGs_auxin = merge(DEGs_auxin, logFC_HDAC_Stress[,c(4,10)], by="external_gene_name")
DEGs_auxin = merge(DEGs_auxin, logFC_Stress[,c(3,9)], by="external_gene_name")

DEGs_auxin$`log2FoldChange.x` = paste("LFC HeatStress-SAHA: ", DEGs_auxin$`log2FoldChange.x`, sep="")
DEGs_auxin$`log2FoldChange.y` = paste("LFC HeatStress: ", DEGs_auxin$`log2FoldChange.y`, sep="")

DEGs_auxin <- DEGs_auxin[!duplicated(DEGs_auxin$external_gene_name),]
DEGs_auxin <- merge(DEGs_auxin, ensembl_c_annotations[,1:2], by="external_gene_name", all.x=TRUE)
DEGs_auxin <- merge(DEGs_auxin, ncbi_c_annotations[,2:3], by="external_gene_name", all.x=TRUE)



DEGs_auxin<-DEGs_auxin[,c(1,3,26,27,10,11,12,24,25,13:23)]


DEGs_auxin<-DEGs_auxin[,c(1,3,4,5,17,18,6:16)]
write.table(DEGs_auxin, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEGs_subset_interest_hormones_general.tsv", sep="\t", row.names=F, quote=F)



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
write_csv(sDEG_ORA, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\DEG_subset_interest_ORA_upregulated.csv")
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

behaviourA = merge(mappings, behaviourA, by="entrezgene_id")
behaviourA$entrezgene_id <- NULL
colnames(behaviourA)<-"gene"
behaviourA = behaviourA[behaviourA$gene %in% pathways_BP_gs$gene,]

ORA_bA <- clusterProfiler::enricher(gene = behaviourA, minGSSize=2, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs, TERM2NAME = pathways_BP_terms)

ORA_bA_selected <- rbind(ORA_bA[ORA_bA[,1] %in% hormone_terms$ID,], ORA_bA[ORA_bA[,1] %in% epigenetic_terms_count$ID,])

deg_all <- c()
count=0
for ( i in 1:dim(ORA_bA_selected)[1] ){
  newlist <- unlist_deg_genes(ORA_bA_selected[i,8])
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


write.table(all_annot, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_DOWN_PE_UP_VM.tsv", sep="\t", row.names=F, quote=F)
write.table(all_annot[!grepl("NA",all_annot$log2FoldChange.x),c(1,2,3,6,10:20)], "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_DOWN_PE_UP_VM_DEGS.tsv", sep="\t", row.names=F, quote=F)


# Behaviour B: positively correlated in PE, negatively correlated in VM: blue 

turquoisemodule = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\SignNet_genes_modules_interest-turquoise.txt", header=F)
brownmodule = read.csv("C:\\Users\\naata\\Documents\\MASTER\\TFM\\SignNet_genes_modules_interest-brown.txt", header=F)
behaviourB = data.frame(entrezgene_id = unique(c(turquoisemodule[,1], brownmodule[,1])))

behaviourB = merge(mappings, behaviourB, by="entrezgene_id")
behaviourB$entrezgene_id <- NULL
colnames(behaviourB)<-"gene"
behaviourB = behaviourB[behaviourB$gene %in% pathways_BP_gs$gene,]

ORA_bB <- clusterProfiler::enricher(gene = behaviourB, minGSSize=2, pvalueCutoff=0.05, pAdjustMethod="BH", gson=NULL, TERM2GENE = pathways_BP_gs, TERM2NAME = pathways_BP_terms)


# SAVE GENES FROM OVERREPRESENTED CATEFORIES

ORA_bB_selected <- rbind(ORA_bB[ORA_bB[,1] %in% hormone_terms$ID,], ORA_bB[ORA_bB[,1] %in% epigenetic_terms_count$ID,])

list_term_gene = list()
deg_all <- c()
count=0
for ( i in 1:dim(ORA_bB_selected)[1] ){
  newlist <- unlist_deg_genes(ORA_bB_selected[i,8])
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
write.table(all_annot, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_UP_PE_DOWN_VM.tsv", sep="\t", row.names=F, quote=F)
write.table(all_annot[!grepl("NA",all_annot$log2FoldChange.x),c(1,2,3,6,10:20)], "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_UP_PE_DOWN_VM_DEGS.tsv", sep="\t", row.names=F, quote=F)
write.table(ORA_bB[,], "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_UP_PE_DOWN_VM_GO_analysis.tsv", sep="\t", row.names=F, quote=F)
write.table(ORA_bA[,], "C:\\Users\\naata\\Documents\\MASTER\\TFM\\WGCNA_behaviour_DOWN_PE_UP_VM_GO_analysis.tsv", sep="\t", row.names=F, quote=F)



