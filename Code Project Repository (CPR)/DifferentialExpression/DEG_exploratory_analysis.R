#!/usr/bin/env Rscript

# Visualize and annotate DEGs in the three DEG sets
# Author: Natalia García Sánchez

#######################################################
# load packages
#######################################################

library(tidyverse)
library(gridExtra)
library(reshape)
library(ggplot2)
library(ggpubr)
library("GO.db")
library(ggvenn)
library("ggVennDiagram") # requires c pckg - libudunits2-dev, libgdal-dev, units, sf CRAN packages
library("clusterProfiler")

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

annot_ncbi<-read.csv("/home/famgarcia/Escritorio/ncbi_centric_annotations.tsv", sep="\t")
ensembl_ann<-read.csv("/home/famgarcia/Escritorio/ensembl_centric_annotations.tsv", sep="\t")
panther_ann<-read_excel("/home/famgarcia/Escritorio/PANTHER_proteins.xlsx")

################ Load and preprocess DEG data
# Load DEA data - necessary to filter lists and make coexpressed gene output findings more robust

expr_HDAC <- read.csv("/home/famgarcia/T_vs_C_DEA/Final_Condition_T_vs_C_LFC1_padj05.csv")
expr_HDAC_Stress <- read.csv("/home/famgarcia/VM_vs_HDAC_DEA/HeatStress-SAHA/Condition_T_vs_VM_result_padj_0.05_DevStage.csv")
expr_Stress <- read.csv("/home/famgarcia/Escritorio/Condition_PEall_vs_VM_result_padj_0.05_DevStage.csv")
# expr_Stress <- read.csv("/home/famgarcia/TFM/data/Testillano col/Report RNAseq MV vs PE/DE_ANALYSIS_SARA/woPErep3/Final_condition_PErep_vs_VMrep_padj0.05_woPErep3.csv")
filter_E = read.csv("/home/famgarcia/TFM/data/Testillano col/Report RNAseq MV vs PE/DE_ANALYSIS_SARA/woPErep3/mart_export.txt")
# colnames(filter_E)<-c("ensembl_gene_id", "external_gene_name")
# expr_Stress<-merge(expr_Stress, filter_E, by="ensembl_gene_id")
# write.table(not_common_HS[, !(colnames(not_common_HS) =="X")], "/home/famgarcia/TFM/data/Testillano col/Report RNAseq MV vs PE/DE_ANALYSIS_SARA/woPErep3/Final_condition_PErep_vs_VMrep_padj0.05_woPErep3.tsv", sep="\t", row.names=F)
VST_values = read.table("/home/famgarcia/Escritorio/VST_norm_counts_no_bce.csv")
VST_values$ensembl_gene_id = rownames(VST_values)
VST_values<-merge(VST_values, filter_E, by="ensembl_gene_id")
VST_values$ensembl_gene_id<-NULL

# preprocessing : taking index out, filtering by LFC>1, taking external gene name log2FC and padj only
attr_retrieve = c("external_gene_name", "log2FoldChange", "padj")
expr_HDAC <- expr_HDAC[,colnames(expr_HDAC) %in% attr_retrieve]
expr_HDAC_Stress <- expr_HDAC_Stress[,colnames(expr_HDAC_Stress) %in% attr_retrieve]
expr_Stress <- expr_Stress[,colnames(expr_Stress) %in% attr_retrieve]

# log2FoldChange>2 is considered a DEG
expr_HDAC <- expr_HDAC[abs(expr_HDAC$log2FoldChange) >1, ]
expr_HDAC_Stress <- expr_HDAC_Stress[abs(expr_HDAC_Stress$log2FoldChange) >1, ]
expr_Stress <- expr_Stress[abs(expr_Stress$log2FoldChange) >1, ]

expr_HDAC <- expr_HDAC[!duplicated(expr_HDAC$external_gene_name), ]
expr_HDAC_Stress <- expr_HDAC_Stress[!duplicated(expr_HDAC_Stress$external_gene_name), ]
expr_Stress <- expr_Stress[!duplicated(expr_Stress$external_gene_name), ]


########################################
#  MAIN
###########################################

# Biological question 1: How do my three DEG datasets behave? Common DEGs among:
# - Conventional embryogenesis
# - Embriogenesis + HDACi_SAHA
# - Proembryos with and without SAHA

x <- list(
  PE_SAHA = expr_HDAC[!is.na(expr_HDAC$external_gene_name), 3],
  ME_HeatStress_SAHA = expr_HDAC_Stress[, 3], 
  ME_HeatStress = expr_Stress[ ,3]
)

################# MAKE DEG Venn plot
png("/home/famgarcia/Escritorio/Intersection_DEG_all3subsets.png", width=16, height=12, units="cm", res=300)
ggVennDiagram(x, c( "PE_SAHA", "ME_HeatStress_SAHA","          ME_HeatStress"), edge_size = 0.3,set_size=3, 
              color=c("PE_SAHA" = "lightblue", "ME_HeatStress_SAHA" = "lightcoral", "ME_HeatStress" = "khaki"),
              set_color = c("PE_SAHA" = "steelblue", "ME_HeatStress_SAHA" = "lightcoral", "ME_HeatStress" = "khaki")) +
  scale_fill_gradient(low="white",high = "steelblue") + 
  scale_color_manual(values = c("PE_SAHA" = "lightblue", "ME_HeatStress_SAHA" = "lightcoral", "ME_HeatStress" = "khaki")) +
  theme(legend.position="none", text=element_text(size=4))+ # no legend
  scale_x_continuous(expand = expansion(mult = .5)) # avoid cropping long labels
dev.off()

colnames(expr_HDAC)<-c("lfc_PE_SAHA", "padj_PE_SAHA", "external_gene_name")
colnames(expr_HDAC_Stress)<-c("lfc_HeatStress_SAHA", "padj_PE_SAHA", "external_gene_name")
colnames(expr_Stress)<-c("lfc_HeatStress", "padj_HeatStress", "external_gene_name")

############################################## Write annotated gene intersections in tables

# Discard these annotations, which are redundant or rarely annotated
ann_discard <-c("UniprotAcc","persistent_id", "gene_name", "ONTOLOGY", "PANTHER_SF_ID","name_PC","id_PC" )


##### Common for all
common_all = data.frame(external_gene_name = intersect(intersect(x[[1]], x[[3]]), x[[2]]) )
common_all = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], common_all, by="external_gene_name", all.y=T)
common_all = merge(panther_ann, common_all, by="external_gene_name", all.y=T)
common_all = merge(ensembl_ann, common_all, by="external_gene_name", all.y=T)
common_all = merge(common_all, expr_HDAC_Stress[,c(3,1)], by="external_gene_name", all.x=T)
common_all = merge(common_all, expr_Stress[,c(3,1)], by="external_gene_name", all.x=T)
common_all = merge(common_all, expr_HDAC[,c(3,1)], by="external_gene_name", all.x=T)

common_all_GO= AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(common_all$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
common_all_GO = common_all_GO[common_all_GO$ONTOLOGY=="BP",]
common_all_GO= common_all_GO[!(common_all_GO$GO == "NA"),]
xx <- as.list(common_all_GO$GO)

list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
common_all_GO$GOTERM = list_terms
common_all_GO= common_all_GO[!is.na(common_all_GO$GO),]
common_all = merge(common_all_GO, common_all, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
common_all = common_all[,!(colnames(common_all) %in% ann_discard)]
common_all$lfc_PE_SAHA = paste("LFC PE-SAHA : ", round(common_all$lfc_PE_SAHA,2))
common_all$lfc_HeatStress_SAHA = paste("LFC HeatStress-SAHA : ", round(common_all$lfc_HeatStress_SAHA,2))
common_all$lfc_HeatStress = paste("LFC HeatStress: ", round(common_all$lfc_HeatStress,2))
common_all = merge(common_all, VST_values, by="external_gene_name", all.x=TRUE)
write.table(common_all, "/home/famgarcia/Escritorio/intersec_all3subsets.tsv", sep="\t", row.names=F, quote=F)



##### Common for HDACi but not ME

common_HDACi = intersect(x[[1]], x[[2]])
common_HDACnoME = data.frame(external_gene_name = common_HDACi[!(common_HDACi %in% x[[3]])])
common_HDACnoME = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], common_HDACnoME, by="external_gene_name", all.y=T)
common_HDACnoME = merge(panther_ann, common_HDACnoME, by="external_gene_name", all.y=T)
common_HDACnoME = merge(ensembl_ann, common_HDACnoME, by="external_gene_name", all.y=T)
common_HDACnoME = merge(common_HDACnoME, expr_HDAC_Stress[,c(3,1)], by="external_gene_name", all.x=T)
common_HDACnoME = merge(common_HDACnoME, expr_HDAC[,c(3,1)], by="external_gene_name", all.x=T)

common_HDACnoME_GO = AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(common_HDACnoME$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
common_HDACnoME_GO = common_HDACnoME_GO[common_HDACnoME_GO$ONTOLOGY=="BP",]
common_HDACnoME_GO= common_HDACnoME_GO[!(common_HDACnoME_GO$GO == "NA"),]

xx <- as.list(common_HDACnoME_GO$GO)

# Get the TERMS for the first element of xx

list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
common_HDACnoME_GO$GOTERM = list_terms
common_HDACnoME_GO= common_HDACnoME_GO[!is.na(common_HDACnoME_GO$GO),]
common_HDACnoME = merge(common_HDACnoME_GO, common_HDACnoME, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
common_HDACnoME = common_HDACnoME[,!(colnames(common_HDACnoME) %in% ann_discard)]
common_HDACnoME$lfc_HeatStress_SAHA = paste("LFC HeatStress-SAHA : ", round(common_HDACnoME$lfc_HeatStress_SAHA,2))
common_HDACnoME$lfc_PE_SAHA = paste("LFC PE-SAHA : ", round(common_HDACnoME$lfc_PE_SAHA,2))
common_HDACnoME = merge(common_HDACnoME, VST_values, by="external_gene_name", all.x=TRUE)
write.table(common_HDACnoME, "/home/famgarcia/Escritorio/intersec_SAHA_treatment_outersec_HeatStress.tsv", sep="\t", row.names=F, quote=F)


####### Common in embryogenesis but not in differentially expressed genes from treated and untreated proembryos


common_HDACi = intersect(x[[2]], x[[3]])
common_MEnoHDAC = data.frame(external_gene_name = common_HDACi[!(common_HDACi %in% x[[1]])])
common_MEnoHDAC = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], common_MEnoHDAC, by="external_gene_name", all.y=T)
common_MEnoHDAC = merge(panther_ann, common_MEnoHDAC, by="external_gene_name", all.y=T)
common_MEnoHDAC = merge(ensembl_ann, common_MEnoHDAC, by="external_gene_name", all.y=T)
common_MEnoHDAC = merge(common_MEnoHDAC, expr_HDAC_Stress[,c(3,1)], by="external_gene_name", all.x=T)
common_MEnoHDAC = merge(common_MEnoHDAC, expr_Stress[,c(3,1)], by="external_gene_name", all.x=T)

common_MEnoHDAC_GO = AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(common_MEnoHDAC$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
common_MEnoHDAC_GO = common_MEnoHDAC_GO[common_MEnoHDAC_GO$ONTOLOGY=="BP",]
common_MEnoHDAC_GO= common_MEnoHDAC_GO[!(common_MEnoHDAC_GO$GO == "NA"),]

xx <- as.list(common_MEnoHDAC_GO$GO)

# Get the TERMS for the first element of xx

list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
common_MEnoHDAC_GO$GOTERM = list_terms
common_MEnoHDAC_GO= common_MEnoHDAC_GO[!is.na(common_MEnoHDAC_GO$GO),]
common_MEnoHDAC = merge(common_MEnoHDAC_GO, common_MEnoHDAC, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
common_MEnoHDAC = common_MEnoHDAC[,!(colnames(common_MEnoHDAC) %in% ann_discard)]
common_MEnoHDAC$lfc_HeatStress_SAHA = paste("LFC HeatStress-SAHA : ", round(common_MEnoHDAC$lfc_HeatStress_SAHA,2))
common_MEnoHDAC$lfc_HeatStress = paste("LFC HeatStress : ", round(common_MEnoHDAC$lfc_HeatStress,2))
common_MEnoHDAC = merge(common_MEnoHDAC, VST_values, by="external_gene_name", all.x=TRUE)
write.table(common_MEnoHDAC, "/home/famgarcia/Escritorio/intersec_HeatStress_Embryogenesis_outersec_PE_treatment.tsv", sep="\t", row.names=F, quote=F)



##############################################################################################################

# Biological question 2: What is the differernce in embryogenesis taking and not taking into account SAHA treatment?

################################################
# a) not segmentating by repression / activation
########################### MAKE VENN
x <- list(
  HeatStress_SAHA = expr_HDAC_Stress[, 3], 
  HeatStress = expr_Stress[ , 3]
)

png("/home/famgarcia/Escritorio/Intersection_DEG_Embryogenesis_notsegmented.png", width=16, height=12, units="cm", res=300)
par(mar=c(0,3,1.5,1))
ggVennDiagram(x, label_alpha = 0, edge_size = 0, set_size=3,
              set_color=c("HeatStress_SAHA" = "lightcoral", "HeatStress" = "gold")
) + scale_fill_gradient(low="white",high = "tan") + 
  theme(legend.position="none",) + # no legend
  scale_x_continuous(expand = expansion(mult = .2)) # avoid cropping long labels

dev.off()

############################################## Write annotated gene intersections in tables

common_ME = data.frame(external_gene_name = intersect(x[[1]],x[[2]]))
common_ME = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], common_ME, by="external_gene_name", all.y=T)
common_ME = merge(panther_ann, common_ME, by="external_gene_name", all.y=T)
common_ME = merge(ensembl_ann, common_ME, by="external_gene_name", all.y=T)
common_ME = merge(common_ME, expr_HDAC_Stress[,c(3,1)], by="external_gene_name", all.x=T)
common_ME = merge(common_ME, expr_Stress[,c(3,1)], by="external_gene_name", all.x=T)

columns(org.Brassicanapus.eg.db)
common_ME_GO = AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(common_ME$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
common_ME_GO = common_ME_GO[common_ME_GO$ONTOLOGY=="BP",]
common_ME_GO= common_ME_GO[!(common_ME_GO$GO == "NA"),]
xx <- as.list(common_ME_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
common_ME_GO$GOTERM = list_terms
common_ME_GO= common_ME_GO[!is.na(common_ME_GO$GO),]
common_ME = merge(common_ME_GO, common_ME, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
common_ME = common_ME[,!(colnames(common_ME) %in% ann_discard)]
common_ME$lfc_HeatStress_SAHA = paste("LFC HeatStress-SAHA : ", round(common_ME$lfc_HeatStress_SAHA,2))
common_ME$lfc_HeatStress = paste("LFC HeatStress : ", round(common_ME$lfc_HeatStress,2))
common_ME = merge(common_ME, VST_values, by="external_gene_name", all.x=TRUE)
write.table(common_ME, "/home/famgarcia/Escritorio/intersec_Embryogenesis_all.tsv", sep="\t", row.names=F, quote=F)

# save common and not common subsets

# present in SAHA treated embryogenesis but not in conventional embryogenesis
not_common_HDAC = expr_HDAC_Stress[!expr_HDAC_Stress$external_gene_name %in% intersect(x[[1]],x[[2]]),] 
not_common_HDAC = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], not_common_HDAC, by="external_gene_name", all.y=T)
not_common_HDAC = merge(panther_ann, not_common_HDAC, by="external_gene_name", all.y=T)
not_common_HDAC = merge(ensembl_ann, not_common_HDAC, by="external_gene_name", all.y=T)
not_common_HDAC = not_common_HDAC[order(not_common_HDAC$lfc_HeatStress_SAHA),]

columns(org.Brassicanapus.eg.db)
not_common_HDAC_GO = AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(not_common_HDAC$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
not_common_HDAC_GO = not_common_HDAC_GO[not_common_HDAC_GO$ONTOLOGY=="BP",]
not_common_HDAC_GO= not_common_HDAC_GO[!(not_common_HDAC_GO$GO == "NA"),]
xx <- as.list(not_common_HDAC_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
  } 
not_common_HDAC_GO$GOTERM = list_terms
not_common_HDAC_GO= not_common_HDAC_GO[!is.na(not_common_HDAC_GO$GO),]
not_common_HDAC = merge(not_common_HDAC_GO, not_common_HDAC, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
not_common_HDAC = not_common_HDAC[,!(colnames(not_common_HDAC) %in% ann_discard)]
not_common_HDAC$lfc_HeatStress_SAHA = paste("LFC : ", round(not_common_HDAC$lfc_HeatStress_SAHA,2))
not_common_HDAC = merge(not_common_HDAC, VST_values, by="external_gene_name", all.x=TRUE)

write.table(not_common_HDAC, "/home/famgarcia/Escritorio/outersec_Embryogenesis_in_HeatStress_SAHA.tsv", sep="\t", row.names=F, quote=F)

##################### NOT COMMON BUT PRESENT IN CONVENTIONAL EMBRYOGENESIS - 32ºC

not_common_HS = expr_Stress[!expr_Stress$external_gene_name %in% intersect(x[[1]],x[[2]]),] 
not_common_HS = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], not_common_HS, by="external_gene_name", all.y=T)
not_common_HS = merge(panther_ann, not_common_HS, by="external_gene_name", all.y=T)
not_common_HS = merge(ensembl_ann, not_common_HS, by="external_gene_name", all.y=T)
not_common_HS = not_common_HS[order(not_common_HS$lfc_HeatStress),]
not_common_HS_GO = AnnotationDbi::select(org.Brassicanapus.eg.db, as.character(not_common_HS$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
not_common_HS_GO = not_common_HS_GO[not_common_HS_GO$ONTOLOGY=="BP",]
not_common_HS_GO= not_common_HS_GO[!(not_common_HS_GO$GO == "NA"),]
xx <- as.list(not_common_HS_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 

not_common_HS_GO$GOTERM = list_terms
not_common_HS_GO= not_common_HS_GO[!is.na(not_common_HS_GO$GO),]
not_common_HS = merge(not_common_HS_GO, not_common_HS, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
not_common_HS$lfc_HeatStress = paste("LFC : ", round(not_common_HS$lfc_HeatStress,2))
not_common_HS = merge(not_common_HS, VST_values, by="external_gene_name", all.x=TRUE)

write.table(not_common_HS, "/home/famgarcia/Escritorio/outersec_Embryogenesis_in_HeatStress.tsv", sep="\t", row.names=F, quote=F)


# Biological question 2: What is the differernce in embryogenesis taking and not taking into account SAHA treatment?
#################################################
# b) segmentating by repression / activation
################################ MAKE VENN DEG
x <- list(
  "HeatStress-SAHA UP" = expr_HDAC_Stress[expr_HDAC_Stress$log2FoldChange> 0, ]$external_gene_name, 
  "HeatStress-SAHA DOWN"  = expr_HDAC_Stress[expr_HDAC_Stress$log2FoldChange < 0, ]$external_gene_name, 
  "HeatStress DOWN" = expr_Stress[(expr_Stress$log2FoldChange < 0) , ]$external_gene_name,
  "HeatStress UP" = expr_Stress[expr_Stress$log2FoldChange > 0 ,]$external_gene_name
)

svg("/home/famgarcia/Escritorio/Intersection_DEG_Embryogenesis_segmented.svg", width=10, height=7)
par(mar=c(0,3,1.5,1))
ggVennDiagram(x, label_alpha = 0, edge_size = 0.3, set_size=6,
              set_color=c("ME_HDACi_UP" = "lightcoral", "ME_HDACi_DOWN" = "green4", "ME_DOWN" = "steelblue", "ME_UP" = "purple")
) + scale_fill_gradient(low="white",high = "tan") + 
  theme(legend.position="none",) + # no legend
  scale_x_continuous(expand = expansion(mult = .2)) # avoid cropping long labels

dev.off()

write.table(data.frame(common_up=intersect(x[[1]], x[[4]])), "/home/famgarcia/Descargas/intersec_up.csv", row.names=F)
write.table(not_common_HDAC, "/home/famgarcia/Descargas/outersec_DEG_ME_HDAC.csv", row.names=F)

################### GET DEGs COMMON AMONG TWO CONDITIONS, BUT SEGMENTED BY different TYPE OF EXPRESSION

############################################## Write annotated gene intersections in tables

##### ME_HDACi DOWN, ME UP - ONLY NON-ZERO DIFFERENT BEHAVIOUR SUBSET

not_common_dExupME = expr_HDAC_Stress[expr_HDAC_Stress$external_gene_name %in% intersect(x[[2]],x[[4]]),] 
not_common_dExupME = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], not_common_dExupME, by="external_gene_name", all.y=T)
not_common_dExupME = merge(not_common_dExupME, expr_Stress, by="external_gene_name", all.x=T)
not_common_dExupME = merge(panther_ann, not_common_dExupME, by="external_gene_name", all.y=T)
not_common_dExupME = merge(ensembl_ann, not_common_dExupME, by="external_gene_name", all.y=T)
not_common_dExupME$lfc_HeatStress_SAHA = paste("LFC HeatStress-SAHA : ", round(not_common_dExupME$lfc_HeatStress_SAHA,2))
not_common_dExupME$lfc_HeatStress = paste("LFC HeatStress : ", round(not_common_dExupME$lfc_HeatStress,2))
not_common_dExupME = merge(not_common_dExupME, VST_values, by="external_gene_name", all.x=TRUE)
not_common_dExupME<-not_common_dExupME[,!(colnames(not_common_dExupME) %in% ann_discard)]

write.table(not_common_dExupME, "/home/famgarcia/Escritorio/intersec_Embryogenesis_DOWN_HeatStress_SAHA_UP_HeatStress_all.tsv", sep="\t", row.names=F, quote=F)
