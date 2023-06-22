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

# Making annotation package
ah <- AnnotationHub()
#Ensembl 56
query(ah, c( "OrgDb", "Brassica napus"))
org.Brassicanapus.eg.db <- ah[["AH107389"]]


################ Load and preprocess DEG data
# Load DEA data - necessary to filter lists and make coexpressed gene output findings more robust

expr_HDAC <- read.csv("/home/famgarcia/T_vs_C_DEA/Final_Condition_T_vs_C_LFC1_padj05.csv")
expr_HDAC_Stress <- read.csv("/home/famgarcia/Escritorio/Condition_T_vs_VM_result_padj_0.05_DevStage.csv")
expr_Stress <- read.csv("/home/famgarcia/Escritorio/Condition_PEall_vs_VM_result_padj_0.05_DevStage.csv")

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
  ME_HDACi = expr_HDAC_Stress[, 3], 
  ME = expr_Stress[ ,3], 
  HDAC = expr_HDAC[!is.na(expr_HDAC$external_gene_name), 3]
)

################# MAKE DEG Venn plot
png("/home/famgarcia/Escritorio/IntersectionDEGgenes_all.png", width=16, height=12, units="cm", res=300)
ggVennDiagram(x, c("ME_HDACi", "ME", "HDAC"), edge_size = 0.3,
              color=c("HDAC" = "lightblue", "ME_HDACi" = "lightcoral", "ME" = "khaki"),
              set_color = c("HDAC" = "steelblue", "ME_HDACi" = "lightcoral", "ME" = "khaki")) +
  scale_fill_gradient(low="white",high = "steelblue") + 
  scale_color_manual(values = c("HDAC" = "lightblue", "ME_HDACi" = "lightcoral", "ME" = "khaki")) +
  theme(legend.position="none")+ # no legend
  scale_x_continuous(expand = expansion(mult = .2)) # avoid cropping long labels
dev.off()


############################################## Write annotated gene intersections in tables
##### Common for HDACi but not ME
common_all = annot_ncbi[annot_ncbi$external_gene_name %in% intersect(intersect(x[[1]], x[[3]]), x[[2]]) , ]
common_all_GO = select(org.Brassicanapus.eg.db, as.character(common_all$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
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

##### Common for HDACi but not ME

common_HDACi = intersect(x[[1]], x[[3]])
common_HDACnoME = annot_ncbi[annot_ncbi$external_gene_name %in% common_HDACi[!(common_HDACi %in% xx[[2]])], ]
common_HDACnoME_GO = select(org.Brassicanapus.eg.db, as.character(common_HDACnoME$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
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






##############################################################################################################

# Biological question 2: What is the differernce in embryogenesis taking and not taking into account SAHA treatment?

################################################
# a) not segmentating by repression / activation
########################### MAKE VENN
x <- list(
  ME_HDACi = expr_HDAC_Stress[, 3], 
  ME = expr_Stress[ , 3]
)

png("/home/famgarcia/Escritorio/IntersectionDEGgenesME_notsegment.png", width=16, height=12, units="cm", res=300)
par(mar=c(0,3,1.5,1))
ggVennDiagram(x, label_alpha = 0, edge_size = 0,
              set_color=c("ME_HDACi" = "lightcoral", "ME" = "gold")
) + scale_fill_gradient(low="white",high = "tan") + 
  theme(legend.position="none",) + # no legend
  scale_x_continuous(expand = expansion(mult = .2)) # avoid cropping long labels

dev.off()

############################################## Write annotated gene intersections in tables
# save common and not common subsets
not_common_HDAC = expr_HDAC_Stress[!expr_HDAC_Stress$external_gene_name %in% intersect(x[[1]],x[[2]]),] 
not_common_HDAC = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], not_common_HDAC, by="external_gene_name")
not_common_HDAC = not_common_HDAC[order(not_common_HDAC$log2FoldChange),]

not_common_ME = expr_Stress[!expr_HDAC_Stress$external_gene_name %in% intersect(x[[1]],x[[2]]),] 
not_common_ME = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], not_common_HDAC, by="external_gene_name")
not_common_ME = not_common_HDAC[order(not_common_HDAC$log2FoldChange),]

columns(org.Brassicanapus.eg.db)
not_common_HDAC_GO = select(org.Brassicanapus.eg.db, as.character(not_common_HDAC$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
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



# Biological question 2: What is the differernce in embryogenesis taking and not taking into account SAHA treatment?
#################################################
# b) segmentating by repression / activation
################################ MAKE VENN DEG
x <- list(
  ME_HDACi_UP = expr_HDAC_Stress[expr_HDAC_Stress$log2FoldChange > 0, 3], 
  ME_HDACi_DOWN = expr_HDAC_Stress[expr_HDAC_Stress$log2FoldChange < 0, 3], 
  ME_DOWN = expr_Stress[(expr_Stress$log2FoldChange < 0) , 3],
  ME_UP = expr_Stress[expr_Stress$log2FoldChange > 0 ,3]
)

png("/home/famgarcia/Escritorio/IntersectionDEGgenes.png", width=16, height=12, units="cm", res=300)
par(mar=c(0,3,1.5,1))
ggVennDiagram(x, label_alpha = 0, edge_size = 0.3,
              set_color=c("ME_HDACi_UP" = "lightcoral", "ME_HDACi_DOWN" = "green4", "ME_DOWN" = "steelblue", "ME_UP" = "purple")
) + scale_fill_gradient(low="white",high = "tan") + 
  theme(legend.position="none",) + # no legend
  scale_x_continuous(expand = expansion(mult = .2)) # avoid cropping long labels

dev.off()

write.table(data.frame(common_up=intersect(x[[1]], x[[4]])), "/home/famgarcia/Descargas/intersec_up.csv", row.names=F)
write.table(not_common_HDAC, "/home/famgarcia/Descargas/outersec_DEG_ME_HDAC.csv", row.names=F)

################### GET DEGs COMMON AMONG TWO CONDITIONS, BUT SEGMENTED BY different TYPE OF EXPRESSION

############################################## Write annotated gene intersections in tables
##### ME_HDACi DOWN, ME_UP

not_common_dExdownME = expr_HDAC_Stress[expr_HDAC_Stress$external_gene_name %in% intersect(x[[1]],x[[3]]),] 
not_common_dExdownME = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], not_common_dExdownME, by="external_gene_name")
not_common_dExdownME = not_common_dExdownME[order(not_common_dExdownME$log2FoldChange),]

not_common_dExdownME_GO = select(org.Brassicanapus.eg.db, as.character(not_common_dExdownME$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
not_common_dExdownME_GO = not_common_dExdownME_GO[not_common_dExdownME_GO$ONTOLOGY=="BP",]
xx <- as.list(not_common_dExdownME_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
not_common_dExdownME_GO$GOTERM = list_terms
not_common_dExdownME_GO= not_common_dExdownME_GO[!is.na(not_common_dExdownME_GO$GO),]
not_common_dExdownME_GO = merge(not_common_dExdownME_GO, not_common_dExdownME, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
write.table(not_common_dExdownME_GO, "/home/famgarcia/Descargas/intersec_DEG_ME_HDAC_UP_ME_DOWN.csv", row.names=F)


##### ME_HDACi UP, ME_DOWN

not_common_dExupME = expr_HDAC_Stress[expr_HDAC_Stress$external_gene_name %in% intersect(x[[2]],x[[4]]),] 
not_common_dExupME = merge(annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], not_common_dExupME, by="external_gene_name")
not_common_dExupME = not_common_dExupME[order(not_common_dExupME$log2FoldChange),]

not_common_dExupME_GO = select(org.Brassicanapus.eg.db, as.character(not_common_dExupME$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
not_common_dExupME_GO = not_common_dExupME_GO[not_common_dExupME_GO$ONTOLOGY=="BP",]
xx <- as.list(not_common_dExupME_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
not_common_dExupME_GO$GOTERM = list_terms
not_common_dExupME_GO= not_common_dExupME_GO[!is.na(not_common_dExupME_GO$GO),]
not_common_dExupME_GO = merge(not_common_dExupME_GO, not_common_dExupME, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
write.table(not_common_dExupME_GO, "/home/famgarcia/Descargas/intersec_DEG_ME_HDAC_DOWN_ME_UP.csv", row.names=F)

write.table(common_all, "/home/famgarcia/Descargas/intersec_all.csv", row.names=F)
