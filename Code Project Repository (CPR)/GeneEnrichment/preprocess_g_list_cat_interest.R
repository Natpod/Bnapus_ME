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
#logFC_Stress<-read_csv("C:\\Users\\naata\\MASTER\\Enrichment\\Condition_PEall_vs_VM_result_padj_0.05_DevStage.csv")
logFC_Stress<-read.csv("C:\\Users\\naata\\MASTER\\Enrichment\\Final_condition_PErep_vs_VMrep_padj0.05_woPErep3_LFC1filtered.csv")



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








####### ######### ######## RADAR/BARPLOTS



hormones_merged = merge(gsea_Stress_H[,c(1,5,7)],gsea_HDAC_Stress_H[,c(1,5,7)], by="ID", all.x=TRUE, all.y=TRUE)
hormones_merged_2 = merge(hormones_merged, gsea_HDAC_H[,c(1,5,7)], by="ID", all.x=TRUE, all.y=TRUE)

hormones_merged_2$ID
tmppt = merge(data.frame(ID=hormones_merged_2$ID), pathways_BP_terms, by.x="ID", by.y="bp", all.x=T)


# Add merged description and setSize
hormones_merged_2$Description = tmppt$name



################## General non-auxin hormones
data_h_general = subset(hormones_merged_2, ( grepl("biosynthetic", Description) | grepl("auxin", Description) | grepl("signaling", Description) | grepl("biosyn", Description)) & !(grepl("regulation of", Description)))

# get data into a different format ( )
data_h_general = data_h_general[ ,colnames(data_h_general) %in%c("ID","Description","NES.x","NES.y","NES") ]

# reformat description to include GO code and term name with linebreaks
data_h_general$Description = gsub("to ", "to\n", data_h_general$Description)
data_h_general$Description = gsub(" signaling", "\nsignaling", data_h_general$Description)
data_h_general$Description = gsub(" biosynthetic", "\nbiosynthetic", data_h_general$Description)
data_h_general$Description = gsub(" the", "\nthe", data_h_general$Description)
data_h_general$Description = paste(data_h_general$ID, data_h_general$Description, sep = "\n")
data_h_general = data_h_general[ ,colnames(data_h_general) %in%c("Description","NES.x","NES.y", "NES") ]

# transpose
trans_data_h = t(data_h_general)
colnames(trans_data_h) <- trans_data_h[ rownames(trans_data_h) %in% c("Description"), ]
trans_data_h <- trans_data_h[1:3,]
trans_data_h <- as.data.frame(trans_data_h)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
max_NES = max(c( max(data_h_general[!is.na(data_h_general[,1]),1]), max(data_h_general[!is.na(data_h_general[,2]),2]) ) )
min_NES = min(c( min(data_h_general[!is.na(data_h_general[,1]),1]), min(data_h_general[!is.na(data_h_general[,2]),2]) ) )

# plot with default options:
#radarchart(as.data.frame(trans_data_h))

trans_data_h[is.na(trans_data_h)] <- 0
trans_data_h[,] <- as.matrix(trans_data_h[,])
trans_data_h<-as.data.frame(sapply(trans_data_h, as.numeric))
rownames(trans_data_h) <- c("HeatStress", "HeatStress-SAHA", "PE-SAHA")

trans_data_h <- trans_data_h %>% rownames_to_column("group")

svg("C:\\Users\\naata\\Desktop\\Radar_hormone_allterms.svg", width=7, height=7)
par(mar=c(8,8,9,8))
p2<-ggradar(
  trans_data_h, 
  axis.label.size = 3,
  axis.label.offset = 1.45,
  values.radar = c("-2", "0", "NES\n2"),
  grid.min = -2, grid.mid = 0, grid.max = 2.05,
  grid.label.size = 3.5,
  gridline.label.offset = -0.25,
  # Relative plot size
  plot.extent.x.sf = 2,
  plot.extent.y.sf = 1.59,
  # Polygons
  group.line.width = 1.1, 
  group.point.size = 2.4,
  group.colours = c("#E1C16E", "lightcoral", "#d3d3d3"),
  # Background and grid lines
  gridline.mid.linetype = "dashed",
  gridline.min.linetype = "solid",
  gridline.max.linetype = "solid",
  grid.line.width = 0.5,
  axis.line.colour ="#DCDCDC",
  gridline.max.colour = "#a32a31",
  gridline.min.colour = "green3",
  gridline.mid.colour = "lightblue4",
  background.circle.colour = "#fbf9f4",
  legend.position = "bottom",
  legend.text.size = "8"
  
) + theme(plot.background = element_rect(fill='transparent', color=NA))

dev.off()



hormones_merged_barplot = hormones_merged_2[,]
hormones_merged_barplot$ID = paste(hormones_merged_barplot$ID,hormones_merged_barplot$Description, sep = " -- ")
hormones_merged_barplot = hormones_merged_barplot[,colnames(hormones_merged_barplot) %in% c("ID","NES.x", "NES.y", "NES")]
hormones_merged_barplot[is.na(hormones_merged_barplot)]<-0
colnames(hormones_merged_barplot)<-c("ID","HeatStress","HeatStress_SAHA", "PE_SAHA")

df = melt(hormones_merged_barplot)
colnames(df)<- c("ID", "Condition", "NES")
order_used = order( unlist(lapply(df$ID, function(x){ unlist(str_split(x, "--"))[2] } )) )


levels_order = df[order_used,1]


png("C:\\Users\\naata\\Downloads\\hormone_terms_Stress_vs_Stress_SAHA.png", width=23, height=10, units="cm", res=300)
par(mar=c(3,3,5,20))
ggplot(data = df, mapping = aes(y=NES,x=ID, fill=Condition))+
  geom_bar(aes(x = as.factor(ID), group=Condition),stat = "identity", position = 'dodge', width = 1.2) + 
  scale_x_discrete(limits = levels_order) +
  scale_fill_manual("Condition", values = c("HeatStress" = "#E1C16E" , "HeatStress_SAHA" = "lightcoral", "PE_SAHA" = "#d3d3d3")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 7), 
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_blank(),
        plot.margin = margin(0.3,0,0,4, "cm"),
        legend.text = element_text(size = 8))
dev.off()





############### ep


epi_merged = merge(gsea_Stress_EPI[,c(1,5,7)],gsea_HDAC_Stress_EPI[,c(1,5,7)], by="ID", all.x=TRUE, all.y=TRUE)
epi_merged_2 = merge(epi_merged, gsea_HDAC_EPI[,c(1,5,7)], by="ID", all.x=TRUE, all.y=TRUE)

tmppt = merge(data.frame(ID=epi_merged_2$ID), pathways_BP_terms, by.x="ID", by.y="bp", all.x=T)



# Add merged description and setSize
epi_merged_2$Description = tmppt$name


epi_merged_barplot = epi_merged_2[,]
epi_merged_barplot$ID = paste(epi_merged_barplot$ID,epi_merged_barplot$Description, sep = " -- ")
epi_merged_barplot = epi_merged_barplot[,colnames(epi_merged_barplot) %in% c("ID","NES.x", "NES.y", "NES")]
epi_merged_barplot[is.na(epi_merged_barplot)]<-0
colnames(epi_merged_barplot)<-c("ID","HeatStress","HeatStress_SAHA", "PE_SAHA")

df = melt(epi_merged_barplot)
colnames(df)<- c("ID", "Condition", "NES")


#order_used = order( unlist(lapply(df$ID, function(x){ unlist(str_split(x, "--"))[2] } )) )
order_used = c( which(grepl("epigenetic", df$ID)),which(grepl("chromatin", df$ID)), which(grepl("heterochromatin", df$ID)),which(grepl("methya", df$ID)), which(grepl("nucleosome", df$ID)), which(grepl("histone", df$ID)))
   

levels_order = df[order_used,1]


svg("C:\\Users\\naata\\Desktop\\epi_terms_Stress_vs_Stress_SAHA.svg", width=23, height=6, bg="transparent" )
par(mar=c(3,3,5,20))
p1<- ggplot(data = df, mapping = aes(y=NES,x=ID, fill=Condition))+
  geom_bar(aes(x = as.factor(ID), group=Condition),stat = "identity", position = 'dodge', width = 1.2) + 
  scale_x_discrete(limits = levels_order) +
  scale_fill_manual("Condition", values = c("HeatStress" = "#E1C16E" , "HeatStress_SAHA" = "lightcoral", "PE_SAHA" = "grey")) + theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 9), 
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(0.3,0,0,4, "cm"),
        legend.text = element_text(size = 8))
dev.off()




ggsave(filename="C:\\Users\\naata\\Desktop\\epi_terms_Stress_vs_Stress_SAHA.svg", plot = p1, width=23, height=6, bg='transparent')

ggsave(filename="C:\\Users\\naata\\Desktop\\radarplot_transparent.svg", plot = p2, width=7, height=7, bg='transparent')


































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


g_bA_summary = data.frame(table(g_bA$GO))
g_bB_summary = data.frame(table(g_bB$GO))

summary_WGCNA = merge(g_bA_summary, gB_summary, by="Var1", all.x=T, all.y=T)
colnames(summary_WGCNA) <- c("bp","Coexpression_Proembryo", "Coexpression_VM")
summary_WGCNA = merge(summary_WGCNA, pathways_BP_terms, by="bp", all.x=T)
summary_WGCNA<-summary_WGCNA[,c(1,4,2,3)]

summary_WGCNA[which(is.na(summary_WGCNA[,3])),3]<-0
summary_WGCNA[which(is.na(summary_WGCNA[,4])),4]<-0
summary_WGCNA$Total <- rowSums(summary_WGCNA[,c(3,4)])



write.table(summary_WGCNA, "C:\\Users\\naata\\Documents\\MASTER\\TFM\\tfm\\WGCNA_behaviours_summary.tsv", sep="\t", row.names=F, quote=F)




###########################
####### Separate by expression - WGCNA

# gA

# rank HeatStress - get df with + / - genes

# rank HeatStress-SAHA


###### Separate by expression  - DEG subsets.

merge(genes_epi_HDAC_Stress, genes_epi_Stress, by="external_gene_name", all.x=T, all.y=T)



# Saving term frequency statistics in file


############ Proembrión

PE_SAHA_DEGs = logFC_HDAC[,c(9,3)]
PE_SAHA_DEGs = PE_SAHA_DEGs[!duplicated(PE_SAHA_DEGs$external_gene_name),]

PE_SAHA_DEGs = merge(PE_SAHA_DEGs, PANTHER_annot, by="external_gene_name")
PE_SAHA_DEGs_h = merge(PE_SAHA_DEGs,df_terms_h,by.x="external_gene_name", by.y="gene")
PE_SAHA_DEGs_h = merge(PE_SAHA_DEGs_h, hormone_terms[,1:3], by = "ID")
PE_SAHA_DEGs_h = merge(PE_SAHA_DEGs_h, VST_values, by="external_gene_name")

PE_SAHA_DEGs_epi = merge(PE_SAHA_DEGs,df_terms_ep,by.x="external_gene_name", by.y="gene")
PE_SAHA_DEGs_epi = merge(PE_SAHA_DEGs_epi, epi_terms[,1:3], by = "ID")
PE_SAHA_DEGs_epi = merge(PE_SAHA_DEGs_epi, VST_values, by="external_gene_name")

epiterms_sum_PE = merge(data.frame(table(PE_SAHA_DEGs_epi[grepl("-",PE_SAHA_DEGs_epi$log2FoldChange),]$context)), data.frame(table(PE_SAHA_DEGs_epi[!(grepl("-",PE_SAHA_DEGs_epi$log2FoldChange)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(epiterms_sum_PE) <- c("context", "down", "up")

hterms_sum_PE = merge(data.frame(table(PE_SAHA_DEGs_h[grepl("-",PE_SAHA_DEGs_h$log2FoldChange),]$context)), data.frame(table(PE_SAHA_DEGs_h[!(grepl("-",PE_SAHA_DEGs_h$log2FoldChange)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(hterms_sum_PE) <- c("context", "down", "up")








#######################################

HS_DEGs = logFC_Stress[,c(11,3)]
HS_DEGs = HS_DEGs[!duplicated(HS_DEGs$external_gene_name),]

HS_DEGs = merge(HS_DEGs, PANTHER_annot, by="external_gene_name")

HS_DEGs_h = merge(HS_DEGs,df_terms_h,by.x="external_gene_name", by.y="gene")
HS_DEGs_h = merge(HS_DEGs_h, hormone_terms[,1:3], by = "ID")
HS_DEGs_h = merge(HS_DEGs_h, VST_values, by="external_gene_name")

HS_DEGs_epi = merge(HS_DEGs,df_terms_ep,by.x="external_gene_name", by.y="gene")
HS_DEGs_epi = merge(HS_DEGs_epi, epi_terms[,1:3], by = "ID")
HS_DEGs_epi = merge(HS_DEGs_epi, VST_values, by="external_gene_name")



epiterms_sum_HS = merge(data.frame(table(HS_DEGs_epi[grepl("-",HS_DEGs_epi$log2FoldChange),]$context)), data.frame(table(HS_DEGs_epi[!(grepl("-",HS_DEGs_epi$log2FoldChange)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(epiterms_sum_HS) <- c("context", "down", "up")

hterms_sum_HS = merge(data.frame(table(HS_DEGs_h[grepl("-",HS_DEGs_h$log2FoldChange),]$context)), data.frame(table(HS_DEGs_h[!(grepl("-",HS_DEGs_h$log2FoldChange)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(hterms_sum_HS) <- c("context", "down", "up")






SAHA_DEGs = logFC_HDAC_Stress[,c(9,4)]
SAHA_DEGs = SAHA_DEGs[!duplicated(SAHA_DEGs$external_gene_name),]

SAHA_DEGs = merge(SAHA_DEGs, PANTHER_annot, by="external_gene_name")

SAHA_DEGs_h = merge(SAHA_DEGs,df_terms_h,by.x="external_gene_name", by.y="gene")
SAHA_DEGs_h = merge(SAHA_DEGs_h, hormone_terms[,1:3], by = "ID")
SAHA_DEGs_h = merge(SAHA_DEGs_h, VST_values, by="external_gene_name")

SAHA_DEGs_epi = merge(SAHA_DEGs,df_terms_ep,by.x="external_gene_name", by.y="gene")
SAHA_DEGs_epi = merge(SAHA_DEGs_epi, epi_terms[,1:3], by = "ID")
SAHA_DEGs_epi = merge(SAHA_DEGs_epi, VST_values, by="external_gene_name")



epiterms_sum_SAHA = merge(data.frame(table(SAHA_DEGs_epi[grepl("-",SAHA_DEGs_epi$log2FoldChange),]$context)), data.frame(table(SAHA_DEGs_epi[!(grepl("-",SAHA_DEGs_epi$log2FoldChange)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(epiterms_sum_SAHA) <- c("context", "down", "up")

hterms_sum_SAHA = merge(data.frame(table(SAHA_DEGs_h[grepl("-",SAHA_DEGs_h$log2FoldChange),]$context)), data.frame(table(SAHA_DEGs_h[!(grepl("-",SAHA_DEGs_h$log2FoldChange)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(hterms_sum_SAHA) <- c("context", "down", "up")





#######################################


HS_DEGS_all = rbind(DEGS_hormone_HS, DEGS_hormone_common[, colnames(DEGS_hormone_common) != "lfc_HeatStress_SAHA"])
HS_DEGS_all$sign = ifelse(grepl("-",HS_DEGS_all$lfc_HeatStress), "downregulated","upregulated")




SAHA_DEGS_all = rbind(DEGS_hormone_SAHA, DEGS_hormone_common[, colnames(DEGS_hormone_common) != "lfc_HeatStress"])
SAHA_DEGS_all$sign = ifelse(grepl("-",SAHA_DEGS_all$lfc_HeatStress_SAHA), "downregulated","upregulated")



hterms_sum_HS = merge(data.frame(table(HS_DEGS_all[grepl("-",HS_DEGS_all$lfc_HeatStress),]$context)), data.frame(table(HS_DEGS_all[!(grepl("-",HS_DEGS_all$lfc_HeatStress)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(hterms_sum_HS) <- c("context", "down", "up")


hterms_sum_SAHA = merge(data.frame(table(SAHA_DEGS_all[grepl("-",SAHA_DEGS_all$lfc_HeatStress_SAHA),]$context)), data.frame(table(SAHA_DEGS_all[!(grepl("-",SAHA_DEGS_all$lfc_HeatStress_SAHA)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(hterms_sum_SAHA) <- c("context", "down", "up")


HS_DEGS_all_epi = rbind(DEGS_epi_HS, DEGS_epi_common[, colnames(DEGS_hormone_common) != "lfc_HeatStress_SAHA"])
HS_DEGS_all_epi$sign = ifelse(grepl("-",HS_DEGS_all_epi$lfc_HeatStress), "downregulated","upregulated")

################################3 +


SAHA_DEGS_all_epi = rbind(DEGS_epi_SAHA, DEGS_epi_common[, colnames(DEGS_hormone_common) != "lfc_HeatStress"])
SAHA_DEGS_all_epi$sign = ifelse(grepl("-",SAHA_DEGS_all_epi$lfc_HeatStress_SAHA), "downregulated","upregulated")



epiterms_sum_HS = merge(data.frame(table(HS_DEGS_all_epi[grepl("-",HS_DEGS_all_epi$lfc_HeatStress),]$context)), data.frame(table(HS_DEGS_all_epi[!(grepl("-",HS_DEGS_all_epi$lfc_HeatStress)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(epiterms_sum_HS) <- c("context", "down", "up")


epiterms_sum_SAHA = merge(data.frame(table(SAHA_DEGS_all_epi[grepl("-",SAHA_DEGS_all_epi$lfc_HeatStress_SAHA),]$context)), data.frame(table(SAHA_DEGS_all_epi[!(grepl("-",SAHA_DEGS_all_epi$lfc_HeatStress_SAHA)),]$context)), by="Var1", all.x=T, all.y=T)
colnames(epiterms_sum_SAHA) <- c("context", "down", "up")

################################## PLOT

epiterms_sum_PE =melt(epiterms_sum_PE)
epiterms_sum_PE$Condition = "PE-SAHA"

hterms_sum_PE =melt(hterms_sum_PE)
hterms_sum_PE$Condition = "PE-SAHA"


epiterms_sum_HS =melt(epiterms_sum_HS)
epiterms_sum_HS$Condition = "HeatStress"

hterms_sum_HS =melt(hterms_sum_HS)
hterms_sum_HS$Condition = "HeatStress"


epiterms_sum_SAHA =melt(epiterms_sum_SAHA)
epiterms_sum_SAHA$Condition = "HeatStress-SAHA"

hterms_sum_SAHA =melt(hterms_sum_SAHA)
hterms_sum_SAHA$Condition = "HeatStress-SAHA"

df_plot_terms = rbind(epiterms_sum_PE, epiterms_sum_HS, epiterms_sum_SAHA)


df_plot_terms$context = gsub(" ", "\n", df_plot_terms$context)
df_plot_terms$variable = gsub("down", "Downregulated DEGs", df_plot_terms$variable)
df_plot_terms$variable = gsub("up", "Upregulated DEGs", df_plot_terms$variable)

svg("C:\\Users\\naata\\Downloads\\barplot_DEG_comparison_epig_Stressncbe.svg", width=13, height=4)
ggplot(df_plot_terms, aes(x=Condition, y=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001),
                                                                                                                                    legend.title = element_blank(), legend.position = "bottom" , axis.title.x=element_blank()) +
                                                                                                                      facet_grid(~context) + ylab("DEG gene count") + 
                                                                                                                        scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
 theme(strip.text.y = element_text(size = 17, angle=0),strip.background = element_rect(fill="#e3eeec", linewidth=NA),
       strip.placement = "outside")
dev.off()




df_plot_terms = rbind(hterms_sum_PE, hterms_sum_HS, hterms_sum_SAHA)
df_plot_terms$context = gsub(" ", "\n", df_plot_terms$context)
df_plot_terms$variable = gsub("down", "Downregulated DEGs", df_plot_terms$variable)
df_plot_terms$variable = gsub("up", "Upregulated DEGs", df_plot_terms$variable)



svg("C:\\Users\\naata\\Downloads\\barplot_DEG_comparison_hormone_Stressnbce.svg", width=10, height=4)
ggplot(df_plot_terms, aes(x=Condition, y=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "bottom", axis.title.x=element_blank()) +
  facet_grid(~context) + ylab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 17, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")

dev.off()



#########################3 

########################## INTERSECTION
HS_DEGs_h_common = merge(HS_DEGs_h, logFC_HDAC_Stress[,c(9,4)], by="external_gene_name")
HS_DEGs_epi_common = merge(HS_DEGs_epi, logFC_HDAC_Stress[,c(9,4)], by="external_gene_name")

hterm_DEGS_common_UP = data.frame(table(HS_DEGs_h_common[!(grepl("-" ,HS_DEGs_h_common$log2FoldChange.x)),]$term))
hterm_DEGS_common_UP$variable = rep("Upregulated DEGs",length(hterm_DEGS_common_UP[,1]))

hterm_DEGS_common_DOWN = data.frame(table(HS_DEGs_h_common[(grepl("-" ,HS_DEGs_h_common$log2FoldChange.x)),]$term))
hterm_DEGS_common_DOWN$variable = rep("Downregulated DEGs",length(hterm_DEGS_common_DOWN[,1]))


hterm_DEGS_common_UP = merge(hterm_DEGS_common_UP, hormone_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T, all.y=F)
hterm_DEGS_common_DOWN = merge(hterm_DEGS_common_DOWN, hormone_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T, all.y=F)


df_plot_terms = rbind(hterm_DEGS_common_UP, hterm_DEGS_common_DOWN)

ggplot(df_plot_terms, aes(y=Var1, x=Freq, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(context)) + xlab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.x = element_text(size = 17, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")


#write.table(df_plot_terms, "C:\\Users\\naata\\Downloads\\freqdegs.tsv", sep="\t", row.names = F, quote=F)
df_plot_terms_ann = read.csv( "C:\\Users\\naata\\Downloads\\freqdegs.tsv", sep="\t")
df_plot_terms_ann$context = gsub(" ", "\n", df_plot_terms_ann$context)

svg("C:\\Users\\naata\\Downloads\\barplot_DEG_intersectionME_hormone_Stressnbce.svg", width=5, height=6.5)
ggplot(df_plot_terms_ann, aes(y=type, x=count, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(context)) + xlab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")
dev.off()

################# EPIGENETIC

HS_DEGs_h_common = merge(HS_DEGs_epi, logFC_HDAC_Stress[,c(9,4)], by="external_gene_name")
HS_DEGs_epi_common = merge(HS_DEGs_epi, logFC_HDAC_Stress[,c(9,4)], by="external_gene_name")

epiterm_DEGS_common_UP = data.frame(table(HS_DEGs_epi_common[!(grepl("-" ,HS_DEGs_epi_common$log2FoldChange.x)),]$term))
epiterm_DEGS_common_UP$variable = rep("Upregulated DEGs",length(epiterm_DEGS_common_UP[,1]))

epiterm_DEGS_common_DOWN = data.frame(table(HS_DEGs_epi_common[(grepl("-" ,HS_DEGs_epi_common$log2FoldChange.x)),]$term))
epiterm_DEGS_common_DOWN$variable = rep("Downregulated DEGs",length(epiterm_DEGS_common_DOWN[,1]))

epiterm_DEGS_common_UP = merge(epiterm_DEGS_common_UP, epi_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T, all.y=F)
epiterm_DEGS_common_DOWN = merge(epiterm_DEGS_common_DOWN, epi_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T, all.y=F)


df_plot_terms = rbind(epiterm_DEGS_common_UP, epiterm_DEGS_common_DOWN)

#write.table(df_plot_terms, "C:\\Users\\naata\\Downloads\\freqdegs_EPI.tsv", sep="\t", row.names = F, quote=F)
df_plot_terms_ann = read.csv( "C:\\Users\\naata\\Downloads\\freqdegs_EPI.tsv", sep="\t")
df_plot_terms_ann$context = gsub(" ", "\n", df_plot_terms_ann$context)


p1 <- ggplot(df_plot_terms_ann[df_plot_terms_ann$context == "chromatin\nremodelling",], aes(y=reorder(term, count), x=count, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(context)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside")
p2 <- ggplot(df_plot_terms_ann[df_plot_terms_ann$context == "Histone\nmodification",], aes(y=reorder(term, count), x=count, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(context)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside")
p3 <- ggplot(df_plot_terms_ann[df_plot_terms_ann$context == "DNA\nmethylation",], aes(y=reorder(term, count), x=count, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(context)) + xlab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside")

svg("C:\\Users\\naata\\Downloads\\barplot_DEG_intersectionME_epig_Stressnbce.svg", width=5, height=5.2)
ggarrange(p1,p2,p3,nrow=3, ncol=1, common.legend = T, legend = "bottom", align="hv", heights = c(0.25,0.4,0.27))
dev.off()




######################### WGCNA
# positive in VM
g_bA_hormone = merge(g_bA, hormone_terms, by.x="ID.y", by.y="term")
g_bA_hormone = g_bA_hormone[!(grepl("NA", g_bA_hormone$log2FoldChange.x)),]

g_bA_epi = merge(g_bA, epi_terms, by.x="ID.x", by.y="term")
g_bA_epi = g_bA_epi[!(grepl("NA", g_bA_epi$log2FoldChange.x)),]

# positive in PE
g_bB_hormone = merge(g_bB, hormone_terms, by.x="ID.y", by.y="term")
g_bB_hormone = g_bB_hormone[!(grepl("NA", g_bB_hormone$log2FoldChange.x)),]


g_bB_epi = merge(g_bB, epi_terms, by.x="ID.x", by.y="term")
g_bB_epi = g_bB_epi[!(grepl("NA", g_bB_epi$log2FoldChange.x)),]


#a_epi_int = merge(g_bA_epi, HS_DEGs_epi_common, by="external_gene_name")
VMstagesum_epi_up = data.frame(table(g_bA_epi[!(grepl(" -", g_bA_epi$log2FoldChange.x)),]$ID.x))
VMstagesum_epi_down = data.frame(table(g_bA_epi[(grepl(" -", g_bA_epi$log2FoldChange.x)),]$ID.x))

#b_epi_int = merge(g_bB_epi, HS_DEGs_epi_common, by="external_gene_name") 
PEstagesum_epi_up = data.frame(table(g_bB_epi[!(grepl(" -", g_bB_epi$log2FoldChange.x)),]$ID.x))
PEstagesum_epi_down = data.frame(table(g_bB_epi[(grepl(" -", g_bB_epi$log2FoldChange.x)),]$ID.x))


VMstagesum_h = data.frame(table(g_bA_hormone$ID.y))
VMstagesum_h_up = data.frame(table(g_bA_hormone[!(grepl(" -", g_bA_hormone$log2FoldChange.x)),]$ID.y))
VMstagesum_h_down = data.frame(table(g_bA_hormone[(grepl(" -", g_bA_hormone$log2FoldChange.x)),]$ID.y))


PEstagesum_h_up = data.frame(table(g_bB_hormone[!(grepl(" -", g_bB_hormone$log2FoldChange.x)),]$ID.y))
PEstagesum_h_down = data.frame(table(g_bB_hormone[(grepl(" -", g_bB_hormone$log2FoldChange.x)),]$ID.y))

# PEstagesum_h = merge(PEstagesum_h_up, PEstagesum_h_down, by="Var1", all.x=T, all.y=T)
# VMstagesum_h = merge(VMstagesum_h_up, VMstagesum_h_down, by="Var1", all.x=T, all.y=T)
PEstagesum_h =PEstagesum_h_down
colnames(PEstagesum_h)<-c("terms", "Downregulated DEGs")
PEstagesum_h = merge(PEstagesum_h, hormone_terms[,c(2,3)], by.x="terms", by.y="term")

VMstagesum_h = VMstagesum_h_up
colnames(VMstagesum_h)<-c("terms", "Upregulated DEGs")

VMstagesum_h = merge(VMstagesum_h, hormone_terms[,c(2,3)], by.x="terms", by.y="term")

PEstagesum_epi = PEstagesum_epi_down
colnames(PEstagesum_epi)<-c("terms", "Downregulated DEGs")

PEstagesum_epi = merge(PEstagesum_epi, epi_terms[,c(2,3)], by.x="terms", by.y="term")
VMstagesum_epi = VMstagesum_epi_up
colnames(VMstagesum_epi)<-c("terms", "Upregulated DEGs")
VMstagesum_epi = merge(VMstagesum_epi, epi_terms[,c(2,3)], by.x="terms", by.y="term")

# VMstagesum_epi = merge(VMstagesum_epi_up, VMstagesum_epi_down, by="Var1", all.x=T, all.y=T)
# PEstagesum_epi = merge(PEstagesum_epi_up, PEstagesum_epi_down, by="Var1", all.x=T, all.y=T)


PEstagesum_h <- melt(PEstagesum_h)

VMstagesum_h <- melt(VMstagesum_h)

VMstagesum_epi <- melt(VMstagesum_epi)

PEstagesum_epi <- melt(PEstagesum_epi)


VMstagesum_h<-VMstagesum_h[,c(1,2,4,5)]
colnames(VMstagesum_h)<-c("terms","context", "variable", "value")


PEstagesum_h$Modules = rep("Coexpressed in PE", length(PEstagesum_h[,1]))
VMstagesum_h$Modules = rep("Coexpressed in VM", length(VMstagesum_h[,1]))
VMstagesum_epi$Modules = rep("Coexpressed in VM", length(VMstagesum_epi[,1]))
PEstagesum_epi$Modules = rep("Coexpressed in PE", length(PEstagesum_epi[,1]))

PEstagesum_h = PEstagesum_h[, c(1,2,3,6,4)]
df_to_plot = rbind(PEstagesum_h, VMstagesum_h)
df_to_plot$Modules = gsub(" ", "\n", df_to_plot$Modules )
df_to_plot$Var

p1<- ggplot(df_to_plot, aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=7.5, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "left", axis.title.y=element_blank()) +
  facet_grid(rows = vars(Modules)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_blank(),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside")



#df_to_plot = rbind(PEstagesum_epi, VMstagesum_epi)
#df_to_plot$Modules  = gsub(" ", "\n", df_to_plot$Modules )

df_to_plot$context = gsub(" ", "\n", df_to_plot$context)

p1<- ggplot(df_to_plot[df_to_plot$context == "auxin\nrelated",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(context ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 

p2<- ggplot(df_to_plot[df_to_plot$context == "abscisic\nacid\nrelated",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(context ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 
p3<- ggplot(df_to_plot[df_to_plot$context == "cytokinin\nrelated",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(context ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 
p4<- ggplot(df_to_plot[df_to_plot$context == "ethylene\nrelated",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(context ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 
p5<- ggplot(df_to_plot[df_to_plot$context == "gibberelin\nrelated",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(context ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 
p6<- ggplot(df_to_plot[df_to_plot$context == "JA\nrelated",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(context ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 
p7<- ggplot(df_to_plot[df_to_plot$context == "SA\nrelated",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(context ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 


svg("C:\\Users\\naata\\Downloads\\auxinwgcna.svg", width=7, height=1.6)
p1
dev.off()


svg("C:\\Users\\naata\\Downloads\\auxinwgcna.svg", width=13, height=7)
p
dev.off()

svg("C:\\Users\\naata\\Downloads\\abawgcna.svg", width=7, height=1.6)
p2
dev.off()

svg("C:\\Users\\naata\\Downloads\\CKwgcna.svg", width=7, height=1.2)
p3
dev.off()

svg("C:\\Users\\naata\\Downloads\\etwgcna.svg", width=5.5, height=1.2)
p4
dev.off()

svg("C:\\Users\\naata\\Downloads\\gibwgcna.svg", width=7, height=1.6)
p5
dev.off()

svg("C:\\Users\\naata\\Downloads\\JAwgcna.svg", width=7, height=1.6)
p6
dev.off()

svg("C:\\Users\\naata\\Downloads\\SAwgcna.svg", width=5.5, height=1.6)
p7
dev.off()

PEstagesum_epi = merge(PEstagesum_epi, epi_terms[,c(2,3)], by.x="terms", by.y="term")
PEstagesum_epi = PEstagesum_epi[,c(1,5,2,3,4)]
df_to_plot = rbind(PEstagesum_epi, VMstagesum_epi)
df_to_plot$Modules  = gsub(" ", "\n", df_to_plot$Modules )
#write.table(df_to_plot, "C:\\Users\\naata\\Downloads\\freqs_wgcna.tsv", sep="\t", row.names=F)
df_to_plot = read.table("C:\\Users\\naata\\Downloads\\freqs_wgcna.tsv", sep="\t", header=T)
df_to_plot$count = gsub(" ", "\n", df_to_plot$count)

p1<- ggplot(df_to_plot[df_to_plot$count == "chromatin\nremodelling",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(count ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 

p2<- ggplot(df_to_plot[df_to_plot$count == "histone\nmodification",], aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(count ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 
p3<- ggplot(df_to_plot[df_to_plot$count == "DNA\nmethylation",], aes(y=reorder(terms, value), x=value, fill=variable))+xlab("DEG count") + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=8, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(count ~Modules) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside") 





svg("C:\\Users\\naata\\Downloads\\WGCNA_epi_histone.svg", width=7, height=1.9)
p2
dev.off()


svg("C:\\Users\\naata\\Downloads\\WGCNA_epi_DNAmet.svg", width=5.5, height=1.35)
p3
dev.off()


p1 <- ggplot(df_to_plot, aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=11, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(cols = vars(Modules)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.x = element_text(size = 12, angle=0, colour = "black", face="bold"),strip.background = element_rect(fill="#A2DAC9",linewidth=NA),
        strip.placement = "outside")

p2 <- ggplot(df_to_plot, aes(y=reorder(terms, value), x=value, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=11, vjust=-0.0000000000000001), axis.title.x=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(cols = vars(Modules)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.x = element_text(size = 12, angle=0, colour = "black", face="bold"),strip.background = element_rect(fill="#A2DAC9",linewidth=NA),
        strip.placement = "outside")


svg("C:\\Users\\naata\\Downloads\\wgcnamerged_epi.svg", width=10, height=4.5)
p1
dev.off()



########################################3333 INTERRSECTION




summary_hormone = merge( as.data.frame(table(DEGS_hormone_HS[,colnames(DEGS_hormone_HS)=="term"])), as.data.frame(table(DEGS_hormone_common[,colnames(DEGS_hormone_common)=="term"])), by="Var1", all.x = T, all.y=T)
summary_hormone = merge( summary_hormone, as.data.frame(table(DEGS_hormone_SAHA[,colnames(DEGS_hormone_SAHA)=="term"])), by="Var1", all.x = T, all.y=T)
colnames(summary_hormone)<-c("GOterm", "HeatStress_freq", "Intersection_freq", "HeatStress_SAHA")
summary_hormone[which(is.na(summary_hormone[,2])),2]<-0
summary_hormone[which(is.na(summary_hormone[,3])),3]<-0
summary_hormone[which(is.na(summary_hormone[,4])),4]<-0
summary_hormone$Total <- rowSums(summary_hormone[,c(2,3,4)])
summary_hormone<- summary_hormone[order(summary_hormone$Total, decreasing=T),]




genes_hormone_HDAC_up = table(genes_hormone_HDAC[!(grepl(" -", genes_hormone_HDAC$log2FoldChange)),]$term)
genes_hormone_HDAC_down = table(genes_hormone_HDAC[(grepl(": -", genes_hormone_HDAC$log2FoldChange)),]$term)
genes_hormone_HDAC_merged = data.frame(genes_hormone_HDAC_up)
colnames(genes_hormone_HDAC_merged)<-c("Var1", "Upregulated DEGs")
genes_hormone_HDAC_merged<- merge(genes_hormone_HDAC_merged, hormone_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T)
genes_hormone_HDAC_merged <- melt(genes_hormone_HDAC_merged)
genes_hormone_HDAC_merged$Condition = rep("PE-SAHA", length(genes_hormone_HDAC_merged[,1]))


genes_hormone_HDAC_Stress_up = table(genes_hormone_HDAC_Stress[!(grepl(" -", genes_hormone_HDAC_Stress$log2FoldChange.x)),]$term)
genes_hormone_HDAC_Stress_down = table(genes_hormone_HDAC_Stress[(grepl(" -", genes_hormone_HDAC_Stress$log2FoldChange.x)),]$term)
genes_hormone_HDAC_Stress_merged = merge(data.frame(genes_hormone_HDAC_Stress_up), data.frame(genes_hormone_HDAC_Stress_down), by="Var1", all.x=T, all.y=T)
colnames(genes_hormone_HDAC_Stress_merged)<-c("Var1", "Upregulated DEGs", "Downregulated DEGs")
genes_hormone_HDAC_Stress_merged<- merge(genes_hormone_HDAC_Stress_merged, hormone_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T)
genes_hormone_HDAC_Stress_merged <- melt(genes_hormone_HDAC_Stress_merged)
genes_hormone_HDAC_Stress_merged$Condition = rep("HeatStress-SAHA", length(genes_hormone_HDAC_Stress_merged[,1]))


genes_hormone_Stress_up = table(genes_hormone_Stress[!(grepl(": -", genes_hormone_Stress$log2FoldChange.y)),]$term)
genes_hormone_Stress_down = table(genes_hormone_Stress[(grepl(": -", genes_hormone_Stress$log2FoldChange.y)),]$term)
genes_hormone_Stress_merged = merge(genes_hormone_Stress_up, genes_hormone_Stress_down, by="Var1", all.x=T, all.y=T)
colnames(genes_hormone_Stress_merged)<-c("Var1", "Upregulated DEGs", "Downregulated DEGs")
genes_hormone_Stress_merged<- merge(genes_hormone_Stress_merged, hormone_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T)
genes_hormone_Stress_merged <- melt(genes_hormone_Stress_merged)
genes_hormone_Stress_merged$Condition = rep("HeatStress", length(genes_hormone_Stress_merged[,1]))

genes_hormone_Stress_merged = genes_hormone_Stress_merged[!is.na(genes_hormone_Stress_merged$value),]
df_to_plot = rbind(genes_hormone_Stress_merged, genes_hormone_HDAC_Stress_merged, genes_hormone_HDAC_merged)
df_to_plot$context = gsub(" ", "\n", df_to_plot$context)
p1_h<- ggplot(df_to_plot[df_to_plot$Condition=="HeatStress",], aes(y=reorder(context,order(variable)), x=value, fill=variable)) + xlab("DEGs count") + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=10, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(rows = vars(Condition)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=90, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")
p2_h<- ggplot(df_to_plot[df_to_plot$Condition=="HeatStress-SAHA",], aes(y=reorder(context, value), x=value, fill=variable)) + xlab("DEGs count") + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(rows = vars(Condition)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=90, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")
p3_h<- ggplot(df_to_plot[df_to_plot$Condition=="PE-SAHA",], aes(y=reorder(context, value), x=value, fill=variable)) + xlab("DEGs count") + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=10, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(rows = vars(Condition)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=90, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")


svg("C:\\Users\\naata\\Downloads\\le_h_genes_c.svg", width=9, height=6.2)
ggarrange(p1_h, p2_h, p3_h, nrows=1, ncols=3, widths = c(0.3, 0.5, 0.3), align ="h" )
dev.off()

genes_epi_HDAC_up = table(genes_epi_HDAC[!(grepl(" -", genes_epi_HDAC$log2FoldChange)),]$term)
genes_epi_HDAC_down = table(genes_epi_HDAC[(grepl(": -", genes_epi_HDAC$log2FoldChange)),]$term)
genes_epi_HDAC_merged = merge(genes_epi_HDAC_up, genes_epi_HDAC_down, by="Var1", all.x=T, all.y=T)
colnames(genes_epi_HDAC_merged)<-c("Var1", "Upregulated DEGs", "Downregulated DEGs")
genes_epi_HDAC_merged<- merge(genes_epi_HDAC_merged, epi_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T)
genes_epi_HDAC_merged <- melt(genes_epi_HDAC_merged)
genes_epi_HDAC_merged$Condition = rep("PE-SAHA", length(genes_epi_HDAC_merged[,1]))


genes_epi_HDAC_Stress_up = table(genes_epi_HDAC_Stress[!(grepl(": -", genes_epi_HDAC_Stress$log2FoldChange.x)),]$term)
genes_epi_HDAC_Stress_down = table(genes_epi_HDAC_Stress[(grepl(": -", genes_epi_HDAC_Stress$log2FoldChange.x)),]$term)
genes_epi_HDAC_Stress_merged = merge(genes_epi_HDAC_Stress_up, genes_epi_HDAC_Stress_down, by="Var1", all.x=T, all.y=T)
colnames(genes_epi_HDAC_Stress_merged)<-c("Var1", "Upregulated DEGs","Downregulated DEGs")
genes_epi_HDAC_Stress_merged<- merge(genes_epi_HDAC_Stress_merged, epi_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T)
genes_epi_HDAC_Stress_merged <- melt(genes_epi_HDAC_Stress_merged)
genes_epi_HDAC_Stress_merged$Condition = rep("HeatStress-SAHA", length(genes_epi_HDAC_Stress_merged[,1]))


genes_epi_Stress_up = table(genes_epi_Stress[(grepl(" -", genes_epi_Stress$log2FoldChange.y)),]$term)
genes_epi_Stress_down = table(genes_epi_Stress[(grepl(" -", genes_epi_Stress$log2FoldChange.y)),]$term)
genes_epi_Stress_merged = merge(genes_epi_Stress_up, genes_epi_Stress_down, by="Var1", all.x=T, all.y=T)
colnames(genes_epi_Stress_merged)<-c("Var1", "Upregulated DEGs", "Downregulated DEGs")
genes_epi_Stress_merged<- merge(genes_epi_Stress_merged, epi_terms[,c(2,3)], by.x="Var1", by.y="term", all.x=T)
genes_epi_Stress_merged <- melt(genes_epi_Stress_merged)
genes_epi_Stress_merged$Condition = rep("HeatStress", length(genes_epi_Stress_merged[,1]))


df_to_plot = rbind(genes_epi_Stress_merged, genes_epi_HDAC_Stress_merged, genes_epi_HDAC_merged)
df_to_plot$context = gsub(" driven by", "\ndriven by", df_to_plot$context)
df_to_plot$context = gsub(" and", "\nand", df_to_plot$context)
df_to_plot$context = gsub("dependent ", "dependent\n", df_to_plot$context)

p1_h<- ggplot(df_to_plot[df_to_plot$Condition=="HeatStress",], aes(y=reorder(context, value), x=value, fill=variable)) + xlab("DEGs count") + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(rows = vars(Condition)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 9, angle=90, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside")
p2_h<- ggplot(df_to_plot[df_to_plot$Condition=="HeatStress-SAHA",], aes(y=reorder(context, value), x=value, fill=variable)) + xlab("DEGs count") + geom_bar(stat="identity", position="stack")+ scale_x_continuous(breaks=c(0,10,20,30,40,50)) + theme_bw() + 
  theme(axis.text.y = element_text(size=10, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank()) +
  facet_grid(rows = vars(Condition)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 9, angle=90, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside")
p3_h<- ggplot(df_to_plot[df_to_plot$Condition=="PE-SAHA",], aes(y=reorder(context, value), x=value, fill=variable)) + xlab("DEGs count") + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.y = element_text(size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(Condition)) + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=90, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#e3eeec",linewidth=NA),
        strip.placement = "outside")

svg("C:\\Users\\naata\\Downloads\\le_pesaha_hgenes.svg", width=8, height=1.85)
p3_h
dev.off()


svg("C:\\Users\\naata\\Downloads\\le_epi_genes_c.svg", width=12, height=9)
ggarrange(p1_h, p2_h, p3_h, nrows=1, ncols=3, widths = c(0.55, 0.60, 0.6), align ="h" )
dev.off()

svg("C:\\Users\\naata\\Downloads\\le_epi_genes.svg", width=23, height=10.5)
ggarrange(p1_h, p2_h, p3_h, nrows=3, widths=c(0.45, 0.5, 0.4), common.legend = T, legend="right" )
dev.off()


table(genes_hormone_HDAC_Stress$context)

table(genes_hormone_Stress$context)















################### unique distinct SAHA genes

SAHA_DEGs_h_nc = SAHA_DEGs_h[!(SAHA_DEGs_h$external_gene_name %in% HS_DEGs_h$external_gene_name),]
SAHA_DEGs_h_nc_UP = data.frame(table(SAHA_DEGs_h_nc[!(grepl("-" ,SAHA_DEGs_h_nc$log2FoldChange)),]$term))
SAHA_DEGs_h_nc_DOWN = data.frame(table(SAHA_DEGs_h_nc[(grepl("-" ,SAHA_DEGs_h_nc$log2FoldChange)),]$term))
SAHA_DEGs_h_nc_UP = merge(SAHA_DEGs_h_nc_UP, hormone_terms[,c(2,3)], by.x="Var1", by.y="term")
SAHA_DEGs_h_nc_DOWN = merge(SAHA_DEGs_h_nc_DOWN, hormone_terms[,c(2,3)], by.x="Var1", by.y="term")



SAHA_DEGs_h_nc_UP$variable = rep("Upregulated DEGs",length(SAHA_DEGs_h_nc_UP[,1]))
SAHA_DEGs_h_nc_DOWN$variable = rep("Downregulated DEGs",length(SAHA_DEGs_h_nc_DOWN[,1]))

df_plot_terms = rbind(SAHA_DEGs_h_nc_UP, SAHA_DEGs_h_nc_DOWN)
write.table(df_plot_terms, "C:\\Users\\naata\\Downloads\\freqdegs_hnc.tsv", sep="\t", row.names = F, quote=F)




SAHA_DEGs_epi_nc = SAHA_DEGs_epi[!(SAHA_DEGs_epi$external_gene_name %in% HS_DEGs_epi$external_gene_name),]

SAHA_DEGs_epi_UP = data.frame(table(SAHA_DEGs_epi_nc[!(grepl("-" ,SAHA_DEGs_epi_nc$log2FoldChange)),]$term))
SAHA_DEGs_epi_DOWN = data.frame(table(SAHA_DEGs_epi_nc[(grepl("-" ,SAHA_DEGs_epi_nc$log2FoldChange)),]$term))
SAHA_DEGs_epi_UP = merge(SAHA_DEGs_epi_UP, epi_terms[,c(2,3)], by.x="Var1", by.y="term")
SAHA_DEGs_epi_DOWN = merge(SAHA_DEGs_epi_DOWN, epi_terms[,c(2,3)], by.x="Var1", by.y="term")

SAHA_DEGs_epi_UP$variable = rep("Upregulated DEGs",length(SAHA_DEGs_epi_UP[,1]))
SAHA_DEGs_epi_DOWN$variable = rep("Downregulated DEGs",length(SAHA_DEGs_epi_DOWN[,1]))
df_plot_terms = rbind(SAHA_DEGs_epi_UP, SAHA_DEGs_epi_DOWN)

write.table(df_plot_terms, "C:\\Users\\naata\\Downloads\\freqdegs_nc.tsv", sep="\t", row.names = F, quote=F)


df_plot_terms_ann = read.csv( "C:\\Users\\naata\\Downloads\\freqdegs_nc.tsv", sep="\t")
df_plot_terms_ann$context = gsub(" ", "\n", df_plot_terms_ann$context)

svg("C:\\Users\\naata\\Downloads\\barplot_DEG_SAHAdistinct_h.svg", width=5, height=8)
ggplot(df_plot_terms_ann, aes(y=term, x=Freq, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(context)) + xlab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 12, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")
dev.off()

#write.table(df_plot_terms, "C:\\Users\\naata\\Downloads\\freqdegs_nc.tsv", sep="\t", row.names = F, quote=F)
df_plot_terms_ann = read.csv( "C:\\Users\\naata\\Downloads\\freqdegs_nc.tsv", sep="\t")
df_plot_terms_ann$term = gsub(" ", "\n", df_plot_terms_ann$term)


#svg("C:\\Users\\naata\\Downloads\\barplot_DEG_distinctSAHA.svg", width=5, height=6.5)
p1 <- ggplot(df_plot_terms_ann[df_plot_terms_ann$term=="chromatin\nremodelling",], aes(y=reorder(context, Freq), x=Freq, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank(), axis.title.x=element_blank()) +
  facet_grid(rows = vars(term)) + xlab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")

p2<- ggplot(df_plot_terms_ann[df_plot_terms_ann$term=="Histone\nmodification",], aes(y=reorder(context, Freq), x=Freq, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "none", axis.title.y=element_blank(), axis.title.x=element_blank()) +
  facet_grid(rows = vars(term)) + xlab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")

p3<-ggplot(df_plot_terms_ann[df_plot_terms_ann$term=="DNA\nmethylation",], aes(y=reorder(context, Freq), x=Freq, fill=variable)) + geom_bar(stat="identity", position="stack")  + theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=9, vjust=-0.0000000000000001), legend.title = element_blank(), legend.position = "bottom", axis.title.y=element_blank()) +
  facet_grid(rows = vars(term)) + xlab("DEG gene count") + 
  scale_fill_manual(values=c("Downregulated DEGs" = "#00A36C", "Upregulated DEGs" ="coral")) +                                                             # Change font size
  theme(strip.text.y = element_text(size = 10, angle=0, colour = "#a8b1b0", face="bold"),strip.background = element_rect(fill="#f0e9c9",linewidth=NA),
        strip.placement = "outside")
#dev.off()

svg("C:\\Users\\naata\\Downloads\\barplot_DEG_SAHAdistinct_epi.svg", width=5, height=5.2)
ggarrange(p1,p2,p3,nrow=3, ncol=1, common.legend = T, legend = "bottom", align="hv", heights = c(0.25,0.4,0.27))
dev.off()


# save new lists
write.table(SAHA_DEGs_h_nc, "C:\\Users\\naata\\Downloads\\SAHA_distinct_hormone.tsv", sep="\t", row.names=F, quote=F)
write.table(SAHA_DEGs_epi_nc, "C:\\Users\\naata\\Downloads\\SAHA_distinct_epi.tsv", sep="\t", row.names=F, quote=F)
write.table(DEGS_hormone_common, "C:\\Users\\naata\\Downloads\\common_hormone.tsv", sep="\t", row.names=F, quote=F)
write.table(DEGS_epi_common, "C:\\Users\\naata\\Downloads\\common_epi.tsv", sep="\t", row.names=F, quote=F)
