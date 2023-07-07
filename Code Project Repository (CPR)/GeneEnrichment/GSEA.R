#!/bin/usr/env Rscript

# Author: Natalia García Sánchez
# Description : index crafted gene category association file for terms of interest, get gene-category count.
# Objective: generate GSEA results for categories of interest


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





