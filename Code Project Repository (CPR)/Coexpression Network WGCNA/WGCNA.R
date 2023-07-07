# script to perform WGCNA
# Code adapted from kpatel427 Youtube tutorials repository - WGCNA.R script https://github.com/kpatel427/YouTubeTutorials/blob/main/WGCNA.R


# Install biomcmanager dependencies

install.packages("BiocManager")
BiocManager::install(version = "3.17")
#install.packages("org.Hs.eg.db", repos="http://bioconductor.org/packages/3.2/data/annotation")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("biomartr")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("minet")
remotes::install_github("jtlovell/limmaDE2") # requires install.packages("venneuler")  as a dependency

#remotes::install_github("kevinblighe/CorLevelPlot")
#install.packages("WGCNA", dependencies = TRUE)

         # allow multi-threading (optional)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(biomartr)
library(biomaRt)
library(CorLevelPlot)
library(gridExtra)
library(reshape)
library(ggplot2)
library(ggpubr)
library(igraph)
library("ggVennDiagram") # requires c pckg - libudunits2-dev, libgdal-dev, units, sf CRAN packages
library("AnnotationDbi")
library("AnnotationHub")
allowWGCNAThreads() 
######################################################
# Load bnapus object
######################################################

# Making annotation package
ah <- AnnotationHub()
#Ensembl 56
query(ah, c( "OrgDb", "Brassica napus"))
org.Brassicanapus.eg.db <- ah[["AH107389"]]

#############################################################
# Load ncbi annotations
annot_ncbi=read.table("/home/famgarcia/Escritorio/ncbi_centric_annotations.tsv", sep="\t", header=TRUE)

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

################
# MAIN
################




# 1. Fetch Data ------------------------------------------------
data <- read.csv("/home/famgarcia/Escritorio/genes_all_experiments_batch_correction.csv", sep=",")
data$ensembl_id = rownames(data)
rownames(data) = NULL

# 2. Filter by having Entrez ID

filter <- read.csv("/home/famgarcia/Descargas/mart_export.csv", sep=",")
colnames(filter)<-c("ensembl_id","gid")
filter = filter[!is.na(filter$gid),]

# Inner join
data_eid = merge(data, filter, by = "ensembl_id")



# Remove genecount duplicates - transcript isoforms discarded for study (already took genes with longest transcripts in DEA analysis)

data_eid = data_eid[!duplicated(data_eid$gid),]
data_eid$ensembl_id = NULL

# prepare data
rownames(data_eid) = data_eid$gid
data_eid$gid = NULL

ggplot(melt(data_eid), aes(x=variable, y=value)) + geom_boxplot()

hist(data_eid)
summary(data_eid)

# joined metadata from :
# - projects NVUK2022082921-EU-ES-CSICCIB-RNAseq-6-12G-Pilar-WBI - SAHA 
# - projects NVUK2022012116-EU-ES-CSIC-CIB-RNAseq-6-12G-Pilar-WBI - Stress

phenoData <- data.frame(title = colnames(data_eid),
                        Treatment = c(rep("Heat_Stress_32C", 3), rep("Stress_SAHA", 3), rep("Heat_Stress_32C", 5)),
                        DevelopmentalStage = c(rep("Proembryo", 6), rep("VacuolatedMicrospore", 3), rep("Proembryo", 2)),
                        report_time = c(rep("2022-10-06", 6), rep("2022-03-09", 5)),
                        cultivar_cell_line = c(rep("Topas_DH4079")))
rownames(phenoData) = colnames(data_eid)

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data_eid))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data_eid[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
par(mar=c(0,4,4,8))
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)
pca.dat$Treatment = c("Heat Stress","Heat Stress","Heat Stress","Heat Stress + SAHA","Heat Stress + SAHA","Heat Stress + SAHA","Heat Stress","Heat Stress","Heat Stress","Heat Stress","Heat Stress")
pca.dat$Dev_Stage = c("Proembryo","Proembryo","Proembryo","Proembryo","Proembryo","Proembryo","Vacuolated Microspore","Vacuolated Microspore","Vacuolated Microspore","Proembryo","Proembryo")

svg("/home/famgarcia/TFM/BATCH_CORRECTION_PCA_preWGCNA_prefiltering.svg")
ggplot(pca.dat, aes(PC1, PC2, color=Dev_Stage, shape=Treatment)) +
  geom_point() +
  geom_text(label = rownames(pca.dat), vjust=1.5) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' % Explained Variance'),
       y = paste0('PC2: ', pca.var.percent[2], ' % Explained Variance')) + 
  scale_x_continuous(limits = c(-370000,480000))
dev.off()



# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

colData <- phenoData 
data.subset <- data

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


# create dds to perform a Variance Stabilizing Transformation Normalization from the counts
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (11*0.75=9.105)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 9,]
nrow(dds75) # 6998 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices
sft.data$color=ifelse(sft.data$Power ==12, "red", "black")

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power, color=color)) +
  geom_point() +
  geom_text(size=2.5,nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = bquote('Scale free topology model fit, signed '~R^2) ) +
  scale_color_manual(values=c('black','red')) +
  theme(legend.position = "none") +
  theme_classic(base_size = 7) + guides(color="none")


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power, color=color)) +
  geom_point( show.legend = FALSE ) +
  geom_text(size=2.5, nudge_y = 0.1, vjust = -2, show.legend = FALSE) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  scale_color_manual(values=c('black','red')) +
  theme(legend.position = "none") +
  theme_classic(base_size = 7)

png("/home/famgarcia/Escritorio/NTA_softthreshold_options.png", width=24, height=7, units="cm", res=300)
grid.arrange(a1, a2, ncol = 2)
dev.off()


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 12
temp_cor <- cor
cor <- WGCNA::cor



# ONE-STEP NETWORK CONSTRUCTION AND MODULE DETECTION memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          minModuleSize = 30,
                          maxPOutliers = 0.05,
                          maxBlockSize = 8000,
                          networkType = "signed",
                          TOMType = "signed",
                          corType="bicor",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          saveTOMs = TRUE,
                          #saveTOMFileBase = "/home/famgarcia/TFM/TOM_DS_COV",
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations


# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(SAHA_Proembryos = ifelse(grepl('SAHA', Treatment), 1, 0)) %>%
  mutate(Ctrl_Proembryos = c(1,1,1,0,0,0,0,0,0,1,1)) %>%
  mutate(Vacuolated_Micrp = ifelse(grepl('Vacuolated', DevelopmentalStage), 1, 0))

traits<-traits[,c(6,7,8)]



# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
textMatrix = paste( signif( module.trait.corr[,], 2 ), "  ;  (" , # Pasting corr to pvalues
                    signif( module.trait.corr.pvals[,], 1), 
                    ifelse(module.trait.corr.pvals[,]<0.05, "*",""), 
                    ifelse(module.trait.corr.pvals[,]<0.01, "*",""), 
                    ifelse(module.trait.corr.pvals[,]<0.001, "*",""),")", sep="")
dim(textMatrix)=dim(module.trait.corr[,])
textMatrix=as.matrix(textMatrix)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


sign_modules = rownames(module.trait.corr.pvals[rowSums(module.trait.corr.pvals < 0.05)>0,])
sign_colors = gsub("ME", "",sign_modules)


x11(type="cairo")
svg("/home/famgarcia/Escritorio/modcorr_signednet.svg", width=15, height=23)
par(mar=c(16,16,16,16))
labeledHeatmap(Matrix=module.trait.corr[,],
               xLabels=t(as.data.frame(names(heatmap.data)[c(21,22,23)])),
               yLabels=names(heatmap.data)[1:20],
               ySymbols=gsub("ME","",names(heatmap.data)[1:20]),
               colorLabels = FALSE,
               cex.lab = 2.5,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins = FALSE,
               cex.text=1.9,
               zlim=c(-1,1)) + geom_raster() + scale_fill_identity() + theme(plot.margin = c(8,8,8,8))
dev.off()

# lightyellow module is of biological interest to analyze SAHA treatment effects on corregulation network (significant module associated to that trait)
# grey60 module is of biological interest to give insights of coexpression in Control Proembryos
# brown, grey60, turquoise blue yellow might be informative to analyze coexpression clusters in Vacuolated Microspores before embryogenesis


# Control Proembryos
module.gene.mapping <- as.data.frame(bwnet$colors)
grey60_genes = module.gene.mapping %>% 
  filter(`bwnet$colors` == 'grey60') %>% 
  rownames()

grey60_genes = as.data.frame(grey60_genes)
colnames(grey60_genes) ="entrezgene_id"
merge(grey60_genes,annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], by="entrezgene_id")


# SAHA Proembryos
lightyellow_genes = module.gene.mapping %>% 
  filter(`bwnet$colors` == 'lightyellow') %>% 
  rownames()

lightyellow_genes = as.data.frame(lightyellow_genes)
colnames(lightyellow_genes) ="entrezgene_id"
merge(lightyellow_genes,annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], by="entrezgene_id")


# Vacuolated Microspores
brown_genes = module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()

brown_genes = as.data.frame(brown_genes)
colnames(brown_genes) ="entrezgene_id"
merge(brown_genes,annot_ncbi[!duplicated(annot_ncbi$entrezgene_id),], by="entrezgene_id")





#  GENE SIGNIFICANCE MODULE MEMBERSHIP ----------------------
# or each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profil


# 6B. Intramodular analysis: Identifying driver genes ---------------



# I. MODULE MEMBERSHIP. Calculate the module membership and the associated p-values-------------------------

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

# Gene significance analysis - not for traits defined but for modules - module eigengenes.

geneTraitSignificance = as.data.frame(cor(norm.counts, module_eigengenes, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

module.membership.measure.pvals[1:10,1:10]



# II. GENE SIGNIFICANCE  -------------------------------------------
# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.
# Calculate the gene significance and associated p-values

################### SAHA TREATMENT

# Plot for modules of interest?


gene.signf.corr <- cor(norm.counts, traits$SAHA_Proembryos, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)
colnames(gene.signf.corr.pvals)=c("pvalue")
significant_genes_treatment = subset(gene.signf.corr.pvals, pvalue<=0.01)

# Plot correlation between module membership and trait significance in hub genes
modNames = rownames(module.membership.measure)
module = "lightyellow"
  column = match(paste("ME",module,sep=""), modNames);
  associated_mod_genes <- module.gene.mapping %>% 
    filter(`bwnet$colors` == 'lightyellow') %>% 
    rownames()
png("/home/famgarcia/TFM/wgcna/mm_traitCPE_lightyellow.png", width=15, height=10, units="cm", res=250)
par(mar=c(5,5,2.6,3))
verboseScatterplot(abs(module.membership.measure[column, associated_mod_genes]),
                   abs(gene.signf.corr[associated_mod_genes, 1]),
                   xlab = "Module Membership in lightyellow module",
                   ylab = "Gene significance \nfor Ctrl Proembryos",
                   main = paste("Module membership vs. trait gene significance\n"),
                   cex.main = 1, cex.lab = 0.9, cex.axis = 0.9, col = '#FBEC5D') +theme(axis.title.x = element_text("Module Membership in grey60 module", size=7))+ theme_bw(base_size=4)
dev.off()



################## DEVELOPMENTAL STAGE - Ctrl proembryos

gene.signf.corr <- cor(norm.counts, traits$Ctrl_Proembryos, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)
colnames(gene.signf.corr.pvals)=c("pvalue")
significant_genes_stagedev_PE = subset(gene.signf.corr.pvals, pvalue<=0.01)
significant_genes_stagedev_PE$entrezgene_id = rownames(significant_genes_stagedev_PE)
significant_genes_stagedev_PE =significant_genes_stagedev_PE[order(significant_genes_stagedev_PE$pvalue),]

modNames = rownames(module.membership.measure)
module = "grey60"
  column = match(paste("ME",module,sep=""), modNames);
  associated_mod_genes <- module.gene.mapping %>% 
    filter(`bwnet$colors` == 'grey60') %>% 
    rownames()

png("/home/famgarcia/TFM/wgcna/mm_traitCPE_grey60.png", width=15, height=10, units="cm", res=250)
par(mar=c(5,5,2.6,3))
verboseScatterplot(abs(module.membership.measure[column, associated_mod_genes]),
                   abs(gene.signf.corr[associated_mod_genes, 1]),
                   xlab = "Module Membership in grey60 module",
                   ylab = "Gene significance \nfor Ctrl Proembryos",
                   main = paste("Module membership vs. trait gene significance\n"),
                   cex.main = 1, cex.lab = 0.9, cex.axis = 0.9, col = module) +theme(axis.title.x = element_text("Module Membership in grey60 module", size=7))+ theme_bw(base_size=4)
dev.off()



################## DEVELOPMENTAL STAGE - VM

gene.signf.corr <- cor(norm.counts, traits$Vacuolated_Micrp, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)
colnames(gene.signf.corr.pvals)=c("pvalue")
significant_genes_stagedev_VM = subset(gene.signf.corr.pvals, pvalue<=0.01)

modNames = rownames(module.membership.measure)
module = "brown"
  column = match(paste("ME",module,sep=""), modNames);
  associated_mod_genes <- module.gene.mapping %>% 
    filter(`bwnet$colors` == module) %>% 
    rownames()
  
png("/home/famgarcia/mme_traitVM_brown.png", width=15, height=10, units="cm", res=250)
verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                   abs(gene.signf.corr[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Vacuolated Microspore",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

modNames = rownames(module.membership.measure)
module = "turquoise"
  column = match(paste("ME",module,sep=""), modNames);
  associated_mod_genes <- module.gene.mapping %>% 
    filter(`bwnet$colors` == module) %>% 
    rownames()
  
png("/home/famgarcia/mme_traitVM_turquoise.png", width=15, height=10, units="cm", res=250)
verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                     abs(gene.signf.corr[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Vacuolated Microspore",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

modNames = rownames(module.membership.measure)
module = "yellow"
  column = match(paste("ME",module,sep=""), modNames);
  associated_mod_genes <- module.gene.mapping %>% 
    filter(`bwnet$colors` == module) %>% 
    rownames()


png("/home/famgarcia/mme_traitVM_yellow.png", width=15, height=10, units="cm", res=250)
verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                   abs(gene.signf.corr[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Vacuolated Microspore",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

modNames = rownames(module.membership.measure)
module = "blue"
  column = match(paste("ME",module,sep=""), modNames);
  associated_mod_genes <- module.gene.mapping %>% 
    filter(`bwnet$colors` == module) %>% 
    rownames()
  

  png("/home/famgarcia/mme_traitVM_blue.png", width=15, height=10, units="cm", res=250)
  verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                     abs(gene.signf.corr[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Vacuolated Microspore",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()

  
  modNames = rownames(module.membership.measure)
  module = "lightyellow"
    column = match(paste("ME",module,sep=""), modNames);
    associated_mod_genes <- module.gene.mapping %>% 
      filter(`bwnet$colors` == module) %>% 
      rownames()
    
    verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                       abs(gene.signf.corr[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for Vacuolated Microspore",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
    
    modNames = rownames(module.membership.measure)
    module = "grey60"
      column = match(paste("ME",module,sep=""), modNames);
      associated_mod_genes <- module.gene.mapping %>% 
        filter(`bwnet$colors` == module) %>% 
        rownames()
      
      verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                         abs(gene.signf.corr[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance for Vacuolated Microspore",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      
  
  write.table(rownames(significant_genes_stagedev_VM), "hub_VM.tsv", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
  write.table(rownames(significant_genes_stagedev_PE), "hub_PEC.tsv", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
  write.table(rownames(significant_genes_treatment), "hub_SAHA.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
  
# Annotate results


Normalized_val = read.csv("/home/famgarcia/Escritorio/VST_norm_counts.csv", sep=" ")
Normalized_val$entrezgene_id = rownames(Normalized_val)
rownames(Normalized_val) <- NULL


# NCBI GENE DESCRIPTIONS

# Expr log2fc among SAHA treated PE - Control PE
significant_genes_treatment$entrezgene_id = rownames(significant_genes_treatment)
rownames(significant_genes_treatment) = NULL
hub_PEall = as.data.frame(significant_genes_stagedev_PE)
putative_hubgenes_involved_development = merge(hub_PEall,annot_ncbi, by='entrezgene_id')
putative_hubgenes_involved_development = putative_hubgenes_involved_development[order(putative_hubgenes_involved_development$pvalue),]

tmp_GO = select(org.Brassicanapus.eg.db, as.character(putative_hubgenes_involved_development$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
tmp_GO = tmp_GO[tmp_GO$ONTOLOGY=="BP",]
xx <- as.list(tmp_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
tmp_GO$GOTERM = list_terms
tmp_GO= tmp_GO[!is.na(tmp_GO$GO),]
putative_hubgenes_involved_development = merge(tmp_GO, putative_hubgenes_involved_development, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
expr_phgenes_treatment = merge(putative_hubgenes_involved_development, expr_HDAC_Stress, by="external_gene_name")
expr_phgenes_treatment = expr_phgenes_treatment[order(expr_phgenes_treatment$pvalue),]

expr_phgenes_treatment = merge(expr_phgenes_treatment, Normalized_val, by.x="ENTREZID", by.y="entrezgene_id")
write.table(expr_phgenes_treatment, "/home/famgarcia/Escritorio/hub_Stress_PEall.tsv", sep="\t", row.names=FALSE, quote=FALSE)

######################## VM hub genes

hub_VM = as.data.frame(significant_genes_stagedev_VM)
hub_VM$entrezgene_id = rownames(hub_VM)
putative_hubgenes_involved_development = merge(hub_VM,annot_ncbi, by='entrezgene_id')
putative_hubgenes_involved_development = putative_hubgenes_involved_development[order(putative_hubgenes_involved_development$pvalue),]

tmp_GO = select(org.Brassicanapus.eg.db, as.character(putative_hubgenes_involved_development$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
tmp_GO = tmp_GO[tmp_GO$ONTOLOGY=="BP",]
xx <- as.list(tmp_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
tmp_GO$GOTERM = list_terms
tmp_GO= tmp_GO[!is.na(tmp_GO$GO),]
putative_hubgenes_involved_development = merge(tmp_GO, putative_hubgenes_involved_development, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)


expr_phgenes_dev = merge(putative_hubgenes_involved_development, expr_Stress, by="external_gene_name")
expr_phgenes_dev = expr_phgenes_dev[order(expr_phgenes_dev$pvalue),]
expr_phgenes_treatment = merge(expr_phgenes_dev, Normalized_val, by.x="ENTREZID", by.y="entrezgene_id")

write.table(expr_phgenes_treatment, "/home/famgarcia/Escritorio/hub_Stress_VM.tsv", sep="\t", row.names=FALSE, quote=FALSE)


######## SAHA treatment

hub_SAHA = as.data.frame(significant_genes_treatment)
hub_SAHA$entrezgene_id = rownames(hub_SAHA)
rownames(hub_SAHA) <- NULL
putative_hubgenes_involved_treatment = merge(hub_SAHA,annot_ncbi, by='entrezgene_id')


tmp_GO = select(org.Brassicanapus.eg.db, as.character(putative_hubgenes_involved_treatment$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
tmp_GO = tmp_GO[tmp_GO$ONTOLOGY=="BP",]
xx <- as.list(tmp_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
tmp_GO$GOTERM = list_terms
tmp_GO= tmp_GO[!is.na(tmp_GO$GO),]
putative_hubgenes_involved_treatment = merge(tmp_GO, putative_hubgenes_involved_treatment, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
expr_phgenes_treatment = merge(putative_hubgenes_involved_treatment, expr_HDAC, by="external_gene_name")
expr_phgenes_treatment = expr_phgenes_treatment[order(expr_phgenes_treatment$pvalue),]
expr_phgenes_treatment = merge(expr_phgenes_treatment, Normalized_val, by.x="ENTREZID", by.y="entrezgene_id")

write.table(expr_phgenes_treatment, "/home/famgarcia/Escritorio/hub_Stress_SAHA.tsv", sep="\t", row.names=FALSE, quote=FALSE)

### Hub genes per module

hub_genes_modules = chooseTopHubInEachModule(norm.counts, bwnet$colors, power=12,type = "signed")
write.table(hub_genes_modules,"/home/famgarcia/Escritorio/hub_genes_per_module.csv", sep=",")
hub_genes_modules=read.table("/home/famgarcia/Escritorio/hub_genes_per_module.csv", sep=",")
hub_genes_modules$module = rownames(hub_genes_modules)
colnames(hub_genes_modules) = "entrezgene_id"
hub_genes_modules=merge(hub_genes_modules, annot_ncbi, by="entrezgene_id", all.x=TRUE)


tmp_GO = select(org.Brassicanapus.eg.db, as.character(hub_genes_modules$entrezgene_id), c("GO" , "ONTOLOGY", "GENENAME"), "ENTREZID")
tmp_GO = tmp_GO[tmp_GO$ONTOLOGY=="BP",]
xx <- as.list(tmp_GO$GO)
list_terms = c()
for (i in 1:length(xx)){
  if (!is.na(xx[[i]])) { 
    list_terms=c(list_terms,  as.character( Term( xx[[i]] ) )) }
  else{list_terms =c(list_terms, "")}
} 
tmp_GO$GOTERM = list_terms
tmp_GO= tmp_GO[!is.na(tmp_GO$GO),]
hub_genes_modules = merge(tmp_GO, hub_genes_modules, by.x="ENTREZID", by.y="entrezgene_id", all.y=TRUE)
hub_genes_modules = merge(hub_genes_modules, Normalized_val, by.x="ENTREZID", by.y="entrezgene_id")

write.table(hub_genes_modules, "/home/famgarcia/Escritorio/Tophub_all_modules.tsv", sep="\t", row.names=FALSE, quote=FALSE)


###############################################################################################
library(limmaDE2)
graph<-wgcna2igraph(net = bwnet, datExpr = norm.counts[],
                    modules2plot = c("black", "greenyellow", "brown", "turquoise", "midnightblue","blue"),
                    colors2plot = c("gold","lightgreen", "brown","cyan","cornflowerblue", "coral"),
                    kME.threshold = 0.5, adjacency.threshold = 0.1,
                    adj.power = pow, verbose = T,
                    node.size = 0, frame.color = NA, node.color = NA,
                    edge.alpha = .5, edge.width =1)
plot(graph)



##############################################################################################
####################### EXPORT TO CYTOSCAPE
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

table(bwnet$colors[bwnet$colors %in% sign_colors])


# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(norm.counts, power = 12);

# Select module
module = "brown";
# Select module probes
probes = names(norm.counts)
inModule = (bwnet$colors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
# vis = exportNetworkToVisANT(modTOM,
#                             file = paste("VisANTInput-", module, ".txt", sep=""),
#                             weighted = TRUE,
#                             threshold = 0,
#                             probeToGene = data.frame(annot_ncbi$entrezgene_id, annot_ncbi$ALIAS) )
# 


nTop = 30;
IMConn = softConnectivity(norm.counts[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  GET TOM - WEIGHTED ADJACENCY M FROM INTERESTING MODULES - CONVERT TO EDGE/NODE FILES FOR CYTOSCAPE
#
#=====================================================================================


# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(norm.counts, power = 12);

# Select modules
#modules = c("grey60", "lightyellow", "yellow", "blue");
modules = c("grey60", "lightyellow", "yellow", "blue");
modules = sign_colors

# Select module probes
probes = colnames(norm.counts)
inModule = is.finite(match(bwnet$colors, modules));
modProbes = probes[inModule];
modGenes = annot_ncbi$GENENAME[match(modProbes, annot_ncbi$entrezgene_id)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/home/famgarcia/Escritorio/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/home/famgarcia/Escritorio/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);


# save hubs from all exp
modules=c("PE_gs_01","PE_SAHA_gs_01","VM_gs_01")
modColors = as.data.frame(rbind( cbind(hub_PEall$entrezgene_id, rep("PE_gs_01",length(hub_PEall$entrezgene_id))),
cbind(hub_SAHA$entrezgene_id, rep("PE_SAHA_gs_01",length(hub_SAHA$entrezgene_id))),
cbind(hub_VM$entrezgene_id, rep("VM_gs_01",length(hub_VM$entrezgene_id))) ))
colnames(modColors)<-c("gid","annot")
rownames(modColors)<-modColors$gid
modColors$gid<-NULL

modProbes = unique(modColors$gid)

modGenes = annot_ncbi$ALIAS[match(modProbes, annot_ncbi$entrezgene_id)];
inModules =  colnames(norm.counts) %in% modProbes
# Select the corresponding Topological Overlap
modTOM = TOM[inModules, inModules];
dimnames(modTOM) = list(modProbes, modProbes)


# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/home/famgarcia/Escritorio/CytoscapeInput-edges_gs_traits-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/home/famgarcia/Escritorio/CytoscapeInput-nodes_gs_traits-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = modColors[!duplicated(modColors$gid),]);


####################################################################################################

# VISUALIZATION


#=====================================================================================
#
#   VISUALIZE module relationships within traits  - EIGENGENE profile similarity
#
#=====================================================================================


# Recalculate module eigengenes
MEs = moduleEigengenes(norm.counts, bwnet$colors)$eigengenes
# Isolate SAHA treated proembryo weight from the traits
weight = as.data.frame(traits$SAHA_Proembryos);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
png("/home/famgarcia/Escritorio/Modulecorr_SAHA_PE.png")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()


# Isolate Ctrl proembryo weight from the traits
weight = as.data.frame(traits$Ctrl_Proembryos);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
png("/home/famgarcia/Escritorio/Modulecorr_Ctrl_PE.png")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

# Isolate weight from the clinical traits
weight = as.data.frame(traits$Vacuolated_Micrp);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
png("/home/famgarcia/Escritorio/Modulecorr_VM.png")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()




#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

moduleColors = bwnet$colors
# Match probes in the data set to the probe IDs in the annotation file 
probes = colnames(norm.counts)
allLLIDs = probes

#probes2annot = match(probes, annot_ncbi$entrezgene_id)
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$ALIAS[probes2annot];
# $ Choose interesting modules

intModules = c("grey60", "lightyellow", "yellow", "blue","brown", "turquoise");

for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("SignNet_genes_modules_interest-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("SignNet_genes_modules_interest-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

#=====================================================================================
#
#  GET TOM - WEIGHTED ADJACENCY M FROM INTERESTING MODULES - CONVERT TO EDGE/NODE FILES FOR CYTOSCAPE
#
#=====================================================================================


# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(norm.counts, power = 12);

hub_
unique(c(hub_SAHA$entrezgene_id, hub_VM$entrezgene_id, hub_PEall$entrezgene_id))


# Select modules
modules = c("lightyellow", "grey60");
modules = sign_colors

# Select module probes
probes = colnames(norm.counts)
inModule = is.finite(match(bwnet$colors, modules));
modProbes = probes[inModule];
modGenes = annot_ncbi$GENENAME[match(modProbes, annot_ncbi$entrezgene_id)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/home/famgarcia/Escritorio/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/home/famgarcia/Escritorio/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);







# Repeat for Development Stage covariate batch adjustment


# 1. Fetch Data ------------------------------------------------
data <- read.csv("/home/famgarcia/Escritorio/genes_all_experiments_batch_correction.csv", sep=",")
data<-data[, c(1:6, 10,11)]
data$ensembl_id = rownames(data)
rownames(data) = NULL

# 2. Filter by having Entrez ID

filter <- read.csv("/home/famgarcia/Descargas/mart_export.csv", sep=",")
colnames(filter)<-c("ensembl_id","gid")
filter = filter[!is.na(filter$gid),]

# Inner join
data_eid = merge(data, filter, by = "ensembl_id")


# Remove genecount duplicates - transcript isoforms discarded for study (already took genes with longest transcripts in DEA analysis)

data_eid = data_eid[!duplicated(data_eid$gid),]
data_eid$ensembl_id = NULL

# prepare data
rownames(data_eid) = data_eid$gid
data_eid$gid = NULL


# joined metadata from :
# - projects NVUK2022082921-EU-ES-CSICCIB-RNAseq-6-12G-Pilar-WBI - SAHA 
# - projects NVUK2022012116-EU-ES-CSIC-CIB-RNAseq-6-12G-Pilar-WBI - Stress

phenoData <- data.frame(title = colnames(data_eid),
                        Treatment = c(rep("Heat_Stress_32C", 3), rep("Stress_SAHA", 3), rep("Heat_Stress_32C", 2)),
                        report_time = c(rep("2022-10-06", 6), rep("2022-03-09", 2)),
                        cultivar_cell_line = c(rep("Topas_DH4079")))
rownames(phenoData) = colnames(data_eid)

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data_eid))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data_eid[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
par(mar=c(0,4,4,8))
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)
pca.dat$Treatment = c("Heat Stress","Heat Stress","Heat Stress","Heat Stress + SAHA","Heat Stress + SAHA","Heat Stress + SAHA","Heat Stress","Heat Stress")
pca.dat$Dev_Stage = c("Proembryo","Proembryo","Proembryo","Proembryo","Proembryo","Proembryo","Proembryo","Proembryo")

ggplot(pca.dat, aes(PC1, PC2, color=Dev_Stage, shape=Treatment)) +
  geom_point() +
  geom_text(label = rownames(pca.dat), vjust=1.5) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' % Explained Variance'),
       y = paste0('PC2: ', pca.var.percent[2], ' % Explained Variance')) + 
  scale_x_continuous(limits = c(-370000,480000))



# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

colData <- phenoData 
data.subset <- data

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


# create dds to perform a Variance Stabilizing Transformation Normalization from the counts
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (9*0.75=6.75)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 7,]
nrow(dds75) # 6998 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1, vjust = 1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor



# ONE-STEP NETWORK CONSTRUCTION AND MODULE DETECTION memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          minModuleSize = 30,
                          maxPOutliers = 0.05,
                          maxBlockSize = 8000,
                          TOMType = "signed",
                          corType="bicor",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "/home/famgarcia/TFM/TOM_DS_COV",
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations


# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(SAHA_Proembryos = ifelse(grepl('SAHA', Treatment), 1, 0))

traits<-as.data.frame(traits[,c(5)])
rownames(traits)<-rownames(colData)
colnames(traits) <-"SAHA_Treatment"


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
textMatrix = paste( signif( module.trait.corr[,], 2 ), "  ;  (" , # Pasting corr to pvalues
                    signif( module.trait.corr.pvals[,], 1), 
                    ifelse(module.trait.corr.pvals[,]<0.05, "*",""), 
                    ifelse(module.trait.corr.pvals[,]<0.01, "*",""), 
                    ifelse(module.trait.corr.pvals[,]<0.001, "*",""),")", sep="")
dim(textMatrix)=dim(module.trait.corr[,])
textMatrix=as.matrix(textMatrix)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


sign_modules = rownames(module.trait.corr.pvals[rowSums(module.trait.corr.pvals < 0.05)>0,])
sign_colors = gsub("ME", "",sign_modules)


x11(type="cairo")
png("/home/famgarcia/Escritorio/modcorr.png", width=18, height=25, units="cm", res=300)
par(mar=c(7,7,0.5,2))
labeledHeatmap(Matrix=module.trait.corr[,],
               xLabels=names(heatmap.data)[28],
               yLabels=names(heatmap.data)[1:27],
               ySymbols=gsub("ME","",names(heatmap.data)[1:27]),
               colorLabels = FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins = FALSE,
               cex.text=1,
               zlim=c(-1,1)) 
dev.off()
png("/home/famgarcia/Escritorio/modcorr_SAHA.png", width=18, height=25, units="cm", res=300)
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[28],
             y = names(heatmap.data)[1:27],
             col = c("blue1", "skyblue", "white", "pink", "red"))
dev.off()


# yellow red


#  GENE SIGNIFICANCE MODULE MEMBERSHIP ----------------------
# or each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profil


# 6B. Intramodular analysis: Identifying driver genes ---------------



# I. MODULE MEMBERSHIP. Calculate the module membership and the associated p-values-------------------------

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

# Gene significance analysis - not for traits defined but for modules - module eigengenes.

# II. GENE SIGNIFICANCE  -------------------------------------------
# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.
# Calculate the gene significance and associated p-values

################### SAHA TREATMENT

gene.signf.corr <- cor(norm.counts, traits$SAHA_Treatment, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)


module.membership.measure.pvals[1:10,1:10]

# Plot for modules of interest - YELLOW

modNames = rownames(module.membership.measure)
module = "yellow"
column = match(paste("ME",module,sep=""), modNames);
associated_mod_genes <- module.gene.mapping %>% 
filter(`bwnet$colors` == 'yellow') %>% 
rownames()

moduleGenes <- which(colnames(module.membership.measure) %in% associated_mod_genes)
sizeGrWindow(7, 7);
par(mfrow = c(1,1));

verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
               abs(gene.signf.corr[moduleGenes, 1]),
               xlab = paste("Module Membership in", module, "module"),
               ylab = "Gene significance for SAHA treatment",
               main = paste("Module membership vs. gene significance\n"),
               cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# RED MODULE

modNames = rownames(module.membership.measure)
module = "red"
column = match(paste("ME",module,sep=""), modNames);
associated_mod_genes <- module.gene.mapping %>% 
filter(`bwnet$colors` == module) %>% 
rownames()

moduleGenes <- which(colnames(module.membership.measure) %in% associated_mod_genes)

verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                   abs(gene.signf.corr[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for SAHA treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

  

colnames(gene.signf.corr.pvals)=c("pvalue")
significant_genes_treatment = subset(gene.signf.corr.pvals, pvalue<=0.05)
write.table(rownames(significant_genes_treatment), "hub_SAHA_ONLYPE_NETWORK.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
  
  

# Annotate results


# NCBI GENE DESCRIPTIONS

º# Expr log2fc among SAHA treated PE - Control PE



hub_SAHA = as.data.frame(significant_genes_treatment)
hub_SAHA$entrezgene_id = rownames(hub_SAHA)
putative_hubgenes_involved_treatment = merge(hub_SAHA,annot_ncbi, by='entrezgene_id')
expr_phgenes_dev_treatment = merge(putative_hubgenes_involved_treatment, expr_HDAC, by="external_gene_name")
expr_phgenes_dev_treatment = expr_phgenes_dev_treatment[order(expr_phgenes_dev_treatment$pvalue.x),]

write.table(expr_phgenes_dev_treatment, "hub_Stress_SAHA_OnlyPE_Network.tsv", sep="\t", row.names=FALSE, quote=FALSE)


library(AnnotationHub)
# Retrieve record for Brassica napus OrgDb (NCBI gene ID annotations based on Brassica napus)
ah <- AnnotationHub()
AnnotationHub::query(ah, c("Brassica", "napus"))
# Bnapus OrgDb will be provided with the object AH107389
Bnapus <- ah[["AH107389"]]
  
ora_analysis_bp <- enrichGO(gene = hub_SAHA$entrezgene_id, 
                            universe = colnames(norm.counts), 
                            OrgDb = Bnapus,  
                            keyType = "ENTREZID",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE, 
                            pool = FALSE)



ora_analysis_bp <- enrichKEGG(gene = hub_SAHA$entrezgene_id, 
                            universe = colnames(norm.counts), 
                            organism = "bna",
                            keyType= 'ncbi-geneid',
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05)




# SAVE GENES IN SIGNIFICANT MODULES
# Match probes in the data set to the probe IDs in the annotation file 
probes = colnames(norm.counts)
allLLIDs = probes
moduleColors = bwnet$colors

#probes2annot = match(probes, annot_ncbi$entrezgene_id)
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$ALIAS[probes2annot];
# $ Choose interesting modules
intModules = c("red", "yellow")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs_OnlyPENetwork-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs_OnlyPENetwork_all-", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)



