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

#remotes::install_github("kevinblighe/CorLevelPlot")
#install.packages("WGCNA", dependencies = TRUE)

allowWGCNAThreads()          # allow multi-threading (optional)
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





#############################################################


# Repeat for Development Stage covariate batch adjustment


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

soft_power <- 12
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
png("/home/famgarcia/Escritorio/modcorr.png", width=18, height=25, units="cm", res=300)
par(mar=c(7,7,0.5,2))
labeledHeatmap(Matrix=module.trait.corr[,],
               xLabels=t(as.data.frame(names(heatmap.data)[c(22,23,24)])),
               yLabels=names(heatmap.data)[1:21],
               ySymbols=gsub("ME","",names(heatmap.data)[1:21]),
               colorLabels = FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins = FALSE,
               cex.text=1,
               zlim=c(-1,1)) + geom_raster() + scale_fill_identity()
dev.off()

# black turquoise and greenyellow modules are of biological interest to compare PEHDAC / PECtrl sign modules
# brown, blue, midnight blue might be informative to analyze coexpression clusters in Vacuolated Microspores before embryogenesis



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

black_gid <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'black') %>% 
  rownames()

brown_gid <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()


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

# Plot for modules of interest?
modNames = rownames(module.membership.measure)
module = "brown"
column = match(paste("ME",module,sep=""), modNames);
associated_mod_genes <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()

moduleGenes <- which(colnames(module.membership.measure) %in% associated_mod_genes)
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                   abs(geneTraitSignificance[moduleGenes, column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)




# II. GENE SIGNIFICANCE  -------------------------------------------
# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.
# Calculate the gene significance and associated p-values

################### SAHA TREATMENT

gene.signf.corr <- cor(norm.counts, traits$SAHA_Proembryos, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)


verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                   abs(gene.signf.corr[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)




colnames(gene.signf.corr.pvals)=c("pvalue")
significant_genes_treatment = subset(gene.signf.corr.pvals, pvalue<=0.01)

################## DEVELOPMENTAL STAGE - Ctrl proembryos

gene.signf.corr <- cor(norm.counts, traits$Ctrl_Proembryos, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                   abs(gene.signf.corr[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)




gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)
colnames(gene.signf.corr.pvals)=c("pvalue")

significant_genes_stagedev_PE = subset(gene.signf.corr.pvals, pvalue<=0.01)
significant_genes_stagedev_PE$entrezgene_id = rownames(significant_genes_stagedev_PE)
significant_genes_stagedev_PE =significant_genes_stagedev_PE[order(significant_genes_stagedev_PE$pvalue),]

################## DEVELOPMENTAL STAGE - Ctrl proembryos

gene.signf.corr <- cor(norm.counts, traits$Vacuolated_Micrp, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


verboseScatterplot(abs(module.membership.measure[column, moduleGenes]),
                   abs(gene.signf.corr[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)
colnames(gene.signf.corr.pvals)=c("pvalue")

significant_genes_stagedev_VM = subset(gene.signf.corr.pvals, pvalue<=0.01)

write.table(rownames(significant_genes_stagedev_VM), "hub_VM.tsv", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(rownames(significant_genes_stagedev_PE), "hub_PEC.tsv", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(rownames(significant_genes_treatment), "hub_SAHA.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)



# Annotate results


significant_genes_stagedev$entrezgene_id = rownames(significant_genes_stagedev)
significant_genes_treatment$entrezgene_id = rownames(significant_genes_treatment)
rownames(significant_genes_stagedev) = NULL
rownames(significant_genes_treatment) = NULL

# NCBI GENE DESCRIPTIONS

º# Expr log2fc among SAHA treated PE - Control PE
hub_PEall = as.data.frame(significant_genes_stagedev_PE)
putative_hubgenes_involved_development = merge(hub_PEall,annot_ncbi, by='entrezgene_id')
putative_hubgenes_involved_development = putative_hubgenes_involved_development[order(putative_hubgenes_involved_development$pvalue),]


expr_phgenes_dev = merge(putative_hubgenes_involved_development, expr_Stress, by="external_gene_name")
expr_phgenes_dev = expr_phgenes_dev[order(expr_phgenes_dev$pvalue.x),]

write.table(expr_phgenes_dev, "hub_Stress_PEall.tsv", sep="\t", row.names=FALSE, quote=FALSE)




hub_VM = as.data.frame(significant_genes_stagedev_VM)
hub_VM$entrezgene_id = rownames(hub_VM)
putative_hubgenes_involved_development = merge(hub_VM,annot_ncbi, by='entrezgene_id')
putative_hubgenes_involved_development = putative_hubgenes_involved_development[order(putative_hubgenes_involved_development$pvalue),]


expr_phgenes_dev = merge(putative_hubgenes_involved_development, expr_Stress, by="external_gene_name")
expr_phgenes_dev = expr_phgenes_dev[order(expr_phgenes_dev$pvalue.x),]

write.table(expr_phgenes_dev, "hub_Stress_VM.tsv", sep="\t", row.names=FALSE, quote=FALSE)




hub_SAHA = as.data.frame(significant_genes_treatment)
hub_SAHA$entrezgene_id = rownames(hub_SAHA)
putative_hubgenes_involved_treatment = merge(hub_SAHA,annot_ncbi, by='entrezgene_id')
expr_phgenes_dev_treatment = merge(putative_hubgenes_involved_treatment, expr_HDAC, by="external_gene_name")
expr_phgenes_dev_treatment = expr_phgenes_dev_treatment[order(expr_phgenes_dev_treatment$pvalue.x),]

write.table(expr_phgenes_dev_treatment, "hub_Stress_SAHA.tsv", sep="\t", row.names=FALSE, quote=FALSE)


###############################################################################################
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
modules = c("black", "greenyellow", "brown", "turquoise", "midnightblue","blue");
modules = sign_colors

# Select module probes
probes = colnames(norm.counts)
inModule = is.finite(match(bwnet$colors, modules));
modProbes = probes[inModule];
modGenes = annot_ncbi$ALIAS[match(modProbes, annot_ncbi$entrezgene_id)];
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




####################################################################################################

# VISUALIZATION
#=====================================================================================
#
#  VISUALIZE NETWORK OF EIGENGENES - ADJACENCY + TOPOLOGICAL OVERLAP PLOT
#
#=====================================================================================

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, bwnet$dendrograms, bwnet$colors, main = "Network heatmap plot, all genes")


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


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
intModules = c("black", "greenyellow","turquoise")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)


############## WCGNA does not support GOEnrichment for Brassica napus, this will be undertaken with clusterProfiler







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



