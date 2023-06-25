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

write.table(data_to_norm,"/home/famgarcia/Escritorio/VST_norm_counts.csv")



################################ NO BATCH EFFECT CORRECTION

# Variance Stabilizing Transformation with DESeq2

data_to_norm <- read.csv("/home/famgarcia/Escritorio/genes_all_experiments.csv", sep=",")
rownames(data_to_norm) <- data_to_norm$gene_id
data_to_norm$gene_id <- NULL
data_to_norm <- data_to_norm[rowSums(data_to_norm)>10,]
data_to_norm <- data_to_norm[,colnames(data_to_norm) != "PErep3"]
# create dds to perform a Variance Stabilizing Transformation Normalization from the counts
dds_tonorm <- DESeqDataSetFromMatrix(countData = data_to_norm,
                                     colData = colData,
                                     design = ~ 1) # not spcifying model


vst_c<-vst(dds_tonorm)
write.table(assay(vst_c),"/home/famgarcia/Escritorio/VST_norm_counts_no_bce.csv")


pca <- prcomp(t(assay(vst_c)))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)
pca.dat$Treatment = c("Heat Stress","Heat Stress","Heat Stress","Heat Stress + SAHA","Heat Stress + SAHA","Heat Stress + SAHA","Heat Stress","Heat Stress","Heat Stress","Heat Stress","Heat Stress")
pca.dat$Dev_Stage = c("Proembryo","Proembryo","Proembryo","Proembryo","Proembryo","Proembryo","Vacuolated Microspore","Vacuolated Microspore","Vacuolated Microspore","Proembryo","Proembryo")

png("/home/famgarcia/VST_nbce.png", width=20, height=10, units="cm", res=300)
ggplot(pca.dat, aes(PC1, PC2, color=Dev_Stage, shape=Treatment)) +
  geom_point() +
  geom_text_repel(label = rownames(pca.dat), vjust=1.5, size=2) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' % Explained Variance'), size=2,
       y = paste0('PC2: ', pca.var.percent[2], ' % Explained Variance'), size=2 ) + 
  scale_x_continuous(limits = c(-250,500))+
  scale_y_continuous(limits = c(-200,150)) + theme(text=element_text(size=10))
dev.off()



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

write.table(data_to_norm,"/home/famgarcia/Escritorio/VST_norm_counts.csv")

