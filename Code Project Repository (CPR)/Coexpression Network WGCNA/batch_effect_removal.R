#!/usr/bin/env Rscript

# Natalia García Sánchez
# Description : Removal of batch effects between two experiments-3 biological conditions(?) with combat


# Install packages
#devtools::install_github("zhangyuqing/sva-devel")
#install.packages("dendextend")

library("dendextend")
library("sva")
library("RUVSeq")
library("DESeq2")
library("ggplot2")
library("dplyr")
library("reshape2")
library("gridExtra")
library("scales")
library("ggpubr")
library("BatchQC")


############################################ FINAL ADJUSTMENT

# Without specifying HDACi - Groups
all_expr_data = read.csv("D:\\other\\genes_all_experiments.csv")

rownames(all_expr_data) <- all_expr_data$gene_id
all_expr_data$gene_id<- NULL

# Filter taking out biological replica discarded in Stress analysis - outlier

all_expr_data<-all_expr_data[,colnames(all_expr_data)!= "PErep3"]
counts<- as.matrix(all_expr_data)

# remove genes with only 0 counts in the subset & in any batch
keep1 <- apply(counts[, batch==1],1,function(x){!all(x==0)})
keep2 <- apply(counts[, batch==2],1,function(x){!all(x==0)})
counts <- counts[keep1 & keep2, ]


### Parameters for batch removal

# Different experiments represented in each batch (protocol for sequencing and getting counts is the same, technical variance)
batch <- c(rep(1, 6), rep(2, 5))

# 0 for vacuolated microspore, 1 for proembryos
Developmental_Stage <- c(1,1,1,1,1,1,0,0,0,1,1)



# Adjust for batch effects

adjusted_counts <- ComBat_seq(counts, batch=batch, group=Developmental_Stage)

# Normalize for library size
cts_norm <- apply(counts, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(adjusted_counts, 2, function(x){x/sum(x)})


# Exploratory analysis of normalized counts

# PCA

col_data <- data.frame(Batch=factor(batch), Group=factor(c("PE","PE","PE","PE_HDACi","PE_HDACi","PE_HDACi","VM","VM","VM","PE","PE")))
rownames(col_data) <- colnames(counts)

seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)

pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group), theme=theme_bw()) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Unadjusted PCA")+ theme(plot.title = element_text(size=12))

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group), theme=theme_bw()) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="ComBat-Seq PCA") + theme(plot.title = element_text(size=12))

plt_PCA_full <- ggarrange(plt,plt_adj, ncol=1, nrow=2, common.legend=TRUE, legend="left")


write.table(pca_obj_adj$data, "D:\\PCA_N_adjusted_batch_effects.tsv", sep="\t", quote=FALSE)
write.table(pca_obj$data, "D:\\PCA_N_with_batch_effects.tsv", sep="\t", quote=FALSE)


# Clustering analysis : hierarchical clustering -average linkage

dend_be <- as.dendrogram(hclust(dist(t(cts_norm)), method = "average"))
# let's add some color:
colors_to_use <- ifelse(as.numeric(batch)-1, "#2596be","#ea906c")
# But sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend_be)]
# Now we can use them
labels_colors(dend_be) <- colors_to_use
# Now each state has a color
dend_be <- dend_be %>% 
  set("branches_k_color", values<-c("#8d7b68","#159895","#b46060"), k = 3)  %>%
  set("labels_cex", 0.58)

par(mar=c(0, 3, 2, 2))
d1 <- ggplot(dend_be, theme = theme_minimal()) + 
  scale_y_continuous(limits=c(-0.0072,0.047)) +
  labs(title="Unadjusted Clustering")+
  theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.line.y.left = element_line(),
  axis.ticks.y.left= element_line(),
  plot.title = element_text(size=12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
)

dend_adj <- as.dendrogram(hclust(dist(t(cts_adj_norm)), method = "average"))
colors_to_use <- ifelse(as.numeric(batch)-1, "#2596be","#ea906c")
colors_to_use <- colors_to_use[order.dendrogram(dend_adj)]
labels_colors(dend_adj) <- colors_to_use
dend_adj <- dend_adj %>% 
  set("branches_k_color", values<-c("#8d7b68","#159895"), k = 2) %>%
  set("labels_cex", 0.58)

d2 <- ggplot(dend_adj, theme = theme_minimal()) + 
  scale_y_continuous(limits=c(-0.013,0.09), breaks = c(0.04,0.03,0.02,0.01,0.00)) +
  labs(title="Combat-Seq Clustering")+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y.left = element_line(),
    axis.ticks.y.left= element_line(),
    plot.title = element_text(size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )


plt_d_full <- ggarrange(d1,d2, ncol=1, nrow=2)+ labs(y="Clustering analysis")

varexp_full <- list(
  unadjusted=batchqc_explained_variation(cpm(counts, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combatseq=batchqc_explained_variation(cpm(adjusted_counts, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation
)

varexp_full_df <- melt(varexp_full)
varexp_full_df$L1 <- factor(varexp_full_df$L1, levels=c("unadjusted","combatseq"))
varexp_full_df$L1 <- plyr::revalue(varexp_full_df$L1, c("unadjusted"="Unadjusted", "combatseq"="ComBat-Seq"))
varexp_full_df$Var2 <- plyr::revalue(varexp_full_df$Var2, c("Full (Condition+Batch)"="Condition+Batch"))

plt_varexp_full <- ggplot(varexp_full_df, aes(x=Var2, y=value), theme=theme_bw()) +
  geom_boxplot() +
  facet_wrap(~L1, nrow=4, ncol=1) +
  labs(y="Explained variation") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

par(mar=c(0, 4, 4, 2))
final_plt <- ggarrange(plt_PCA_full,plt_d_full, plt_varexp_full, ncol=3, widths=c(0.70, 0.5, 0.35))
final_plt

ggsave(file="D:\\Batch_effect_adjustment.svg", plot=final_plt, width=10, height=8)

write.table(as.data.frame(adjusted_counts), "D:\\genes_all_experiments_batch_correction.csv",sep=",", row.names=TRUE, quote=FALSE, col.names=TRUE)


ggplot(melt(as.data.frame(adjusted_counts)), aes(x=variable, y=value), theme=theme_bw()) +geom_boxplot()

ggplot(melt(as.data.frame(cts_adj_norm)), aes(x=variable, y=value), theme=theme_bw()) +geom_boxplot()


rm(cts_adj_norm)
rm(cts_norm)
rm(adjusted_counts)
rm(all_expr_data)
rm(counts)
rm(final_plt)

################ TPM violinplots



df_lgTPM <- melt(as.data.frame(lTPM_before_bc))
colnames(df_lgTPM)<- c("sample","values")
df_lgTPM$Group <- rep(" ",length(df_lgTPM[,1]))

index_PEs<-which(df_lgTPM$sample %in% c("PErep1", "PErep2", "Crep1", "Crep2", "Crep3"))
index_PESAHAs<-which(df_lgTPM$sample %in% c("Trep1", "Trep2", "Trep3"))
index_VMs<-which(df_lgTPM$sample %in% c("VMrep1", "VMrep2", "VMrep3"))

df_lgTPM[index_PEs, 3] <- "Ctrl Proembryos"
df_lgTPM[index_PESAHAs, 3] <- "SAHA Proembryos"
df_lgTPM[index_VMs, 3] <- "Vacuolated Microspores"



vp1 <- ggplot(df_lgTPM, aes(sample, values, fill=Group, color=Group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.5) + 
  geom_boxplot(width=0.07, outlier.size=0.05, outlier.color = '#BCBCBC') +
  ylab(expression(log[2]("TPM + 1"))) +
  labs(title="a2")+
  theme_classic() + 
  scale_fill_manual(values=c('#BCE4D8','#E6A492','tan')) +
  scale_color_manual(values=c('#44AA8A','lightcoral','#BB8E51')) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 7, color='#6B6B6B'),
        axis.title.x = element_blank(),
       legend.position="none"
       ) 



df_lgTPM_bce <- melt(as.data.frame(lTPMval_bc))
colnames(df_lgTPM_bce)<- c("sample","values")
df_lgTPM_bce$Group <- rep(" ",length(df_lgTPM_bce[,1]))

index_PEs<-which(df_lgTPM_bce$sample %in% c("PErep1", "PErep2", "Crep1", "Crep2", "Crep3"))
index_PESAHAs<-which(df_lgTPM_bce$sample %in% c("Trep1", "Trep2", "Trep3"))
index_VMs<-which(df_lgTPM_bce$sample %in% c("VMrep1", "VMrep2", "VMrep3"))

df_lgTPM_bce[index_PEs, 3] <- "Ctrl Proembryos"
df_lgTPM_bce[index_PESAHAs, 3] <- "SAHA Proembryos"
df_lgTPM_bce[index_VMs, 3] <- "Vacuolated Microspores"

vp2<- ggplot(df_lgTPM_bce, aes(sample, values, fill=Group, color=Group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.5) + 
  geom_boxplot(width=0.07, outlier.size=0.05, outlier.color = '#BCBCBC') +
  ylab(expression(log[2]("TPM + 1"))) +
  labs(title="b2")+
  theme_classic() + 
  scale_fill_manual(values=c('#BCE4D8','#E6A492','tan')) +
  scale_color_manual(values=c('#44AA8A','lightcoral','#BB8E51')) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 7, color='#6B6B6B'),
        axis.title.x = element_blank())
        # legend.position="none",
        # legend.background = element_rect(fill="#F2F2F2",
        #                                  size=0.5),
        # legend.text = element_text(size=7),
        # legend.title = element_text(size=8)


plt_varexp_full <- ggarrange(vp1,vp2, ncol=1, nrow=2, legend = FALSE)


par(mar=c(0, 4, 4, 2))
final_plt <- ggarrange(plt_d_full, plt_varexp_full,plt_PCA_full, ncol=3, widths=c(0.4, 0.45, 0.70))
final_plt

ggsave(file="D:\\Batch_effect_adjustment_no_covariates_woPErep.svg", plot=final_plt, width=10, height=6.5)

