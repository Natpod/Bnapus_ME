# Coexpression network analysis

As an additional informative way to get insight from the transcriptomic data, a Weighted Gene Coexpresion Analysis (WGCNA) with batch-adjusted normalized high signal data from both experiments. The purpose of the code is cluster genes by coexpression in an unsupervised way, so that we extrapolate coexpression top putative corregulation in a biological network comprised of modules, and associate the former modules to a trait in the expression data (Developmental Stage or HDACi treatment)

Code in this section

- batch_correction.R : analyses and use of Combat-Seq to perform batch correction adjusting for a biological replicate.
- WGCNA.R : performs cleaning preprocessing steps from adjusted Combat-seq raw count, normalization and weighted conetwork analysis --> selection of hub genes and modules related to VM/PE/SAHA trait