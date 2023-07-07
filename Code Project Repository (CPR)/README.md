# Code Project Repository (CPR)
---
## Description preface

---
## Overview of tools used
- QC : MultiQC and fastqc
- Trimming : fastp
- Mapping : HiSAT2 using BW algorithm with graph based indexing(open source binary available in following [link with install instructions](http://daehwankimlab.github.io/hisat2/download/)) 
- Sorting bam : samtools
- Quantification of expression : HTSeq - HTSeqCount, NormalizationPipeline.sh (R source)
- Differential expression analysis : DeSeq2
- Gene enrichment : clusterProfiler, bioMaRt
- Batch correction : Combat-seq
- Coexpression Analysis : WGCNA


---
## General Information

| Document Summary||
|-----------------|--------|
| Project| Bnapus_ME TFM |
| Date | - |
|Summary | A collection of the methods to perform RNAseq analysis pipeline|

### List of Authors

| Author name| Email|
|-----------------:|-----------|
|Natalia García Sánchez | natalia.garcia.sanchez@alumnos.upm.es |


## Code summary

| Folder  | Code's aim |
| - |-|
| AnnotationRetrieval | Retrieval of annotations of miRNA in genomic gff, ATG autophagy related genes, Oxidative stress genes |
| RNAseq_pipeline_first_steps| Scripts relative to trimming, QC, mapping to a reference genome, read count and normalization |
| DifferentialExpression | Scripts relative to applying statistical inference and methods for inference of Differential Expression  |
| GeneEnrichment | Applying algorithms to get a measure of statistical significance of overrepresentation of genes belonging to a curated biological category through (ORA, GSEA, GSVA) |
| Coexpression Network WGCNA | Scripts relative to inference of coexpression among genes and conversion into a regulatory network separating genes by coexpression modules. Batch effect correction of expression data prior to this process is also included in this folder| 
| List generation | A collection of hormone pathways /Epigenetic mechanism terms used to filter genes |


