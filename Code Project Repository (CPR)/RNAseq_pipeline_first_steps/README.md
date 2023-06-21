The contents of this folder include a series of utility scripts and general scripts to conduct the RNASeqPipeline with Paired End (PE) sequencing data

| Code  | Pipeline step  | Description |  Usage  |
|---|---|---|---|
| `QC_trimming.sh`  | Quality Checks & Trimming | Performs quality check reports with fastp/fastqc before and after a minimum filtering with fastqc  | Call in Unix shell `./QC_trimming.sh` + arg in User interaction  |
| `download_map_assemble.sh` | Mapping to reference genome + assemble  | Downloads _Brassica napus_ GFF3 gene model annotation files from Ensembl FTP endpoint and converts format with ACAT, gets exon and splicesites, performs and reports splice aware mapping with HiSAT2, creates bam indexes, assembles de novo with StringTie | `./download_map_assemble.sh <trimmed_fq_read_directory_path> <numsamples_per_condition>`  |
| `read_count.sh` |  Quantification of gene expression  |  Creates bam index (optional), reports and performs an union-exon based strategy for gene count using `download_map_assemble.sh` bam files | `./read_count.sh`   |
| `NormalizationPipeline` | Quantification of gene expression  |  Calculates gene length from GTF genome annotation file and calls R normalization script for TPM / FPKM (chosen in arguments) on all counts from a provided directory | `./NormalizationPipeline.sh` <countsdir> <GTFfile.gtf> <mode [TPM/FPKM]>`  |

Utility functions

| Code  | Description |
|---|---|---|
| `QC.sh`  | Calls fastqc to perform quality check of fastq files | 
| `mergeCounts.sh`  |  Merges counts before normalization  |  
| `NormalizeCountsR` |  Uses gene length file and merged counts csv to calculate TPM / FPKM  |  
