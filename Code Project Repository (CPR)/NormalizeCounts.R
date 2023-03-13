#!/usr/bin/env Rscript

# Author : Natalia García Sánchez
# Date : 3/03/2023
# Description : script to do TPM calculation
# Usage : NormalizeCounts.R merged.counts genelengthfile mode
# ---


# Print usage message
writeLines(sprintf("\n NormalizeCounts.R \n script usage: \n TPMcalc.R mergedcountsfile \n ---"))


args = commandArgs(trailingOnly = TRUE)

# Read arguments and import tables
filepath <- args[1]
genelengthpath <- args[2]
mode <- args[3]


mergedcounts <-read.table(filepath, header = TRUE)
genelength <-read.csv(genelengthpath, header = TRUE, sep=",")



# Define calc functions with normalization methods  FPKM/TPM

FPKM_calc <- function(mergedcount_dataset){
  
  # FPKM calculation
  
  # Count up the total reads in a sample and divide that number by 1,000,000 – to calculate scaling factor.
  # Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
  # Divide the RPM values by the length of the gene, in kilobases.
  
  countscol <- ncol(mergedcount_dataset)
  scalingfactor <- colSums(mergedcount_dataset)/(10**6)
  
  mergedcounts$gene_id <- rownames(mergedcounts)

  mgdataset <- merge(x = mergedcounts, y = genelength, by.x="gene_id", by.y = "gene_id") 
  normalizedcounts <- mgdataset[,2:(countscol+1)]
  
  for (i in 1:countscol){
    normalizedcounts[,i] <- (normalizedcounts[,i])/scalingfactor[i]
  }
  normalizedcounts <- (normalizedcounts*1000)/mgdataset$length
  
  normalizedcounts$gene_id <- mgdataset$gene_id
  return(normalizedcounts)
}

TPM_calc <- function(mergedcount_dataset){
  
  # TPM calculation
  
  # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
  # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
  # Divide the RPK values by the “per million” scaling factor. 
  
  countscol <- ncol(mergedcount_dataset)

  mergedcounts$gene_id <- rownames(mergedcounts)
  mgdataset <- merge(x = mergedcounts, y = genelength, by.x="gene_id", by.y = "gene_id")
  normalizedcounts <- mgdataset[,2:(countscol+1)]
  
  normalizedcounts <-  (normalizedcounts*1000)/mgdataset$length
  scalingfactor <- colSums(normalizedcounts)/(10**6)
    
  for (j in 1:countscol){
    normalizedcounts[,j] <- (normalizedcounts[,j])/scalingfactor[j]
  }
  normalizedcounts$gene_id <- mgdataset$gene_id
    return(normalizedcounts)
}

if (toupper(mode) == "FPKM"){
  sprintf("Using NormalizeCounts.R to calculate Fragments Per Kilobase of transcript per Million mapped reads")
  FPKM_table <- FPKM_calc(mergedcounts)
  write.table(FPKM_table, "FPKM_values.csv", sep=",", dec = ".")
}

if (toupper(mode) == "TPM"){
  sprintf("Using NormalizeCounts.R to calculate Transcripts Per Million")
  TPM_table <- TPM_calc(mergedcounts)
  write.table(TPM_table, "TPM_values.csv", sep=",", dec = ".")
}
