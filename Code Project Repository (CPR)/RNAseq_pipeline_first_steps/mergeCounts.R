#!/usr/bin/env Rscript

# Author : Natalia García Sánchez
# Date : 3/03/2023
# Description : script to do outermerge on files to get summary of HTSeqCounts
# Also shows read count statistics
# Usage : mergeCounts.R countdir
# ---

# Print usage message
writeLines(sprintf("Using mergeCounts. script usage: \n mergeCounts.R countdirectory \n For the script to recognize count files, they require the .counts extension \n ---"))

args = commandArgs(trailingOnly = TRUE)

# Create an object with the directory containing HTseq counts:
directory <- args[1]

#list.files(directory)
sampleFiles <- list.files(directory, pattern = ".*counts")
sampleFiles <- strsplit(sampleFiles, " ")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- read.table(paste(directory,sampleFiles[1], sep ="/", collapse=NULL))
colnames(dataset) <- c("gene_id", strsplit(as.character(sampleFiles[[1]]),".counts"))

#had to specify columns to get rid of the total column
for (i in 2:length(sampleFiles)){
  temp_data <- read.table(paste(directory,sampleFiles[i], sep ="/", collapse=NULL)) #each file will be read in
  colnames(temp_data) <- c("gene_id", strsplit(as.character(sampleFiles[[i]]),".counts"))

  dataset <- merge(x = dataset, y = temp_data, by = "gene_id", all = TRUE) #for each iteration, bind the new data to the building dataset
}

# Print statistics about count per file
print('Merged HTSeqCount stats')
print(dataset[1:5,])

# Only get counts
countsTable<-dataset[6:length(dataset[,1]),]
rownames(countsTable)<-countsTable$gene_id
countsTable$gene_id<-NULL

# Write to csv
print('Creating output directory \"merged\"...')
dir.create(paste(directory,"merged/", sep = "/"))
print('Printing merged counts onto \"merged.counts\" directory nested in count dir')
write.table(countsTable, file=paste(directory,"merged/merged.counts", sep ="/", collapse=NULL), sep=" ", , dec=".", row.names=TRUE, col.names=TRUE)



