#!/bin/bash

# Code for retrieval of types of RNA in old annotation gtf file

# Author: Natalia García Sánchez
# Date 21/02/2023
# v.0.1
# ---

mv ~/TFM/owndata

# Retrieval of RNA types in gff genome file

cat genome.gtf  | grep -oP "gene_biotype\s\".+RNA\";\s" | cut -d" " -f2 | sort -u |uniq >RNA_types.txt

# Creation of Annotation Final File

touch ncRNA_annotations

# Detect if string has mirna/antisense/siRNA/piRNA/miscRNA as RNA type 

for line in $(cat RNA_types.txt)
do
	if [[ $line =~ "antisense" ]] || [[ $line =~ "misc" ]] || [[ $line =~ [spm]iRNA ]]
	then
		grep "$line" genome.gtf >>ncRNA_annotations
	fi

done

num_ann_id=$(cat ncRNA_annotations | grep -oP "gene_biotype\s\"\w+\"[;\s$]" | wc -l)
num_ann_type=$(cat ncRNA_annotations | grep -oP "gene_id\s\"\w+\"[;\s$]" | wc -l)

if [ $num_ann_id -eq $num_ann_type ];
then

	## Getting gene id and gene biotype lines from new annotation file
	cat ncRNA_annotations | grep -oP "gene_id\s\"\w+\"[;\s$]" | cut -d" " -f2 | sed 's/[\";]//g' >ncRNA_gID_tmp
	cat ncRNA_annotations | grep -oP "gene_biotype\s\"\w+\"[;\s$]" | cut -d" " -f2 | sed 's/[\";]//g' >ncRNA_type_tmp


	## Pasting both files (all annotations will have a g_id and g_biotype) 
	##  and put them into a "non coding RNA genes of interest"  ncRNA_genes_interest tsv file

	### print header for file
	echo -e "gene_id\tgene_biotype" >ncRNA_genes_interest.tsv
	paste ncRNA_gID_tmp ncRNA_type_tmp | sort -n | uniq >>ncRNA_genes_interest.tsv


	## Remove tmp files
	rm ncRNA_gID_tmp
	rm ncRNA_type_tmp

else
	# Only take gene id, query them in downstream analysis if gene biotype is needed

	## Getting gene id and gene biotype lines from new annotation file

	### print header for file
	echo -e "gene_id" >ncRNA_genes_interest.tsv
	cat ncRNA_annotations | grep -oP "gene_id\s\"\w+\"[;\s$]" | cut -d" " -f2 | sed 's/[\";]//g' >>ncRNA_genes_interest.tsv
	
fi

# Remove temporal annotation and RNA types files
rm ncRNA_annotations
rm RNA_types.txt