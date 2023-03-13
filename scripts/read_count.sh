#!/bin/bash

# Description: Creating count files of undirected sequenced library
# Author : Natalia García Sánchez
# Date : 24/02/2023
# ---

# RETRIEVING AND CHECKING INPUT DATA

# Asking for dir path with input bam files and checking if dir is empty

echo "Enter directory path to bam files: \n>"
path_to_bam_files=read()

if [ -d "$path_to_bam_files" ]
then
	echo "Bam file directory exists"
	
	if [ $(ls *.bam | wc -l) -gt 0 ]; 
	then 
		echo "Bam files found in directory. Proceeding to read count";
	
	else
		echo "No bam files found. \nExiting..."
		exit(1)

	fi

else
	echo "Bam file directory is empty. Try again with a proper dir. \nExiting..."
	exit(1)

fi


# Ask if GTF annotation file is present. Otherwise Download latest annotation in GTF format

echo "Do you have an annotation file available (y/n)"
read response
case $response in
	Y|y|[Yy]es)
		echo "Provide GTF filepath"
		read path_to_genome
	;;
	N|n|[Nn]o)
	    echo "Latest annotation GTF file will be downloaded"
	  	echo "Downloading B napus latest annnotation GFF file"
		wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/3708/102/GCF_020379485.1_Da-Ae/GCF_020379485.1_Da-Ae_genomic.gtf.gz
		gzip -d GCF_020379485.1_Da-Ae_genomic.gtf.gz
		mv GCF_020379485.1_Da-Ae_genomic.gtf.gz genome.gtf
		path_to_genome=$(pwd)"/genome.gtf"
	;;
	  *)
	    echo "Irrelevant response. Exitting"
	    exit(0)
	;;
	esac


# SORTING AND MAKING REPORTS OF BAM FILES

# Reading bam file with samtools

conda source activate
conda activate TFM
mkdir sorted_bam
mkdir flagstat_reports

for file in $path_to_bam_files/*.bam; do 
	# Sorting by position
	samtools sort $file >sorted_bam/$file
done

cd sorted_bam

#Performing flagstat 
for file in *.bam; do samtools flagstat $file >../flagstat_reports/${file%.bam}.flagstat_report.txt; done

cd ..

# take first column from flagstat                                               
cat $(ls flagstat_reports/*flagstat_report*.txt | head -n1) |  sed 's/\s0\s/\s0\t/g' | cut -d"$(echo -e "\t")" -f2 | sed 's/(.*)//g' >flagstat_reports/flagstat_report

# Join all flagstat reports into a single document and remove the single files  

for file in flagstat_reports/*.flagstat_report.txt; do 

	# Get only information about reads into tmp_file
	cat $file |  sed 's/\s0\s/\s0\t/g' | cut -d"$(echo -e "\t")" -f1 >flagstat_reports/tmp_file

	join flagstat_reports/flagstat_report flagstat_reports/tmp_file; 

	done
rm flagstat_reports/*.flagstat_report.txt


# Adding header
printf "MAPPING\t" >flagstat_reports/flagstat_report
echo "$(ls flagstat_reports/*flagstat_report*.txt)" | tr "\n" "\t" >>flagstat_reports/flagstat_report



# INDEXING FILES PRE-READ COUNT

for file in $path_to_bam_files/*.bam; do samtools index $file; done

conda deactivate
conda source activate
conda activate __

# --- Performing htseq count for genes

for file in $path_to_bam_files/*.bam; 
do 
	htseq-count -i gene_id -s no -r pos -f bam $file $path_to_genome -m union -t exon >${file%.bam}".counts"; 

done


# --- Joining all counts

# Printing first part of the header
printf "Gene\t" >merged.counts.tsv

# Getting the rest of the header and sample order
for file in *.counts; 
do 
	printf %s "$file" >>merged.counts.tsv
	printf "\t" >>merged.counts.tsv
done
fileorder=($(ls *.counts))

# joining the contents of the count files one by one in order
join -t "$(echo -e "\t")" ${fileorder[0]} ${fileorder[1]} | join -t
"$(echo -e "\t")" - ${fileorder[2]}| join -t "$(echo -e "\t")" -
${fileorder[3]} >>merged.counts.tsv