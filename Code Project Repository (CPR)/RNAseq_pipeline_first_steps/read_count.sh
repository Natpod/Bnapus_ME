#!/bin/bash

# Description: Sorting and reporting bam, creating count files of undirected sequenced library
# Author : Natalia García Sánchez
# Date : 24/02/2023
# ---

# Bash output color parameter for user interaction
GREEN='\033[0;32m'
NC='\033[0m'

echo -e "\n${GREEN} ~~~~~~~~~~~~ read_count.sh ~~~~~~~~~~~~\n${NC}"
echo -e "A script for read count using HISAT read alignment output BAM files"
echo -e "\n\n ${GREEN} USAGE : ${NC} bash download_map_assemble.sh ${GREEN} \n\tpath_to_bam and genome files are provided in user input.${NC}"

help(){
echo -e echo -e "\n${GREEN} ~~~~~~~~~~~~ read_count.sh ~~~~~~~~~~~~\n${NC} USAGE : path_to_bam and genome files are provided in user input\n bam directory required can optionally contain custom indexes."
exit 1
}

# Loading necessary modules
echo -e "${GREEN}\n------------------------------------\nLoading SAMTOOLS, AGAT, PERL, GCC, HTSeq dependencies\n------------------------------------${NC}"
module load Python/3.8.6-GCCcore-10.2.0
module load GCCcore/10.2.0
module load Perl/5.32.0-GCCcore-10.2.0-minimal
module load gzip/1.10-GCCcore-8.2.0
module load HTSeq/0.11.1-foss-2019a
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load AGAT/0.9.2-GCC-11.2.0

# RETRIEVING AND CHECKING INPUT DATA

# Asking for dir path with input bam files and checking if dir is empty

echo -e "${GREEN}\n------------------------------------\n1) Checking bam and genome files provided in user input \n-------------------${NC}"
echo -e "${GREEN}Enter directory path to bam and bam.bai (index) files: \n>${NC}"
read path_to_bam_files

if [ -d "$path_to_bam_files" ]
then
	echo "Bam file directory exists"
	
	if [[ $(ls "$path_to_bam_files"/*.bam | wc -l) -gt 0 ]] 
	then 
		echo "Bam files found in directory. Proceeding to checking index files."
		
		if [[ $(ls "$path_to_bam_files"/*.bam.bai | wc -l) -gt 0 ]]
		then
			echo -e "\n-Found index files in directory. Proceeding to read count..."
		else 
			echo -e "\n-No index files were found. Proceeding to indexing files and read count"
			
			for file in "$path_to_bam_files"/*.bam
			do
				ifilename=$(basename "$file")
				echo "Indexing $ifilename"	 
				samtools index "$file" >"$ifilename".bai
			done
		fi

	else
		echo -e "No bam files found. \nExiting..."
		help

	fi

else
	echo -e "Bam file directory is empty. Try again with a proper dir. \nExiting..."
	help

fi


# Ask if GTF annotation file is present. Otherwise Download latest annotation in GTF format

echo -e "${GREEN}Do you have an annotation file available (y/n)${NC}"
read response
case "$response" in
	Y|y|[Yy]es)
		echo -e "${GREEN}Provide GTF filepath\n>${NC}"
		read path_to_genome
	;;
	N|n|[Nn]o)
	    	echo -e "${GREEN}\n------------------------------------\nDownloading Brassica napus Ensembl Plants GFF3 file annotation and converting into GTF file \n--------------------------------${NC}"
                wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/brassica_napus/Brassica_napus.AST_PRJEB5043_v1.56.gff3.gz
		gzip -d Brassica_napus.AST_PRJEB5043_v1.56.gff3.gz
		echo "Converting GFF3 file to GTF3"
		agat_convert_sp_gff2gtf.pl --gff Brassica_napus.AST_PRJEB5043_v1.56.gff3 -o genome.gtf 2>tmp.conversion.log

	;;
	  *)
	    echo "Irrelevant response. Exiting"
	    exit 0
	;;
esac

if [[ ! -d read_count_output ]]
then
mkdir read_count_output
fi

# --- Performing htseq count for genes
echo -e "\n${GREEN}\n------------------------------------\n2) Performing read count\n------------------------------------\n${NC}"
for file in "$path_to_bam_files"/*.bam 
do
    cfilename=$(basename "$file")
    srun -p long -t 1:00:00 --mem 128GB htseq-count -f bam -i gene_id -s no -r pos -m union -t exon "$file" "$path_to_genome" >read_count_output/"${cfilename%.bam}".counts 

done





