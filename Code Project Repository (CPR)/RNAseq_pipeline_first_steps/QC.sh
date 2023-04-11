#!/bin/bash

# Description: Helper script to perform first Quality Check (QC) of reads from raw read sequence files with fastqc and multiqc
# 
# Author : Natalia García Sánchez
# Date : 23/02/2023


check_Idir_exists(){

# checks if input directory exists

if [ -d "$1" ] 
	then 
		echo "Found input directory"
	else 
		echo "$1 directory not found. Exiting..."
		exit 1
	fi
		}


# MAIN

main(){

GR='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'


echo -e "\n${GR}Enter the directory where your read files can be found${NC}"

read path_to_raw_reads
check_Idir_exists "$path_to_raw_reads"


echo -e "Files will be found in ${RED} QC_output ${NC} output directory in working directory. \nCreating output directory..."

mkdir QC_output
mkdir QC_output/logs


echo -e "\nPerforming QC and developing joint report..."

for file in "$path_to_raw_reads"/*; do fastqc "$file" -o QC_output
multiqc QC_output/*.zip 2>QC_output/logs
}

main
