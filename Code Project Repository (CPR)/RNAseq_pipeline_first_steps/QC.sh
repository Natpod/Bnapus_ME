#!/bin/bash

# Description: Helper script to perform first Quality Check (QC) of reads from raw read sequence files with fastqc and multiqc
# Returns files in QC_ouput
# Author : Natalia García Sánchez
# Date : 23/02/2023



help(){
	echo -e "Usage of ${GR} QC.sh ${NC}:\n QC.sh dir_path_to_reads\nOutputs in ${RED} QC_output ${NC}"
}

check_Idir_exists(){

# checks if input directory exists

if [ -d "$1" ] 
	then 
		echo "Found input directory"
	else 
		echo "$1 directory not found. Exiting..."
		help
	fi
		}


# MAIN

main(){

echo -e "\n${GR}Checking the directory where your read files can be found${NC}"

check_Idir_exists "$1"


echo -e "Files will be found in ${RED} QC_output ${NC} output directory in working directory. \nCreating output directory..."

mkdir QC_output
mkdir QC_output/logs


echo -e "\nPerforming QC and developing joint report..."

for file in "$1"/*; do fastqc "$file" -o QC_output
multiqc QC_output/*.zip 2>QC_output/logs ; done 
}


# ----------

GR='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'


main $1
