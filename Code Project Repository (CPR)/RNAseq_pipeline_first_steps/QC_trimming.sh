#!/bin/bash

# Description: 1st stage of RNAseq pipeline script : 
# User input wrapper script for fastp jointly used with fastqc quality checks
# Author : Natalia García Sánchez
# Date : 23/02/2023
# ---
# 
# REQUIREMENTS
#
# Only for Illumina Novaseq platform transcriptome sequencing - paired end SBS
# 1 factor design experiment 
# (2 conditions with sequencing files either starting with C "Control", or T "Treatment")
#
#	- working with raw reads, 
#	- trimming removing
#		* adapters in user input
#		* reads with PhredQ score<5
# 		* reads with more than 10 N bases
#		* reads under uniform length of 150bp
# 	- Performing Quality Control (QC) pre and post filtering


#---

ask_user_input_yn(){
echo -e "\033[0;31mDirectory is not empty\n\tCurrent directory will be overwritten. Want to continue? (y/n)${NC}"; 
				read response1

				case $response1 in
				  Y|y|[Yy]es)
					echo "Overwritting $1 directory"
				  ;;
				  N|n|[Nn]o)
				    echo "Leaving..."
				  	exit 0
				  ;;
				  *)
				    ask_user_input_yn
				  ;;
				esac
				
			}

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

check_Odir_exists(){

# checks if input directory exists
# if not empty, asks for user input to overwrite
# if it does not exist, creates directory

	if [ -d "$1" ] 
		then 
			echo "Found output directory"


		if [ -z "$(ls -A -- "$1")" ]; 
					then 
						ask_user_input_yn
					else 
						echo "Making directories"; 
						rm -r "$1"
						mkdir "$1"
					fi
				

		else
			echo "Making directories"; 
			mkdir "$1"
		fi
}

# MAIN

RED='\033[0;32m'
NC='\033[0m'

wdir=$(pwd)

# ENTER INPUT AND OUTPUT DIRECTORIES

#Asking for sequencing file user input
echo -e "\n-----------------------------------------\n1) Asking for sequencing file user input\n-----------------------------------------\n"

echo -e "\n${RED}Enter the directory where your raw read files can be found${NC}"
read path_to_raw_reads
check_Idir_exists "$path_to_raw_reads"

# Asking for number of samples in raw read files (per condition)
echo -e "\n${RED}Enter number of samples per condition${NC}"
read numbers

echo -e "\n${RED}Enter the filepath where your sequencing adapters in fasta format can be found${NC}"
read path_to_adapters
if [ -f "$path_to_adapters" ]; 
	then 
		echo "Adapter file exists"
	else echo "No adapter file found"
	fi

echo -e "\n-----------------------------------------\n2) Performing output directories control\n-----------------------------------------\n"


echo -e "\n${RED}Enter the directory where your want to decompress your raw read files${NC}"
read path_to_draw_reads
check_Odir_exists "$path_to_draw_reads"

echo -e "\n${RED}Enter the directory where you want your Quality Check (QC) files and trimming reports to be found${NC}"
read path_to_fq
check_Odir_exists "$path_to_fq"

echo -e "\n${RED}Enter the directory where your want to put your clean read files${NC}"
read path_to_clean_reads
check_Odir_exists "$path_to_clean_reads"


# Making fq reads output directory
if [ -d "$path_to_fq" ]; then rm -r "$path_to_fq"; fi
mkdir "$path_to_fq"
mkdir "$path_to_fq"/pre_filter_fastqc
mkdir "$path_to_fq"/post_filter_fastqc
mkdir "$path_to_fq"/trimmingreport

# QC AND TRIMMING


# Building an array that establishes the conditions of the experiment ('Control':'C', 'Treatment-DHDAC' : 'T')
array_conditions=('C' 'T')


echo -e "\n-----------------------------------------\n3) Performing first QC and trimming\n-----------------------------------------\n"



# For every F/R pair of read .fq files in the paired sequencing, sample and condition: 
# 1) Files are decompressed
# 2) FQ reports are made for each file

for (( samplenum=1; samplenum<="$numbers"; samplenum++ ));
	do
		echo "Iterating over sample num $samplenum"

		for condition in "${array_conditions[@]}";
		do 
			

			# Find Forward and Reverse sequence files names, put them in a temporary array

			filebatch=( )
			ls -A "$path_to_raw_reads"/"$condition"rep"$samplenum"_* | while IFS= read -r line
				do
					echo -e "saving file $line\n"
					filebatch+=( "$line" )
				#done



			# Decompress Forward and Reverse sequence files per sample and condition

			#for file in "${filebatch[@]}";
				#do

				file="$line"

					echo -e "\nDecompressing sample number $samplenum in $condition read files\n"; 
			
					echo "$file"
					tfilename=$(basename "$file")

					# Making fastqc reports of each read file
					
					# fastqc "$file" -o "$path_to_fq"/pre_filter_fastqc

					# Decompress reads per sample and condition
					
					dfile="${tfilename%.fq.gz}".fq

					
					gzip -dkc "$file" >"$path_to_draw_reads"/"$dfile"


				done

				echo -e "\n\nBUGGGGGGG\n"
 
				# Trimming files per sample and condition
				fastp -m --adapter_fasta "$path_to_adapters" \
				-i "$path_to_draw_reads"/"$condition"rep"$samplenum"_1.fq \
				-I "$path_to_draw_reads"/"$condition"rep"$samplenum"_2.fq \
				-o "$path_to_clean_reads"/"$condition"rep"$samplenum"_1.trimmed.fq \
				-O "$path_to_clean_reads"/"$condition"rep"$samplenum"_2.trimmed.fq \
				-n 15 -q 5 -u 50 -l 150 \
				-j "$path_to_fq"/trimmingreport/HDAC_trimrep.json \
				-h "$path_to_fq"/trimmingreport/HDAC_trimrep.html

					
		

		done

	done

# Making fastqc reports of each read file
echo -e "\n-----------------------------------------\n4) Performing final QC of clean files\n-----------------------------------------\n"
		
fastqc "$path_to_clean_reads"/*.fq -o "$path_to_fq"/post_filter_fastqc

multiqc "$path_to_fq"/post_filter_fastqc/*.zip 2>"$path_to_fq"/post_filter_fastqc

multiqc "$path_to_fq"/pre_filter_fastqc/*.zip 2>"$path_to_fq"/pre_filter_fastqc