#!/bin/bash

# Description: 1st stage of RNAseq pipeline script : 
#	- working with raw reads, 
#	- trimming 
# 	- Performing Quality Control (QC) pre and post filtering

# Author : Natalia García Sánchez
# Date : 23/02/2023
# ---


#---

ask_user_input_yn(){
echo -e "Directory is not empty\n\tCurrent directory will be overwritten. Want to continue? (y/n) "; 
				read response1

				case $edad in
				  Y|y|[Yy]es)
					continue
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
#if input directory exists
if [ -d "$1" ] 
	then 
		echo "$1 directory exists"
	else 
		echo "$1 directory not found. Exiting..."
		exit 1
	fi
		}

check_Odir_exists(){
if [ -d "$1" ] 
	then 
		echo "$1 directory exists"
	
	if [ -z "$(ls -A -- "$1")" ]; 
				then 
					ask_user_input_yn
				else 
					echo "Making directories"; 
					rm -r "$1"
					mkdir "$1"
				fi
			fi
}

# MAIN

# ENTER INPUT AND OUTPUT DIRECTORIES

echo 'Enter the directory where your raw read files can be found'
read path_to_raw_reads
check_Idir_exists "$path_to_raw_reads"

echo 'Enter the filepath where your sequencing adapters in fasta format can be found'
read path_to_adapters
check_Idir_exists "$path_to_adapters"


echo 'Enter the directory where your want to decompress your raw read files'
read path_to_draw_reads
check_Odir_exists "$path_to_draw_reads"

echo 'Enter the directory where you want your fq read files can be found'
read path_to_fq
check_Odir_exists "$path_to_fq"

echo 'Enter the directory where your want to put your clean read files'
read path_to_clean_reads
check_Odir_exists "$path_to_clean_reads"

# Making fq reads output directory
mkdir "$path_to_fq"/pre_filter_fastqc
mkdir "$path_to_fq"/post_filter_fastqc

# Building an array that establishes the conditions of the experiment ('Control':'C', 'Treatment-DHDAC' : 'T')
array_conditions=('C' 'T')


# For every F/R pair of read .fq files in the paired sequencing, sample and condition: 
# 1) Files are decompressed
# 2) FQ reports are made for each file


for sample_num in {1..3};
	do

		for condition in "${array_conditions[@]}";
		do 
			echo "Decompressing sample number $sample_num in $condition read files"; 

			for file in "$path_to_raw_reads"/"$condition"rep_"$samplenum*"
			do
				# Making fastqc reports of each read file
				
				fastqc "$path_to_raw_reads"/"$file" -o "$path_to_fq"/pre_filter_fastqc

				# Decompress reads per sample and condition
				
				gzip -dkc "$path_to_raw_reads"/"$file" >"$path_to_draw_reads"/"${file%.fq.gz}".fq



			done
		done

		for condition in "${array_conditions[@]}";
		do 

			# Trimming files
			fastp -m --include_unmerged --adapter_fasta "$path_to_adapters" \
			-i "$path_to_draw_reads"/"$condition"rep"$sample_num"_1.fq \
			-I "$path_to_draw_reads"/"$condition"rep"$sample_num"_2.fq \
			-o "$path_to_clean_reads"/"$condition"rep"$sample_num"_1.trimmed.fq \
			-O "$path_to_clean_reads"/"$condition"rep"$sample_num"_2.trimmed.fq \
			-n 15 -q 5 -u 50

			
		done
		

	done

# Making fastqc reports of each read file
				
fastqc "$path_to_clean_reads"/*.fq -o "$path_to_fq"/post_filter_fastqc

