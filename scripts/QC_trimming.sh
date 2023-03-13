#!/bin/bash

# 1st stage of RNAseq pipeline script : working with raw reads, trimming and doing Quality Control (QC)
# Author : Natalia García Sánchez
# Date : 23/02/2023
# ---

# Line necessary to run 
conda source activate

#---

# ENTER INPUT AND OUTPUT DIRECTORIES

#echo 'Enter the directory where your raw read files can be found'
#read path_to_raw_reads
path_to_raw_reads='/home/famgarcia/TFM/data/Testillano col/Report Novogene RNAseq Deacetilation Inhibition treatment vs Control/rawdata'

#echo 'Enter the directory where your want to decompress your raw read files'
#read path_to_raw_reads
path_decompress_raw='/home/famgarcia/TFM/owndata/raw'

#echo 'Enter the directory where you want your fq read files can be found'
#read path_to_raw_reads
path_clean_reads='/home/famgarcia/TFM/owndata/clean'

ask_user_input_yn(){
echo -e "Directory is not empty\n\tCurrent directory will be overwritten. Want to continue? (y/n) "; 
				read response1

				case $edad in
				  Y|y|[Yy]es)
					continue
				  ;;
				  N|n|[Nn]o)
				    echo "¡Correcto!"
				  	exit(0)
				  ;;
				  *)
				    ask_user_input_yn
				  ;;
				esac
				
			}

#if directory exists
if [ -d "$path_decompress_raw" ] 
	then 
		echo "Decompressed raw read directory exists"
		if [ -z "$(ls -A -- "$path_decompress_raw")" ]; 
			then 
				ask_user_input_yn
			else 
				echo "Making directories"; 
				mkdir "$path_decompress_raw"
			fi

# Making clean reads output directory
mkdir "$path_clean_reads"
mkdir '~/TFM/owndata/pre_filter_fastqc'

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

			for file in $condition"rep_"$samplenum*
			do
				# Move to raw reads directory and decompress reads per sample and condition
				
				cd $path_to_raw_reads
				gzip -dkc "$file" >"/home/famgarcia/TFM/owndata/raw/${file%.fq.gz}.fq"

				# Making fastqc reports of each read file
				
				conda activate TFM
				cd "$path_decompress_raw"/
				fastqc $file -o '~/TFM/owndata/pre_filter_fastqc'

			done
		for condition in "${array_conditions[@]}";
		do 

			tmp_read_files_FR=($condition"rep_"$sample_num*)

			# Trimming files
			trimmomatic PE -threads 4 -trimlog logtrim.log ${tmp_read_files_FR[0]} ${tmp_read_files_FR[1]}			\
			${tmp_read_files_FR[0]%.fastq}.trimmed.fastq ${tmp_read_files_FR[0]%.fastq}_unpair.trimmed.fastq \
			${tmp_read_files_FR[1]%.fastq}.trimmed.fastq ${tmp_read_files_FR[0]%.fastq}_unpair.trimmed.fastq \
			ILLUMINACLIP:"/home/famgarcia/TFM/data/Testillano col/Report Novogene RNAseq Deacetilation Inhibition treatment vs Control/adapters/adapters_HDAC.fa" \

			# por ver- adaptadores + output directory
		done
		done

	done

# Building an array that establishes the conditions of the experiment ('Control':'C', 'Treatment-DHDAC' : 'T')
array_conditions=('C' 'T')
for sample_num in {1..3}
do
	for condition in "${array_conditions[@]}"
	do
		tmp_read_files_FR=$($condition"rep_"$sample_num*)
		fastp -m --include_unmerged -a --adapter_fasta "/home/famgarcia/TFM/data/Testillano col/Report Novogene RNAseq Deacetilation Inhibition treatment vs Control/adapters/adapters_HDAC.fa" -i ${tmp_read_files_FR[0]} -I ${tmp_read_files_FR[1]} --merged_out ${tmp_read_files_FR[0]}.mergedout.fastqc -n 15 -q 5 -u 50
		done
done




fastp -m --include_unmerged --adapter_fasta "/home/famgarcia/TFM/data/Testillano col/Report Novogene RNAseq Deacetilation Inhibition treatment vs Control/adapters/adapters_HDAC.fa" -i Crep1_1.fq -I Crep1_2.fq --merged_out Crep1.mergedout.fastqc -n 15 -q 5 -u 50
		


#gzip -dkc "$file" >"/home/famgarcia/TFM/owndata/raw/${file%.fq.gz}.clean.fq.gz"
 mkdir fastqc

for file in *.fq*
do 
	fastqc $file -o fastqc

done
for file in fastqc/*.zip; do unzip $file; done

for file in Crep1*; do gzip -dkc "$file" >"1rep_adaptortrimming/"${file%.fq.gz}.clean.fq.gz; done
