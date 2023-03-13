#!/bin/bash

# Author : Natalia Garcia Sanchez
# Date : 5/03/2022
# Description: Script to calculate FPKM/TPM
# Usage : NormalizationPipeline.sh countsdir GTFfile.gtf mode
# mode can either be "FPKM" or "TPM"
# counts have to have .counts extension
# --

# USER INPUT CONTROL FUNCTIONS

help(){
echo -e  "Usage : NormalizationPipeline.sh countsdir GTFfile.gtf mode \n mode can either be FPKM or TPM, counts have to have .counts extension"
}

firstchecks(){
	#This function performs and calls first basic checks of arguments called by user

	echo -e "\n\t1) Performing basic checks"
	local numparam=3
	local paramarray=("$@")
	if [ "${#paramarray[@]}" != "$numparam" ]
		then
		echo "ERROR: 3 param were not called in script"
		help
		exit 1
		else
		local mode="${paramarray[2]}"
		if [[ "${mode^^}" =~ "FPKM" ]] || [[ "${mode^^}" =~ "TPM" ]]
			then
				countdir="${paramarray[0]}"
		        	GTFfile="${paramarray[1]}"
				checkingdirs "$countdir"
				checkingfiles "$GTFfile"
			else
				echo "Incorrect norm mode chosen (not FPKM or TPM). Exiting..."
				help
				exit 1
			fi
		fi
}

checkingdirs(){
	# This function checks if counts from first main argument exists and has counts with .counts extension

	directory="$1"
	if [ -d "$directory" ]
	then
		if [[ ! $(ls "$directory"/*.counts | wc -l) -gt 0 ]];
		then
			echo "No count files found in directory. \nExiting..."
			exit 1
		fi

	else
		echo "Count file directory is empty. Try again with a proper dir. \nExiting..."
		exit 1
	fi
}

checkingfiles(){
	# This function checks if GTF file from second main argument exists

	file="$1"
	echo "$file"
		if [ -e "$file" ];
			then
			echo "found GTF file"
			else
			echo "no GTF file was found"
			exit 1
			fi
}



# MAIN PROCEDURES

joincounts(){
	# This function performs an outer join to merge all counts from all genes in countsdir with .counts extesion
	# Input : countsdir
	# Output : mergedcounts

	echo -e '\n\t2) Joining lines with mergeCounts.R script...'
	# Uses an Rscript that writes the HTSeq statistics and returns a file called merged.counts
	Rscript "mergeCounts.R" "$1"
}

calculatelengths(){
	# This function takes GTF file and calculates gene length to account for normalization factors in later stages
	# Input : GTF file `AnnotationGTF`
	# Output : length_genes.csv with G_ID,length fields created in working directory

	echo -e '\n\t3) Calculating gene lengths from GTF file...'
	local AnnotationGTF="$1"
	cat "$AnnotationGTF" | grep -v "^#" | grep -P "\tgene\t" | cut -d"$(echo -e "\t")" -f 4,5,9 | sed 's/gene_id//g' | sed 's/;.*//g' | sed 's/\"//g' | sed 's/\s\s/\t/g' | tr "$(echo -e "\t")" " " | cut -d" " -f 3,1,2 >tmplengthmap.tsv
	echo -e "gene_id,length" >length_genes.csv
	awk ' FS =" " {print $3 " " $2-$1}' "tmplengthmap.tsv" | tr " " "," >>length_genes.csv 
	
	# Remove tmp files
	rm "tmplengthmap.tsv"
}

normcounts(){
	# This function normalizes merged count file according to mode chosen
	# Input : mergedcounts, csvlengthfile, mode(FPKM|TPM)
	# Output : normalized merged counts csv table in working directory

	echo -e '\n\t4) Normalizing merged gene counts...'
	Rscript "NormalizeCounts.R" "$1" "$2" "$3"
}



# Establish run order. Main function outputs merged counts and normalized counts tables

main() {
	echo -e "\nRunning Normalizing count pipeline..."

	firstchecks "$@"
	joincounts "$1"	
	calculatelengths "$2"
	normcounts "$1/merged/merged.counts" "length_genes.csv" "$3"

	# Remove tmp files
	rm "length_genes.csv"
}


main "$1" "$2" "$3"
