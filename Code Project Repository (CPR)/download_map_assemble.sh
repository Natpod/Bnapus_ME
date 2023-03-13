#!/bin/bash

# Download ensembl assembly into directory of choice and create hisat2 mapping
# Author : Natalia García Sánchez
# Date : 24/02/2023
# --- 

# Download reference genome ensembl assembly in fasta format
echo "Downloading B napus ensembl assembly"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/751/015/GCA_000751015.1_AST_PRJEB5043_v1/GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna.gz
gzip -d GCF_020379485.1_Da-Ae_genomic.fna.gz
mv GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna.gz genome.fa

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/379/485/GCF_020379485.1_Da-Ae/GCF_020379485.1_Da-Ae_genomic.fna.gz
#gzip -d GCF_020379485.1_Da-Ae_genomic.fna.gz
#mv GCF_020379485.1_Da-Ae_genomic.fna genome.fa

# Checking read data 
echo "Enter clean read data directory path : \n>"
path_to_clean_reads=read()

if [ -d "$path_to_clean_reads" ]
then

		if [ $(ls *.fq | wc -l) -gt 0 ]; 
		then 
			echo "Clean read files found in directory. Proceeding to mapping...";
		
		else
			echo "No clean read files found. \nExiting..."
			exit 1
	else
		echo "Clean reads directory is empty. Try again with a proper dir"
	exit 1
fi

# ----------

# Mapping 

# -- Creating genome index wit BW comp algorithm: This will create a series of indexing files needed

bwa index genome.fa

# -- Mapping to index genome and converting formats

array_condition=('C' 'T')

for condition in "${array_condition[@]}"
do
	for sample_num in {1..3}
	do
		# Evaluate mappings to the indexed genome

		bwa mem $path_to_clean_reads/$condition"rep_"$sample_num"_1.clean".fq $path_to_clean_reads/$condition"rep_"$sample_num"_2.clean".fq > $condition"rep"$sample_num.sam
		2>$condition"rep"$sample_num.err

		# Convert SAM mappings into Binary BAM format 
		
		samtools view -b -h $condition"rep"$sample_num.sam > $condition"rep"$sample_num.bam
		
		# Remove heavy sam file

		rm $condition"rep"$sample_num.bam

		# Sort BAM for further downstream analyses

		samtools sort $condition"rep"$sample_num.bam > $condition"rep"$sample_num.sam
	
	done
done

# -- Evaluating mapping

# Creating tmp flagstat reports

for file in *.bam; do samtools flagstat $file > $file.flagstat_report.txt ;done

# take first column from flagstat

cat $(ls *flagstat_report* | head -n1) | cut -f1 >flagstat_report

# Join all flagstat reports into a single document and remove the single files

for file in *.flagstat_report.txt; do join flagstat_report $file; done
rm *.flagstat_report.txt






