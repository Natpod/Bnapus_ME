#!/bin/bash

# Download ensembl assembly into directory of choice and create hisat2 mapping
# Author : Natalia García Sánchez
# Date : 24/02/2023
# --- 


# Download reference genome ensembl assembly in fasta format

# Uncomment if you'd rather use the latest ncbi b napus accession assembly
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/379/485/GCF_020379485.1_Da-Ae/GCF_020379485.1_Da-Ae_genomic.fna.gz
#gzip -d GCF_020379485.1_Da-Ae_genomic.fna.gz
#mv GCF_020379485.1_Da-Ae_genomic.fna genome.fa

echo "1) Downloading B napus ensembl assembly"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/751/015/GCA_000751015.1_AST_PRJEB5043_v1/GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna.gz
gzip -d GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna.gz
mv GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna.gz genome.fa

# Download ensembl annotation in GFF format and convert 
echo "Ensembl annotation GTF file will be downloaded"
echo "Downloading B napus ensembl plants annnotation GFF3 file"
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/brassica_napus/Brassica_napus.AST_PRJEB5043_v1.56.gff3.gz
gzip -d Brassica_napus.AST_PRJEB5043_v1.56.gff3.gz
echo "Converting GFF3 file to GTF3"
agat_convert_sp_gff2gtf.pl --gff Brassica_napus.AST_PRJEB5043_v1.56.gff3.gz -o genome.gtf

# Extract genomic sequences corresponding to exon/ss sequences
echo "2) Extracting splice sites and exon sequences"
echo "write directory path to hisat2 scripts hisat2_extract_splice_sites and hisat2_extract_exons"
path_to_hisat=read()
"$path_to_hisat"/hisat2_extract_splice_sites.py genome.gtf > genome.ss
"$path_to_hisat"/hisat2_extract_exons.py genome.gtf > genome.exon


# Checking read data 
echo "> Enter clean read data directory path : \n>"
path_to_clean_reads=read()

# Checking read data 
echo "> Enter data directory path where you want to store your indexing and bam files: \n>"
path_to_bam=read()

# Checking read data 
echo "> Enter data directory path where you want to store bam files flagstat_reports: \n>"
path_to_flagstat=read()

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

hisat2-build -p 16 --exon genome.exon --ss genome.ss genome.fa "$path_to_bam"/genome_tran.idx

# -- Mapping to index genome and converting formats

array_condition=(C T)

for condition in "${array_condition[@]}"
do
	for sample_num in {1..3}
	do
		# Evaluate mappings to the indexed genome
		# Convert SAM mappings into Binary BAM format and not creating SAM file by pipe streaming
		# Evaluating mapping with 
		# Sort BAM for further downstream analyses

		hisat2 -x genome_tran.idx \
		-1 "$path_to_clean_reads"/"$condition"rep_"$sample_num"_1.trimmed.fq \
		-2 "$path_to_clean_reads"/"$condition"rep_"$sample_num"_2.trimmed.fq | \
		  tee >(samtools flagstat - > "$path_to_flagstat"/"$condition"rep_"$sample_num".flagstat) | \
		  samtools sort -O BAM | \
		  tee "$path_to_bam"/"$condition"rep_"$sample_num".bam | \
		  samtools index - "$path_to_bam"/"$condition"rep_"$sample_num".bam.bai
	done
done





