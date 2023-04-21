#!/bin/bash

# Download ensembl assembly into directory of choice and create hisat2 mapping
# Author : Natalia García Sánchez
# Date : 24/02/2023
# ---


# Bash output color parameter for user interaction
GREEN='\033[0;32m'
NC='\033[0m'

echo -e "\n${GREEN} ~~~~~~~~~~~~ download_map_assemble.sh ~~~~~~~~~~~~\n${NC}"
echo -e "A script for HISAT2 read alignment to Brassica napus reference genome and assembly"
echo -e "\n\n ${GREEN} USAGE : ${NC} bash download_map_assemble.sh ${GREEN} <trimmed_fq_read_directory_path>* <numsamples_per_condition> ${NC} \n\n\t (*) : OPTIONAL arg - IF NOT PROVIDED IN COMMAND LINE, ASKED IN USER INPUT"


# Loading necessary modules
echo -e "${GREEN}\n------------------------------------\nLoading Python, HISAT2, SAMTOOLS, AGAT, PERL, GCC, StringTie dependencies\n------------------------------------${NC}"
module load Python/3.8.6-GCCcore-10.2.0
module load GCCcore/10.2.0
module load Perl/5.32.0-GCCcore-10.2.0-minimal
module load gzip/1.10-GCCcore-8.2.0
module load HTSeq/0.11.1-foss-2019a
module load HISAT2/2.1.0-foss-2019a
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load StringTie/2.0.5-GCCcore-8.2.0
module load AGAT/0.9.2-GCC-11.2.0

# Download reference genome ensembl assembly in fasta format

echo -e "${GREEN}\n------------------------------------\n1) Downloading B napus ensembl assembly genome\n------------------------------------${NC}"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/751/015/GCA_000751015.1_AST_PRJEB5043_v1/GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna.gz
gzip -d GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna.gz
mv GCA_000751015.1_AST_PRJEB5043_v1_genomic.fna genome.fa

# Download ensembl annotation in GFF format and convert 

echo -e "${GREEN}\n------------------------------------\n2)Downloading Brassica napus Ensembl Plants GFF3 file annotation and converting into GTF file \n------------------------------------${NC}"
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/brassica_napus/Brassica_napus.AST_PRJEB5043_v1.56.gff3.gz
gzip -d Brassica_napus.AST_PRJEB5043_v1.56.gff3.gz
echo "Converting GFF3 file to GTF3"
agat_convert_sp_gff2gtf.pl --gff Brassica_napus.AST_PRJEB5043_v1.56.gff3 -o genome.gtf 2>tmp.conversion.log


# Extract genomic sequences corresponding to exon/ss sequences

echo -e "${GREEN}3) Extracting splice sites and exon sequences${NC}"
path_to_hisat="$(whereis hisat2 | grep -oP "/.+/bin")/hisat2"

echo -e "${GREEN}\n------------------------------------\nDownloading and using python helper scripts to extract splice sites and exons in reference genome fasta file\n------------------------------------${NC}"

wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_splice_sites.py
chmod 755 hisat2_extract_splice_sites.py
python3 hisat2_extract_splice_sites.py genome.gtf > genome.ss

wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_exons.py
chmod 755 hisat2_extract_exons.py
python3 hisat2_extract_exons.py genome.gtf > genome.exon

echo -e "${GREEN}\n------------------------------------\n4)Align to reference genome with splice sensitive HISAT2 alignment\n------------------------------------${NC}"

# Checking read data 

if [ $# -eq 0 ]; then  
	
	echo -e "\n------------------------------------\n No command line arguments found, asking for user input... \n------------------------------------\n"
	echo -e "\n${GREEN}> Enter clean read data directory path : \n>${NC}"
	read path_to_clean_reads

	echo -e "\n${GREEN}> Enter number of samples : \n >${NC}"
	read fnumbers
fi

# Making output directories

echo -e "${GREEN}\n------------------------------------\nMaking bam_output assembly_output and flagstat_output directories${NC}\n------------------------------------" 
if ! [ -d "bam_output" ] || ! [ -d "flagstat_output" ] || ! [ -d "assembly_output" ];then
mkdir bam_output
mkdir assembly_output
mkdir flagstat_output
fi

path_to_bam="bam_output"
path_to_flagstat="flagstat_output"
path_to_assembly="assembly_output"


if [ -d "$path_to_clean_reads" ]
then

	if [ $(ls "$path_to_clean_reads"/*.fq | wc -l) -gt 0 ]; 
		then 
			echo "Clean read files found in directory. Proceeding to mapping...";
		
		else
			echo "No clean read files found. \nExiting..."
			exit 1
		fi
	else
		echo "Clean reads directory is empty. Try again with a proper dir"
	exit 1
fi

# ----------

# Mapping 

# -- Creating genome index wit BW comp algorithm: This will create a series of indexing files needed

echo -e "${GREEN}\n------------------------------------\n 5) Creating index files accounting for splicesites and exonic sequences in reference genome \n------------------------------------${NC}"

srun -p medium -t 00:40:00 --mem 200GB hisat2-build -p 16 --exon genome.exon --ss genome.ss genome.fa "$path_to_bam"/genome_tran.idx

# -- Mapping to index genome indexing and converting SAM into a sorted bam in a single nested for loop

echo -e "${GREEN}\n------------------------------------\n 6) Read alignment to reference genome and conversion into sorted BAM output \n------------------------------------${NC}"

array_condition=(C T)

for condition in "${array_condition[@]}"
do
	for (( sample_num=1; sample_num<="$fnumbers"; sample_num++ ))
	do
		# Evaluate mappings to the indexed genome
		# Convert SAM mappings into Binary BAM format and not creating SAM file by pipe streaming
		# Evaluating mapping with 
		# Sort BAM for further downstream analyses

		srun -p long -t 1:00:00 --mem 128GB hisat2 -x "$path_to_bam"/genome_tran.idx \
		-1 "$path_to_clean_reads"/"$condition"rep"$sample_num"_1.trimmed.fq \
		-2 "$path_to_clean_reads"/"$condition"rep"$sample_num"_2.trimmed.fq | \
		  samtools sort -O BAM | \
		  tee "$path_to_bam"/"$condition"rep"$sample_num".bam | \
		  samtools index - "$path_to_bam"/"$condition"rep"$sample_num".bam.bai
	done
done

# Make assessments of HISAT2 alignment in each biological replicate

echo -e "${GREEN}\n------------------------------------\n 7) Make read alignment reports \n------------------------------------${NC}"

for file in "$path_to_bam"/*.bam; 
do 
	filename=$(basename "$file")
	samtools flagstat "$file" >"$path_to_flagstat"/${filename%.bam}.flagstat 2>"$path_to_flagstat"/${filename%.bam}.flagstat.log
done


# Assemble mapped reads into GTF format taking a reference

echo -e "${GREEN}\n------------------------------------\n 8) Assemble and annotate mappings into GTF format using the Ensembl Plants GTF reference functional annotation file \n------------------------------------${NC}"

stringtie -B -G genome.gtf -o "$path_to_assembly"/Bna_final_assemble_HDAC.gtf "$path_to_bam"





