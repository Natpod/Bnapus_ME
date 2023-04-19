#!/bin/bash

# Description : rnkgen.sh converts a differential gene expression csv file (CSV) into a rank file (RNK)
# Author : Natalia García Sánchez
# Date : 12/04/2023
# ----


help(){
	echo -e "${GREEN}\n~~~~~~~~ create_rnk_file.sh ~~~~~~~~~~~\n\tUsage:${NC} create_rnk_file.sh${GREEN} <DeSeq2_wth_external_gene_name_and_stat_output.csv> ${NC}"
	exit 1
}

check_DE_file_exists(){
if [ -z "$1" ]; 
	then 
		echo "Processing DeSeq2 Output file..."
	else 
		echo -e "No DEG Deseq2 csv file found"
		help
	fi
}

main(){

check_DE_file_exists "$1"
cat "$1"| sed 's/\"//g' | awk -v I="$ID" -v P="$P" 'FS = "," {print $I"\t"$P}' \
| awk '$1!="NA" && $2!="NA"' | sort -k2gr   >bna_HDAC_deseq2stat.rnk

}

# Bash user input console colors 
	
GREEN='\033[0;32m'
NC='\033[0m'

#Specify the external gene ID column
ID=$2
#Specify the DeSeq2 wald stat "stat" column num
P=$3
main "$1" "$2" "$3"

