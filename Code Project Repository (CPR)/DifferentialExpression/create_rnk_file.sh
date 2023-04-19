#!/bin/bash

# Description : rnkgen.sh converts a differential gene expression csv file (CSV) into a rank file (RNK)
# Author : Natalia García Sánchez
# Date : 12/04/2023
# ----


help(){
	echo -e "${GREEN}\n~~~~~~~~ create_rnk_file.sh ~~~~~~~~~~~\n\tUsage:${NC} create_rnk_file.sh${GREEN} <DeSeq2_wth_external_gene_name_and_stat_output.csv> <external_gene_name_colnum_in_csv> <stat_colnum_in_csv>${NC}"
	exit 1
}

check_DE_file_exists(){
if [ -f "$1" ]; 
	then 
		echo "Processing DeSeq2 Output file..."
	else 
		echo -e "No DEG Deseq2 csv file found"
		help
	fi
}

main(){

check_DE_file_exists "$1"
cat "$1" | awk -v I="$ID" -v P="$P" 'FS = "," {print $I","$P}' | awk '$1!="NA" && $2!="NA"' | sed 's/\"//g' | sort -t, -k2gr   >bna_HDAC_deseq2stat.rnk

}

# Bash user input console colors 
	
GREEN='\033[0;32m'
NC='\033[0m'

#Specify the external gene ID column
ID=$2
#Specify the DeSeq2 wald stat "stat" column num
P=$3

main "$1" "$2" "$3"

