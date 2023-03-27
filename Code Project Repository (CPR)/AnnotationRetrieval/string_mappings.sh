#!/bin/bash
# Author : Natalia García Sánchez
# Description : String aliases file parsing to map ensembl_gene_identifiers to string identifiers

#
help(){
    echo -e "\n-----\nUsage: stringdb_mappings.sh file_with_id_to_map.csv\n---------\nNote that identifiers need to be in the first column of the file"
}

check_id_file(){
    if [ ! -f "$1" ]
    then 
        echo -e "File not found."
        help()
        exit 1
    fi
}

main(){
    check_id_file "$1"
    # Download Brassica napus STRING aliases file
    wget https://stringdb-static.org/download/protein.aliases.v11.5/3708.protein.aliases.v11.5.txt.gz

    gzip -d 3708.protein.aliases.v11.5.txt.gz


    # FILTERING by alias source
    echo "Possible aliases sources found"
    cat 3708.protein.aliases.v11.5.txt | cut -d"$(echo -e "\t")" -f3 | sort | uniq -c

    # Enter alias to map with file
    echo -e "Enter one of these possible alias to map with file"
    
    # In our case our identifiers will work with BLAST_Uniprot_GN_ORFNames alias, so we will adjust this parameter. This line can be uncommented if another filtering is needed
    filter_source="BLAST_UniProt_GN_ORFNames"

    # Optionally, one can uncomment these lines of code if another filter source is needeed
    # for user input:
    # filter_source=read()
    # for custom change of filter
    # filter_source=

    cat 3708.protein.aliases.v11.5.txt | grep -a "$filter_source" | awk 'FS= " "{print $2","$1}' | sort >tmpmappingfile.csv
    
    sort "$1" >tmpsortedid.csv
    
    # header
    echo "gene_id,stringdb_id" >stringdb_mappings.csv
    # join to id file
    join -t, tmpsortedid.csv tmpmappingfile.csv >stringdb_mappings.csv

    # rm temporary files
    rm tmpsortedid.csv
    rm tmpmappingfile.csv

    }

main $1