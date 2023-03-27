#!/bin/bash
# The recipes below use wget to download files 
# Please change these env variables to use other tools ie curl
EXE="wget"
ARGSDEF=" -c " 
ARGSTOFILE=" -O "
ARGSTDOUT=" --quiet $ARGSTOFILE - "

# set servers & division
SERVER="ftp://ftp.ensemblgenomes.org/pub"
DIV=plants
BIOMARTSERVICE="http://plants.ensembl.org/biomart/martservice"

# get Ensembl Plants current release number
SUMFILE="${SERVER}/${DIV}/current/summary.txt"
RELEASE=$($EXE $ARGSTDOUT $SUMFILE | \
	perl -lne 'if(/Release (\d+) of Ensembl/){ print $1 }')

echo "Ensembl Plants current release number used in script: $RELEASE"

# work out Ensembl Genomes release
EGRELEASE=$((RELEASE - 53));

# optional arguments, if any
#OPTARG=$1

echo "EGRELEASE=${EGRELEASE} OPTARG=${OPTARG}"
echo

# set example species
SPECIES=Brassica_napus

UNIPTSV="${SPECIES}*.uniprot.tsv.gz"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/tsv/${SPECIES,,}/$UNIPTSV"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

gzip -d "${SPECIES}*.uniprot.tsv.gz"


cat "${SPECIES}*.uniprot.tsv.gz" | cut -d"$(echo -e "\t")" -f5 | sort | uniq -c
