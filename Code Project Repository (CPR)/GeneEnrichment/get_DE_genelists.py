#!/usr/bin/python

# Author: Natalia García Sánchez
# Date : 20/04/2023
# Description : Get genes from DeSeq2 result, retrieve three gene lists 
# ---
# Background genes : all genes detected after a minimum expression filtering by DeSeq2 (e.g. log10(count)>1)
# Disregulated genes : : all following padj and LFC thresholds
# 2 Gene lists with specific differential expression (up and downregulated lists) 

print('Usage : get_DE_genelists.py <DeSeq_totalresult_output.csv>')

import sys
import os
import pandas as pd

print('Reading files...')
filepath = str(sys.argv[1])

dirpath = str(os.path.dirname(filepath))
DE_results = pd.read_csv(filepath)

print('Filtering and generating gene lists in .txt in input folder...')

# Background genes 
bfname = dirpath+'/background_genes_Stress.txt'
DE_results['ensembl_gene_id'].to_csv(bfname, header=False, index=False)

# Disregulated genes
dfname = dirpath+'/genes_dysregulated_Stress.txt'
DE_results = DE_results[abs(DE_results['log2FoldChange']) >= 1]
DE_results = DE_results[DE_results['padj'] < 0.05]
DE_results['ensembl_gene_id'].to_csv(dfname, header=False, index=False)

# Up and downregulated genes
up = DE_results[DE_results['log2FoldChange'] > 0]
down = DE_results[DE_results['log2FoldChange'] < 0]

ufname = dirpath+'/genes_upregulated_Stress.txt'
wfname = dirpath+'/genes_downregulated_Stress.txt'

up['ensembl_gene_id'].to_csv(ufname, header=False, index=False)
down['ensembl_gene_id'].to_csv(wfname, header=False, index=False)
