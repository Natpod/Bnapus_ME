#!/usr/bin/env python3


# Author: Natalia García Sánchez
# Date : 20/04/2023
# Description : Use as helper script to merge File with EnsemblID -- Other_externalID mappings from biomart API service download


print('Usage : get_DE_genelists.py <DeSeq_totalresult_output.csv> <ID_mappings_list>')

import sys
import os
import pandas as pd

print('Reading files...')
filepath_res = str(sys.argv[1])
filepath_map = str(sys.argv[2])
dirpath = str(os.path.dirname(filepath_res))

DE_results = pd.read_csv(filepath_res)
mappings = pd.read_csv(filepath_map)
mappings.rename(columns={mappings.columns[0] : 'ensembl_gene_id', mappings.columns[1] : 'external_gene_name'}, inplace=True)

res_annot = pd.merge(DE_results, mappings, how='inner', on='ensembl_gene_id')

res_annot.to_csv(dirpath+'/Final_condition_PErep_vs_VMrep_total_woPErep3_anot.csv', index=False)


