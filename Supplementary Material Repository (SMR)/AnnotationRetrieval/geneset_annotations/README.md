# Geneset annotations

Annotations pertaining to crafted list of gene - functional set annotations jointly represented in PANTHERdb and Ensembl last release (Feb 2023v.56) in RDS/csv files

| file | Representing |
|----------|:-------------:|
| `full_bp_genesets_noorth.rds` |  GO Biological Process (BP) | 
| `full_cc_genesets.rds`  |  GO Cellular Component (CC) |  
| `full_mf_genesets.rds`  | GO Molecular Function (MF)|
| `KEGG_gs_pairs.csv`  | KEGG |

In addition, a series of csv representing the Gene Ontology/KEGG set category terms is found in this folder, namely in the following filepaths:

| file | Representing |
|----------|:-------------:|
| `full_bp_terms.csv` |  GO Biological Process (BP) | 
| `full_cc_terms.csv`  |  GO Cellular Component (CC) |  
| `full_mf_terms.csv`  | GO Molecular Function (MF)|
| `KEGG_gs_descriptions.csv`  | KEGG |

The lists use ENA identifiers, as they are more interoperable
