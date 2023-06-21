# Annotation retrieval

Series of Programmatic Access and File parsing scripts used throughout the pipeline
| Code  |  Description | Used in  |  
|---|---|---|
| `Bna_ensembl_ID_gene_mappings.R`  | R Script to retrieve programmatically via Ensembl Biomart API any attribute specified in the arguments from _Brassica napus_ ensembl genes `bnapus_egenes` BiomaRt object  |  Output annotation  | 
| `PANTHER_API.ipynb` |  Notebook to access PANTHER REST API _Brassica napus_ Panther protein annotations and merge them with PANTHER gene-protein sequence classifications obtained in the FTP server via curl |  Output annotation  | 
| `get_genesets.Rmd` | Script to retrieve an exhaustive list of biological knowledge categories (KEGG, GO) associated to _Brassica napus_ genes and _Arabidopsis thaliana_ _Brassica napus_ confidence orhtologs of interest belonging to categories of study |  Gene Set Enrichment - Functional Analysis  | 
| `string_mappings.sh` | Programmatic access to STRINGdb to get high confidence DEG PPI interactions from _Brassica napus_ |  PPI  | 
| `uniprot_mappings.sh` | get _Brassica napus_ gene to protein mappings prior to STRINGdb data acquisition|  PPI  | 
| `retrievalRNA.sh` | parser to get ncRNA annotations from Ensembl Plants gene model annotation file  |  MVvsPE exploration of miRNA  | 
