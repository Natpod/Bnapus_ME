This folder contains a set of files accounting for gene annotations used throughout the analysis

| Tables   | Annotations Used |  Source | Acquisition through | License|
|----------|:-------------:|------:|------|------|
| `ensembl_centric_annotations.tsv` | Gene model annotation descriptions, Uniprot code mappings, Entrez Gene ID mappings  |  Ensembl release v.56 | programmatic access to BiomaRt R API in `get_genesets.Rmd`, bnapus_eg_gene mart | |
| `ncbi_centric_annotations.tsv` | Entrez Gene ID mappings, gene aliases and Description | dataprovider: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ (2022/10/26) NCBI/UniProt source gene ID based annotations about Brassica napus  | programmatic access through AnnotationHub _Brassica napus_ OrgDB AH107389 object with AnnotationDbi and AnnotationHub | |
| `PANTHER_proteins.xlsx` | PANTHER Protein class, Family class, Superfamily class | PANTHER v17 | programmatic access to PANTHER v17 RESTdb API through python requests library + gene protein sequence classifications in PANTHERdb v17 FTP server | |
| `high_confidence_homologues` | _Brassica napus - Arabidopsis thaliana_ reciprocal high confidence orthologs and %id between gene sequence pairs categorized in Ensembl Compara |  (Ensembl Plants release v.56) | programmatic access to BiomaRt R API in `get_genesets.Rmd`, bnapus_eg_gene mart | |

    