# Enrichment Analysis

### Enrichment of DEG results with 
* Overrepresentation Analysis (ORA) using ShinyGO Web Tool to apply statistical inference for representation of gene sets from multiple curated vocabularies (GO - BP/CC/MF, KEGG Pathways, CycGSAD, LitGSAD, TF-GSAD). FDR q-val threshold : 0.05 $\right_arrow()$ _Does not take expression data into account, only DEG gene list_

* Gene Set Analysis (GSAD): Takes expression data into account with a total gene list with minimum filtering. Ranks genes and applies a statistical method to score the gene sets for a more exhaustive gene set control.