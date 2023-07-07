# Gene Set Enrichment Analysis

The purpose of the following scripts is to perform a functional analysis of the expression data through curated categories taken from existing knowledge of biological processes.
This is done in order to uncover mechanisms underlying in our expression data, more specifically by answering the following questions with the corresponding tools that are highlighted:

- In my list of Differentially Expressed Genes (DEG), which category is overrepresented ? **(ORA)**
- In my full list of genes ranked by expression,  which category is overrepresented ? **(GSEA)**
- How present is expression of genes in a term/pathway for all my list of heterogeneous condition data ? **(GSVA)**


This section contains the following scripts
| Code  | Pipeline Stage  | Description  |
|---|---|---|
| `GSEA`  | GeneEnricment  | Performs GSEA in the three lists, for the selected hormone and epigenetic terms  |   
| `gsva.Rmd`  | Functional Enrichment Analysis  |  Performs Gene Set Variation Analysis (GSVA) on batch adjusted matrix expression data from both experiments | 
| `GSVA_viz.ipynb`  | Functional Enrichment Analysis  | Visualizes GSVA resuts | 
| `ORA.Rmd`  | Functional Enrichment Analysis  | Performs ORA over list of HDACi/Stress DEGs, visualizes the result |  

Prior to running these scripts, a collection of genesets must be performed. The following table contains the code needed for this, as well as for the selection of full gene ranked expression data for the GSEA analysis

| Code  | Pipeline Stage  | Description  |
|---|---|---|
| `get_genesets.Rmd` | Data collection pre GSE | Obtaining KEGG and Gene Ontology (GO) associations between genes and Celullar Component (CC) / Molecular Function categories. Getting _A- thaliana - B. napus _ homologs of interest in specific categories and integrating them in the merged gene set categories obtainesd from GO Biological Process (BP) |
| `create_rnk_file.sh` | Data preprocessing pre GSEA | Utility bash scripts to rank and sort genes by expression data scoring metric (Wald Statistic) |
| `get_DE_genelists.sh` | Data preprocessing pre ORA | Utility bash script to get Differentially Expressed Genes for (ORA) Gene Set Enrichment (GSE) |

More specifically, it will deal with functional terms that are related with :
- hormonal response (including but not limited to auxin)
- epigenetic mechanisms (chromatin modifications, nucleosome, histone modifications, DNA methylation mediated or not mediated by RdRp )
- developmental terms
- Other terms with relevance to ME : lipid storage, cell wall modifications, oxidative stress, autophagy

A list of terms used for filtering Enrichment results can be  
