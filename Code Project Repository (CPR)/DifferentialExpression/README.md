Application of Differential Expression Analysis (DEA) on HDACi sequemce data (SAHA treated and untreated Proembryos) and all count data, adjusted for batch effects.

Prupose : Identifying Differentially Expressed Genes (DEGs) from:
1) SAHA treated proembryo replicates (Trep1, Trep2, Trep3) vs Ctrl proembryos (Crep1, Crep2, Crep3) - *accounts for SAHA effect on proembryos*
2) All adjusted Control Proembryos in the two sequencing projects (Crep1, Crep2, Crep3, PErep1, PErep2) vs Vacuolated Microspore replicates (VMrep1, VMrep2, VMrep3)
3) Adjusted SAHA treated proembryo replicates(Trep1, Trep2, Trep3) vs Vacuolated Microspore replicates (VMrep1, VMrep2, VMrep3)

| Code  | Pipeline stage  | Description  |  
|---|---|---|
| `DESEQ2_DE_analysis.Rmd` |  Differential expression  | DEA analysis for points 2 and 3 |
| `DESEQ_DE_analysis.Rmd` | Differential expression  | DEA analysis for 1  |
| `DEG_exploratory_analysis.R`| Analysis of DEA results | makes Venn Diagrams for all 1/2/3 DEG point combinations, and between 1/2 segmented by type of expression|
