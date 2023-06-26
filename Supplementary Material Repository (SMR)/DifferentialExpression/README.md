# Differential Expression
This folder contains results and visualizations from the **Differential Expression Analysis** performed by DESeq2 scripts, described in the following table:
| Folder   |      Dataset      |  Description |
|----------|:-------------:|------:|
| HeatStress | DEA of Ctrl Proembryos vs Vacuolated Microspore expression obtained in [DESEQ2_PE_HDAC.Rmd](https://github.com/Natpod/Bnapus_ME/blob/main/Code%20Project%20Repository%20(CPR)/DifferentialExpression/DESEQ_PE_HDAC.Rmd) | Annotated DEG lists, Clustering visualization for outlier detection, dispersion plots, volcano plots, MA plots  |
| HeatStress_SAHA |  DEA of SAHA treated Proembryos vs Vacuolated Microspore expression obtained in [DESEQ2_PE_HDAC.Rmd](https://github.com/Natpod/Bnapus_ME/blob/main/Code%20Project%20Repository%20(CPR)/DifferentialExpression/DESEQ_PE_HDAC.Rmd)  | Annotated DEG lists, Clustering visualization for outlier detection, dispersion plots, volcano plots, MA plots |
| PE_SAHA | DEA of SAHA treated Proembryos vs Vacuolated Microspore expression obtained in [DESEQ2_DE_analysis.Rmd](https://github.com/Natpod/Bnapus_ME/blob/main/Code%20Project%20Repository%20(CPR)/DifferentialExpression/DESEQ2_DE_analysis.Rmd) | Annotated DEG lists, Clustering visualization for outlier detection, dispersion plots, volcano plots, MA plots |


In addition, the `EDA_DEG_Subsets` folder contains results and Visualizations from the following Exploratory Data Analysis (EDA) inspection of intersecting DEGs in the three DEG lists, performed in [DEG_exploratory_analysis.R](https://github.com/Natpod/Bnapus_ME/blob/main/Code%20Project%20Repository%20(CPR)/DifferentialExpression/DEG_exploratory_analysis.R):

- DEGs in HeatStress
- DEGs in HeatStress + SAHA treatment
- DEGs in Proembryos + SAHA treatment

