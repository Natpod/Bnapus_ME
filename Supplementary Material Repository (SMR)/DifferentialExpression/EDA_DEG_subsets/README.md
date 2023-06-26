This folder contains a series of Venn Diagram representations of Differentially Expressed Genes (DEG) intersections among the three lists described in the parent folder `DifferentialExpression`.

- `Intersection_DEG_Embryogenesis_notsegmented.png`: VM vs SAHA Treated PE and VM vs Ctrl PE DEG lists intersections
- `Intersection_DEG_Embryogenesis_segmented.png`:VM vs SAHA Treated PE and VM vs Ctrl PE DEG lists intersections, separated as well by type of expression (increased or decreased)
- `Intersection_DEG_all3subsets.png`: all three DEG lists

![alt text](https://github.com/Natpod/Bnapus_ME/blob/main/Supplementary%20Material%20Repository%20(SMR)/DifferentialExpression/EDA_DEG_subsets/Intersection_DEG_Embryogenesis_segmented.png)

It also contains a series of tsv with annotations of each intersecting subset of DEGs of interest


| Tables   |      Description      | Approach to biological question |
|----------|:-------------:|------|
| `intersec_all3subsets.tsv` | DEGs intersecting across the three lists (Proembryo SAHA treatment, ME with and without SAHA treatment)  | conserved DEGs |
| `intersec_HeatStress_Embryogenesis_outersec_PE_treatment.tsv` | DEGs intersecting across ME DEG lists, but not PE SAHA treatment DEG list  | conserved DEGs only in Embryogenesis |
| `intersec_SAHA_treatment_outersec_HeatStress.tsv` | DEGs intersecting across SAHA treatment DEG lists, but not conventional Embryogenesis "HeatStress" DEG list | conserved DEGs in SAHA treatment |
| `intersec_Embryogenesis_all.tsv` | DEGs intersecting across Embryogenesis lists "HeatStress" and "HeatStress-SAHA" | conserved DEGs in Embryogenesis |
| `outersec_Embryogenesis_in_HeatStress.tsv` | DEGs in Conventional Embryogenesis list "HeatStress" that are not observed with SAHA treatment in "HeatStress-SAHA" DEG list | SAHA effect in Embryogenesis |
| `outersec_Embryogenesis_in_HeatStress_SAHA.tsv` |  DEGs in Embryogenesis with SAHA treatment that are not observed with Conventional HeatStress-only treatment in "HeatStress" DEG list | SAHA effect in Embryogenesis |
| `intersec_Embryogenesis_DOWN_HeatStress_SAHA_UP_Heatstress_all.tsv` | DEGs that are downregulated in Embryogenesis with SAHA treatment but upregulated in Conventional Embryogenesis | SAHA effect in Embryogenesis |
    
