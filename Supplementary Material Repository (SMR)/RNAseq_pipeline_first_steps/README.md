Reports related to the fist stages of the analysis (Quality Check + Trimming, Mapping to reference genome, Read Count) available in the following files:
| File   | Description  |
|---|---|
| `trimming_report_summary.csv`  | A collection of data statistics from the trimming process with fastp - Q20,Q30 values are reported.  |
| `multiqc_report_before_trimming` |  A collection of data statistics from all fastqc files in all replicates and Forward and Reverse read files merged in a single file - after trimming |
| `multiqc_report_before_trimming`  | A collection of data statistics from all fastqc files in all replicates and Forward and Reverse read files merged in a single file - before trimming  |
| `flagstat_reports_summary.ods`  |  Statistics of flagstat returned by samtools in HiSAT2 bam files.  |
| `htseq-count_log_all.csv` | Read count statistic - showing count of ambiguous, low quality read mappings or not found read mappings  |


`adapters_HDAC.fa`  indicates the adapters used for the analysis in the HDACi dataset
Folders with more specific quality report visualizations of per base content,... etc.
