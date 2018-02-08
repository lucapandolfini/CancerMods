# CancerMods
Multiple stats on TCGA data via TCGAbiolinks package

### TCGAbiolinks.R
Access the TCGA dataset through TCGAbiolinks package.

* Get mRNA query
* Get miRNA query
* Run queries
* Parse data in more friendly dataframes:
   * TCGA_curated_dataset.RData (Samples with matched Clinical data, mRNA and miRNA expression)
   * TCGA_curated_dataset.mRNAonly.RData (Samples with matched Clinical data and mRNA expression)

### TCGAcurated.R
Functions for data extraction from the curated datasets (_TCGA_curated_dataset.*.RData_)

`extract_curated(project, genes, mirnas, clinical_data)`\
Extract data from specific project for genes, miRNA and clinical data

`extract_all_curated(genes, mirnas, clinical_data)`\
Extract arbitrary data from all the projects available in a single dataframe 

`plotExprCancer(gene, max_y_value=NULL)`\
Plot gene expression profile across all tumor projects

`plotExprCancerVsNormal(gene, max_y_value=NULL)`\
Plot gene expression profile across all projects for paired tumor vs normal samples, when available

`plotExprCancerVsMetastasis(gene, max_y_value=NULL)`\
Plot gene expression profile across all projects for paired primary tumor vs metastasis samples, when available

`doKaplanMeier(cfu, label, stoptime_analysis=NULL, stoptime_plot=NULL)`\
Parse survival data, fit Cox Proportional hazards model and plot results as Kaplan-Meier plot

`geneSurvival(project, gene, percentile, stoptime_analysis=NULL, stoptime_plot=NULL)`\
Wrapper function for _doKaplanMeier()_ to extract percentiles of gene expression and plotting survival resuls

`microSurvival(project, mirna, percentile, stoptime_analysis=NULL, stoptime_plot=NULL)`\
Wrapper function for _doKaplanMeier()_ to extract percentiles of miRNA expression and plotting survival resuls

`plotPairedBiKM(project, gene, mirna, percentile, stoptime_analysis=NULL, stoptime_plot=NULL)`\
Plot survival analysis for a matched gene-miRNA pair in a specific project

`generateAllPairedKMs<-function(gene, micro, percentile, stoptime_analysis=NULL, stoptime_plot=NULL)`\
Plot survival analysis for a matched gene-miRNA pair in all the projects available

## Files available for download upon request:

md5sum                            |  File
--------------------------------- | ----------------------------------------
51d65a22b4f8c5772d97338ec4df8091  |	 enzymes.txt
d647a96306d6c7ba505b612cc64a41e4	|	 TCGA_curated_dataset.RData
fd874ed9e5ee959262f595d31bbd62a3	|  TCGA_curated_dataset.mRNAonly.RData
