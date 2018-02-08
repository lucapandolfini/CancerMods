# CancerMods
Multiple stats on TCGA data via TCGAbiolinks package

`TCGAbiolinks.R`\
Access the TCGA dataset through TCGAbiolinks package.
* Get mRNA query
* Get miRNA query
* Run queries
* Parse data in more friendly dataframes:
   * TCGA_curated_dataset.RData (Samples with matched Clinical data, mRNA and miRNA expression)
   * TCGA_curated_dataset.mRNAonly.RData (Samples with matched Clinical data and mRNA expression)

`TCGAcurated.R`\
Functions for data extraction from the curated datasets (_TCGA_curated_dataset.*.RData_)

* _extract_curated_: 
Extract data from specific project for genes, miRNA and clinical data
* _extract_all_curated_: Extract arbitrary data from all the projects available in a single dataframe 
* _plotExprCancer_: Plot gene expression profile across all tumor projects
* _plotExprCancerVsNormal_: Plot gene expression profile across paired tumor vs normal samples, when available


## Files available for download upon request:

md5sum                            |  File
--------------------------------- | ----------------------------------------
51d65a22b4f8c5772d97338ec4df8091  |	 enzymes.txt
d647a96306d6c7ba505b612cc64a41e4	|	 TCGA_curated_dataset.RData
fd874ed9e5ee959262f595d31bbd62a3	|  TCGA_curated_dataset.mRNAonly.RData
