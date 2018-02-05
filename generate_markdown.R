#args = commandArgs(trailingOnly=TRUE)

setwd("~/LAB/TCGA/CancerMods/")
library(rmarkdown);
source("TCGAcurated.mRNAonly.R")
library(stringr)
library(bsselectR)
library(knitr)
library(magrittr)
library(DT)
library(plyr)

generateGenePage<-function(ens_id, gof_id, descr){
	rmarkdown::render("geneTab.Rmd", output_file=paste0(ens_id, ".html"), quiet=T, intermediates_dir=paste0("imd_",ens_id), params = list(gene = ens_id, gene_official= gof_id, gene_descr=descr))
	system(paste0("mkdir ", ens_id))
	system(paste0("mv ", ens_id, "*.png ", ens_id))
	system(paste0("mv ", ens_id, "*.html ", ens_id))
	system(paste0("rm -rf imd_",ens_id))
}

#generateGenePage("ENSG00000165819", "METTL3", "Some description of the gene")
#generateGenePage("ENSG00000037897", "METTL1", "Some description of the gene")


enzymes<-read.csv("enzymes.txt", sep="\t", header=T, stringsAsFactors=F)


##############################################################
## Parallel implementation
##############################################################

library(foreach)
library(doParallel)
cores=detectCores()

registerDoParallel(cores=14)
getDoParWorkers()

foreach(i=1:nrow(enzymes)) %dopar% {

	ensid_proc <- enzymes[i, "ENSEMBL_ID"]
	gene_proc <- enzymes[i, "Gene_name"]
	desc_proc <- enzymes[i, "Description"]
	if(!file.exists(paste0(ensid_proc, "/", ensid_proc, ".html"))){
		#message(ensid_proc)
		generateGenePage(ensid_proc, gene_proc, desc_proc)
	}
}

stopCluster(cl)

rmarkdown::render("index.Rmd", output_file="index.html")
system("rm -rf imd_*")

#ls ENSG00000* | grep .html | wc -l #check status of the make process

#mv ENSG00000* site/
#tar -cvf site.tar site/
#pigz -p 15 site.tar
#scp -r site.tar.gz lp471@cb-head2:/mnt/beegfs/NGS/WebServers/Pipeline/htdocs/data/Luca/


##############################################################
## Chunk processing
##############################################################
#
#a=as.numeric(args[1])
#b=as.numeric(args[2])
#
#for(i in seq(a,b,1)){
#
#	ensid_proc <- enzymes[i, "ENSEMBL_ID"]
#	gene_proc <- enzymes[i, "Gene_name"]
#	desc_proc <- enzymes[i, "Description"]
#	
#	if(!file.exists(paste0(ensid_proc, "/", ensid_proc, ".html"))){
#		generateGenePage(ensid_proc, gene_proc, desc_proc)
#	}
#}
#
#rmarkdown::render("index.Rmd", output_file="index.html")
#
## Rscript generate_markdown.R 1 35 2>&1 | tee A.log & 
## Rscript generate_markdown.R 36 71 2>&1 | tee B.log & 
## Rscript generate_markdown.R 72 107 2>&1 | tee C.log & 
## Rscript generate_markdown.R 108 143 2>&1 | tee D.log & 
## Rscript generate_markdown.R 144 179 2>&1 | tee E.log & 
## Rscript generate_markdown.R 180 215 2>&1 | tee F.log &
## Rscript generate_markdown.R 216 251 2>&1 | tee G.log &
## Rscript generate_markdown.R 252 287 2>&1 | tee H.log &
## Rscript generate_markdown.R 288 310 2>&1 | tee I.log &
