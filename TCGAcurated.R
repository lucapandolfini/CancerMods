library("TCGAbiolinks")
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
setwd("~/LAB/TCGA/")

message("Loading curated TCGA dataset...")
load("TCGA_curated_dataset.RData")
options(stringsAsFactors=T)

#############################################################################################################

## FUNCTIONS FOR DATA EXTRACTION FROM CURATED DATASETS

#curated_dataset[["TCGA-OV"]][["mRNA"]][1:5,1:5]
#curated_dataset[["TCGA-OV"]][["microRNA"]][1:5,1:5]
#curated_dataset[["TCGA-OV"]][["Clinical"]][1:5,1:5]

#curated_dataset[["TARGET-NBL"]][["mRNA"]][1:5,1:5]
#curated_dataset[["TARGET-AML"]][["microRNA"]][1:5,1:5]
#curated_dataset[["TARGET-AML"]][["Clinical"]][1:5,1:5]

#proj="TCGA-OV"
#genes="ENSG00000165819"
#mirnas="hsa-miR-17-5p"
#clinics=c("tumor_stage","days_to_death")

#############################################################################################################
## EXTRACT RELEVANT INFORMATION FROM CURATED DATASETS OF SPECIFIC PROJECT ###################################
#############################################################################################################

extract_curated <-function(proj, genes, mirnas=NULL, clinics){

	mRNA_xp <- curated_dataset[[proj]][["mRNA"]][genes,]
	if(length(genes)>1){identifiers <- colnames(mRNA_xp)}else{identifiers <- names(mRNA_xp)}
	
	if ("shortLetterCode" %in% rownames(curated_dataset[[proj]][["Clinical"]])){
		clin<- curated_dataset[[proj]][["Clinical"]][c("shortLetterCode", clinics), identifiers]
	} else {
		clin <- rbind(rep("TP", length(identifiers)), curated_dataset[[proj]][["Clinical"]][clinics, identifiers])
		rownames(clin) <- c("shortLetterCode", clinics)
	}

		if(length(mirnas)>0){
			microRNA_xp <- curated_dataset[[proj]][["microRNA"]][mirnas, identifiers]
			extracted <- rbind(rep(proj, length(identifiers)), mRNA_xp, microRNA_xp, clin)
		} else {
			extracted <- rbind(rep(proj, length(identifiers)), mRNA_xp, clin)
			}

	rownames(extracted)[1]<-"project"
	if(length(genes)<2){rownames(extracted)[2]<-genes}

	extracted<-t(extracted)

	extracted2<-as.data.frame(matrix(unlist(extracted), nrow=nrow(extracted), ncol=ncol(extracted), byrow=FALSE))

	colnames(extracted2)<-colnames(extracted)
	rownames(extracted2)<-rownames(extracted)
	numeric_columns <- c(grep("ENSG", colnames(extracted2)), grep("hsa-", colnames(extracted2)), grep("days_to_death", colnames(extracted2)))

	for(nc in numeric_columns){
		extracted2[,nc]<-as.numeric(as.vector(extracted2[,nc]))
	}

	return(extracted2)

}

#head(extract_curated(proj="TCGA-OV", genes="ENSG00000037897", mirnas="hsa-let-7e", clinics=c("days_to_death","days_to_last_follow_up","vital_status")))

#############################################################################################################
## EXTRACT RELEVANT INFORMATION FROM CURATED DATASETS OF ALL CURATED PROJECT ################################ 
#############################################################################################################

extract_all_curated<-function(genes, mirnas=NULL, clinics){

	all_extracted<-c()
	for (i in names(curated_dataset)){
		
		#message(i)
		all_extracted<-rbind(all_extracted, extract_curated(i, genes, mirnas, clinics))

	}

	all_extracted$oncologic<-all_extracted$shortLetterCode
	all_extracted$solid<-all_extracted$shortLetterCode
	all_extracted$metastasis<-all_extracted$shortLetterCode
	#> levels(all$shortLetterCode)
	#[1] "TP"  "NT"  "TR"  "TM"  "TAM" "TAP" "TB" 
	levels(all_extracted$oncologic) <-  c(T,F,T,T,T,T,T)
	levels(all_extracted$solid) <-      c(T,T,T,T,T,T,F)
	levels(all_extracted$metastasis) <- c(F,F,F,T,T,F,F)
	all_extracted$oncologic <- as.logical(all_extracted$oncologic)
	all_extracted$solid <- as.logical(all_extracted$solid)
	all_extracted$metastasis <- as.logical(all_extracted$metastasis)

	return(all_extracted)

}

#all <- extract_all_curated(genes=c("ENSG00000165819"), mirnas=c("hsa-miR-17-5p","hsa-miR-17-3p"), clinics=c("tumor_stage","days_to_death"))

#ggplot(all, aes(x=`hsa-miR-17-5p`, y=`hsa-miR-17-3p`)) + geom_point(aes(colour = project)) + scale_x_continuous(trans="log2") + scale_y_continuous(trans="log2")


#############################################################################################################
## PLOT GENE EXPRESSION PROFILE ACROSS ALL TUMOR DATA ####################################################### 
#############################################################################################################

plotExprCancer<-function(gene, max=0){	

	all <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))
	all_oncologic <- subset(all, oncologic)
	p <- ggplot(all_oncologic, aes(x=project, y=all_oncologic[,gene], fill=project)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") + geom_hline(yintercept=mean(all_oncologic[,gene]),linetype="dashed", color = "red") + ylab(gene)

	if(max != 0){ p + ylim(0, max) } else {p}
	
	table(droplevels(all_oncologic[,"project"]))

}

#############################################################################################################
## PLOT GENE EXPRESSION PROFILE ACROSS PAIRED TUMOR vs NORMAL SAMPLES ####################################### 
#############################################################################################################

plotExprCancerVsNormal<-function(gene, max=0){	

	all <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))

	all_with_normaldata <- subset(all, !oncologic)
	all2<-all[as.vector(all$project) %in% levels(droplevels(all_with_normaldata$project)),]
	
	p <- ggplot(all2, aes(x=project, y=all2[,gene], fill=oncologic)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

	if(max != 0){ p + ylim(0, max) } else {p}

	table(droplevels(all2[,c("project","oncologic")]))

}


#############################################################################################################
## PLOT GENE EXPRESSION PROFILE ACROSS PAIRED PRIMARY TUMOR vs METASTASIS SAMPLES ########################### 
#############################################################################################################

plotExprCancerVsMetastasis<-function(gene, max=0){	

	all <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))

	all_with_metastasis <- subset(all, metastasis)
	all2<-all[as.vector(all$project) %in% levels(droplevels(all_with_metastasis$project)),]
	all2<-subset(all2, oncologic)
	
	p <- ggplot(all2, aes(x=project, y=all2[,gene], fill=metastasis)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

	if(max != 0){ p + ylim(0, max) } else {p}

	table(droplevels(all2[,c("project","metastasis")]))

}


#############################################################################################################
## PLOT KAPLAN MEIER SURVIVAL CURVES LAYERED FOR GENE EXPRESSION ACROSS EACH PROJECT ######################## 
#############################################################################################################

library(survival)
library(survminer)

## Check the presence of the required clinical data in all the datasets
#KM_datasets<-c()
#
#for( n in names(curated_dataset)){
#	cat(paste0(c(n,"... ")))
# if(sum(c("days_to_death","days_to_last_follow_up","vital_status") %in% rownames(curated_dataset[[n]][["Clinical"]])) == 3){
#	KM_datasets <- c(KM_datasets, n)
#	message(" YES")
#	} else { message(" NO") }
#}

doKaplanMeier<-function(cfu, name=NULL, stoptime_analysis=NULL, stoptime_plot=NULL){

    cfu$vital_status <- as.vector(cfu$vital_status)
    cfu$days_to_last_follow_up <- as.numeric(as.vector(cfu$days_to_last_follow_up))

   # Set alive death to inf
    if(length(grep("alive",cfu$vital_status,ignore.case = TRUE)) > 0) cfu[grep("alive",cfu$vital_status,ignore.case = TRUE),"days_to_death"]<-(-Inf)

    # Set dead follow up to inf
    if(length(grep("dead",cfu$vital_status,ignore.case = TRUE)) > 0) cfu[grep("dead",cfu$vital_status,ignore.case = TRUE),"days_to_last_follow_up"]<-(-Inf)

    cfu <- cfu[ !(is.na(cfu[,"days_to_last_follow_up"])),]
    cfu <- cfu[ !(is.na(cfu[,"days_to_death"])),]

    cfu$ttime <- as.numeric(cfu[, "days_to_death"])

    cfu$status <- ( cfu$ttime > 0) 

    cfu[!cfu$status, "ttime"] <- as.numeric(cfu[!cfu$status, "days_to_last_follow_up"])

    if(!is.null(stoptime_analysis)){	
	cfu <- subset(cfu, ttime <= stoptime_analysis)
    } 

    cfu<-cfu[is.finite(cfu$ttime),]
    fit2 <- survfit(Surv(ttime, status) ~ expr_code, data=cfu)

    if(!is.null(stoptime_plot)){	
   	 ggsurvplot(fit2, data = cfu, risk.table = TRUE, xlim = c(0, stoptime_plot), pval = TRUE, legend.title = name) #, conf.int = TRUE
    } else {
  	 ggsurvplot(fit2, data = cfu, risk.table = TRUE, pval = TRUE, legend.title = name) #, conf.int = TRUE
    }

}

#project="TCGA-BRCA"
#gene="ENSG00000165819"
#percentile=0.25

geneSurvival<-function(project, gene, percentile, stoptime_analysis=NULL, stoptime_plot=NULL){

	xc<-extract_curated(proj=project, genes=gene, clinics=c("days_to_death","days_to_last_follow_up","vital_status"))
	xc<-xc[!is.na(xc[,3]),]
	lowq_exp=as.vector(quantile(xc[,2],  probs=percentile))
	hiq_exp=as.vector(quantile(xc[,2],  probs=1-percentile))
	xc$expr_code <- rep(NA, nrow(xc))
	xc[xc[,2] > hiq_exp, "expr_code"]<-"HIGH"
	xc[xc[,2] < lowq_exp, "expr_code"]<-"LOW"
	xc <- xc[!is.na(xc$expr_code),]
	
	doKaplanMeier(xc, project, stoptime_analysis, stoptime_plot)

}

#geneSurvival("TCGA-BRCA", "ENSG00000165819", 0.25)
#geneSurvival("TCGA-BRCA", "ENSG00000165819", 0.25, stoptime_plot=2000)
#geneSurvival("TCGA-BRCA", "ENSG00000165819", 0.25, stoptime_analysis=1600)

microSurvival<-function(project, mirna, percentile, stoptime_analysis=NULL, stoptime_plot=NULL){

	xc<-extract_curated(proj=project, genes="ENSG00000165819", mirnas=mirna, clinics=c("days_to_death","days_to_last_follow_up","vital_status"))
	xc<-xc[!is.na(xc[,3]),]
	lowq_exp=as.vector(quantile(xc[,3],  probs=percentile))
	hiq_exp=as.vector(quantile(xc[,3],  probs=1-percentile))
	xc$expr_code <- rep(NA, nrow(xc))
	xc[xc[,3] > hiq_exp, "expr_code"]<-"HIGH"
	xc[xc[,3] < lowq_exp, "expr_code"]<-"LOW"
	xc <- xc[!is.na(xc$expr_code),]
	
	doKaplanMeier(xc, project, stoptime_analysis, stoptime_plot)

}

#microSurvival("TCGA-BRCA", "hsa-miR-17-3p", 0.25)
#microSurvival("TCGA-BRCA", "hsa-miR-92b-3p", 0.25, stoptime_analysis=3000)

gene="ENSG00000165819"
micro="hsa-miR-17-3p"
percentile=0.25
proj="TCGA-BRCA"

plotPairedBiKM<-function(proj, gene, micro, percentile, stoptime_analysis=NULL, stoptime_plot=NULL){
	p1<-geneSurvival(proj, gene, percentile, stoptime_analysis, stoptime_plot)
	p2<-microSurvival(proj, micro, percentile, stoptime_analysis, stoptime_plot)
	ggdraw() + draw_plot(p1[[1]], 0,0.25,0.5,0.75) + draw_plot(p1[[2]], 0,0,0.5,0.25) + draw_plot(p2[[1]], 0.5,0.25,0.5,0.75) + draw_plot(p2[[2]], 0.5,0,0.5,0.25)
}

#plotBiKM("TARGET-AML", "ENSG00000165819", "hsa-miR-17-3p", 0.25)

generateAllPairedKMs<-function(gene, micro, percentile, stoptime_analysis=NULL, stoptime_plot=NULL){

	filename=paste0("KM_", gene, "_", micro, "_", percentile, ".pdf")

	pdf(file=filename, width=8, height=5)
	for (pj in names(curated_dataset)){
		if(! pj %in% c("TCGA-GBM")){	# !! The TCGA dataset has only 4 samples with miRNA info
		message(pj)
		print(plotPairedBiKM(pj, gene, micro, percentile, stoptime_analysis=NULL, stoptime_plot=NULL))
		}
	}
	dev.off()
	
}

#generateAllPairedKMs("ENSG00000165819", "hsa-miR-17-3p", 0.25)
