library("TCGAbiolinks")
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(plyr)
library(MASS)
setwd("~/LAB/TCGA/CancerMods/")

message("Loading curated TCGA dataset (mRNA 0nly)...")
load("TCGA_curated_dataset.mRNAonly.RData")
options(stringsAsFactors = TRUE)

#############################################################################################################

## FUNCTIONS FOR DATA EXTRACTION FROM CURATED DATASETS

#curated_mRNA_dataset[["TCGA-GBM"]][["Clinical"]][1:5,1:5]

#curated_mRNA_dataset[["TARGET-NBL"]][["mRNA"]][1:5,1:5]
#curated_mRNA_dataset[["TARGET-AML"]][["Clinical"]][1:5,1:5]

#proj="TCGA-OV"
#genes="ENSG00000165819"
#clinics=c("tumor_stage","days_to_death")

#############################################################################################################
## EXTRACT RELEVANT INFORMATION FROM CURATED DATASETS OF SPECIFIC PROJECT ###################################
#############################################################################################################

extract_curated <-function(proj, genes, clinics){

	mRNA_xp <- curated_mRNA_dataset[[proj]][["mRNA"]][genes,]
	if(length(genes)>1){identifiers <- colnames(mRNA_xp)}else{identifiers <- names(mRNA_xp)}
	
	if ("shortLetterCode" %in% rownames(curated_mRNA_dataset[[proj]][["Clinical"]])){
		clin<- curated_mRNA_dataset[[proj]][["Clinical"]][c("shortLetterCode", clinics), identifiers]
	} else {
		clin <- rbind(rep("TP", length(identifiers)), curated_mRNA_dataset[[proj]][["Clinical"]][clinics, identifiers])
		rownames(clin) <- c("shortLetterCode", clinics)
	}
	
	extracted <- rbind(rep(proj, length(identifiers)), mRNA_xp, clin)
	rownames(extracted)[1]<-"project"
	if(length(genes)<2){rownames(extracted)[2]<-genes}

	extracted<-t(extracted)

	extracted2<-as.data.frame(matrix(unlist(extracted), nrow=nrow(extracted), ncol=ncol(extracted), byrow=FALSE))

	colnames(extracted2)<-colnames(extracted)
	rownames(extracted2)<-rownames(extracted)
	numeric_columns <- c(grep("ENSG", colnames(extracted2)), grep("days_to_death", colnames(extracted2)))

	for(nc in numeric_columns){
		extracted2[,nc]<-as.numeric(as.vector(extracted2[,nc]))
	}

	return(extracted2)

}

#head(extract_curated(proj="TCGA-OV", genes="ENSG00000165819", clinics=c("days_to_death","days_to_last_follow_up","vital_status")))

#############################################################################################################
## EXTRACT RELEVANT INFORMATION FROM CURATED DATASETS OF ALL CURATED PROJECT ################################ 
#############################################################################################################

extract_all_curated<-function(genes, clinics){

	all_extracted<-c()
	for (i in names(curated_mRNA_dataset)){
		
		#message(i)
		all_extracted<-rbind(all_extracted, extract_curated(i, genes, clinics))

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

#all <- extract_all_curated(genes=c("ENSG00000165819", "ENSG00000145388"), clinics=c("tumor_stage","days_to_death"))

#ggplot(all, aes(x=`ENSG00000165819`, y=`ENSG00000145388`)) + geom_point(aes(colour = project)) + scale_x_continuous(trans="log2") + scale_y_continuous(trans="log2")


#############################################################################################################
## PLOT GENE EXPRESSION PROFILE ACROSS ALL TUMOR DATA ####################################################### 
#############################################################################################################

plotExprCancer<-function(gene, max=0){	

	all <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))
	all_oncologic <- subset(all, oncologic)
	p <- ggplot(all_oncologic, aes(x=project, y=all_oncologic[,gene], fill=project)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1), legend.position="none") + geom_hline(yintercept=mean(all_oncologic[,gene]),linetype="dashed", color = "red") + ylab(paste0(gene, " (RPKM)"))

	if(max != 0){ print(p + ylim(0, max)) } else {print(p)}
	
	table(droplevels(all_oncologic[,"project"]))

}

#plotExprCancer("ENSG00000165819", max=30)

#############################################################################################################
## PLOT GENE EXPRESSION PROFILE ACROSS PAIRED TUMOR vs NORMAL SAMPLES ####################################### 
#############################################################################################################

plotExprCancerVsNormal<-function(gene, max=0){	

	all <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))

	all_with_normaldata <- subset(all, !oncologic)
	all2<-all[as.vector(all$project) %in% levels(droplevels(all_with_normaldata$project)),]
	
	p <- ggplot(all2, aes(x=project, y=all2[,gene], fill=oncologic)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) + ylab(paste0(gene, " (RPKM)"))

	if(max != 0){ print(p + ylim(0, max)) } else {print(p)}

	table(droplevels(all2[,c("project","oncologic")]))

}

#plotExprCancerVsNormal("ENSG00000165819")

#############################################################################################################
## PLOT GENE EXPRESSION PROFILE ACROSS PAIRED PRIMARY TUMOR vs METASTASIS SAMPLES ########################### 
#############################################################################################################

plotExprCancerVsMetastasis<-function(gene, max=0){	

	all <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))

	all_with_metastasis <- subset(all, metastasis)
	all2<-all[as.vector(all$project) %in% levels(droplevels(all_with_metastasis$project)),]
	all2<-subset(all2, oncologic)
	
	p <- ggplot(all2, aes(x=project, y=all2[,gene], fill=metastasis)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab(paste0(gene, " (RPKM)"))

	if(max != 0){ print(p + ylim(0, max)) } else {print(p)}

	table(droplevels(all2[,c("project","metastasis")]))

}

plotAllExpr<-function(gene, max=0){

	filename=paste0("All_Exp_", gene, ".pdf")
	pdf(file=filename, width=8, height=5)
	message("Plot ExprCancer")
	plotExprCancer(gene, max)
	message("Plot ExprCancerVsNormal")
	plotExprCancerVsNormal(gene, max)
	message("Plot ExprCancerVsMetastasis")
	plotExprCancerVsMetastasis(gene, max)
	dev.off()
}

#plotAllExpr("ENSG00000165819")

#############################################################################################################
## PLOT KAPLAN MEIER SURVIVAL CURVES LAYERED FOR GENE EXPRESSION ACROSS EACH PROJECT ######################## 
#############################################################################################################

library(survival)
library(survminer)

## Check the presence of the required clinical data in all the datasets
#KM_datasets<-c()
#
#for( n in names(curated_mRNA_dataset)){
#	cat(paste0(c(n,"... ")))
# if(sum(c("days_to_death","days_to_last_follow_up","vital_status") %in% rownames(curated_mRNA_dataset[[n]][["Clinical"]])) == 3){
#	KM_datasets <- c(KM_datasets, n)
#	message(" YES")
#	} else { message(" NO") }
#}

p2sym <- function(x){return(as.character(symnum(x, corr = FALSE,
                    cutpoints = c(0,  .001,.01,.05,  1),
                    symbols = c("***","**","*"," "))[[1]]))}

doKaplanMeier<-function(cfu, name=NULL, stoptime_analysis=NULL, stoptime_plot=NULL, make_plot=T){

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
    km_pval <- surv_pvalue(fit2, data=cfu)$pval
    num <- fit2$n
    fit_cox <- as.vector(coxph(Surv(ttime, status) ~ expr_code, data=cfu)$coefficients)

    if(make_plot){
       if(!is.null(stoptime_plot)){	
   	   print(ggsurvplot(fit2, data = cfu, risk.table = TRUE, xlim = c(0, stoptime_plot), pval = TRUE, legend.title = name)) #, conf.int = TRUE
       } else {
  	   print(ggsurvplot(fit2, data = cfu, risk.table = TRUE, pval = TRUE, legend.title = name)) #, conf.int = TRUE
       }
    }
    return(c(num[1],fit_cox,km_pval,p2sym(km_pval)))
    
}

#project="TCGA-BRCA"
#gene="ENSG00000165819"
#percentile=0.25

geneSurvival<-function(project, gene, percentile, stoptime_analysis=NULL, stoptime_plot=NULL, make_plot=T){

	xc<-extract_curated(proj=project, genes=gene, clinics=c("days_to_death","days_to_last_follow_up","vital_status"))
	xc<-xc[! as.vector(xc$shortLetterCode) == "NT",] #remove non-tumour samples
	lowq_exp=as.vector(quantile(xc[,2],  probs=percentile))
	hiq_exp=as.vector(quantile(xc[,2],  probs=1-percentile))
	xc$expr_code <- rep(NA, nrow(xc))
	xc[xc[,2] > hiq_exp, "expr_code"]<-"HIGH"
	xc[xc[,2] <= lowq_exp, "expr_code"]<-"LOW"
	xc <- xc[!is.na(xc$expr_code),]
	
	doKaplanMeier(xc, project, stoptime_analysis, stoptime_plot, make_plot)

}

#geneSurvival("TCGA-BRCA", "ENSG00000165819", 0.25)
#geneSurvival("TCGA-BRCA", "ENSG00000165819", 0.25, stoptime_plot=2000)
#geneSurvival("TCGA-BRCA", "ENSG00000165819", 0.25, stoptime_analysis=1600)

generateAllKMs<-function(gene, percentile, stoptime_analysis=NULL, stoptime_plot=NULL, make_plot=T, messages=T){

	stats<-c()
	dnames<-names(curated_mRNA_dataset)
	dnames<-dnames[(!dnames %in% c("TARGET-NBL"))]	# !! This dataset has no days_to_follow_up information!! 

	filename=paste0("KaplanMeier_", gene, "_", percentile, ".pdf")

	if(make_plot){pdf(file=filename, width=5, height=5)}
	for (pj in dnames){

			if(messages){message(pj)}
			suppressWarnings(this_line<-rbind(geneSurvival(pj, gene, percentile, stoptime_analysis=NULL, stoptime_plot=NULL, make_plot)))
			rownames(this_line)<-pj
			stats<-rbind(stats,this_line)
	}
	if(make_plot){dev.off()}
	colnames(stats)<-c("N_surv","cox_coef","pval","signif")
	stats<-as.data.frame(stats, stringsAsFactors=F)
	stats$N_surv<-as.numeric(stats$N_surv)
	stats$cox_coef<-as.numeric(stats$cox_coef)
	stats$pval<-as.numeric(stats$pval)
	stats$col<-"black"
	stats[stats$cox_coef < 0 & stats$pval < 0.05 , "col"]<-"red"
	stats[stats$cox_coef > 0 & stats$pval < 0.05 , "col"]<-"dodgerblue"
	return(stats)
}

#foo<-generateAllKMs("ENSG00000165819", 0.25)

plotExprCancer_surv<-function(gene, percentile, title=""){	
	# adds color to project label if the genes is a associated with better (green) or worse (red) diagnosis
	stats<-generateAllKMs(gene, percentile, make_plot=F, messages=F)
	all <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))
	all_oncologic <- subset(all, oncologic)
	colnames(all_oncologic)[which(colnames(all_oncologic) == gene)]<-"Expr"
	pj_medians<-ddply(all_oncologic, "project", summarise, N = length(Expr), Expr = median(Expr))
	rownames(pj_medians)<-pj_medians$project
	
	pj_medians$glob_medmean <- median(pj_medians$Expr)
	pj_medians$glob_medsd <- sd(pj_medians$Expr)
	pj_medians$log2_delta <- log(pj_medians$Expr / pj_medians$glob_medmean,2)
	
	 ga <- tryCatch(	## Tries to fit a gamma distribution. If it doesn't converge, it fits a normal distribution
        	{
			fitdistr(pj_medians$Expr,"gamma")
		},
		error=function(cond) {
			return("Exception")
		}, silent=T)
		
	if(ga != "Exception"){
		up_CI <- qgamma(0.95, shape=as.vector(ga$estimate["shape"]), rate=as.vector(ga$estimate["rate"]))
		low_CI <-qgamma(0.05, shape=as.vector(ga$estimate["shape"]), rate=as.vector(ga$estimate["rate"]))
		pj_medians$pval_diff <- pgamma(pj_medians$Expr, shape=as.vector(ga$estimate["shape"]), rate=as.vector(ga$estimate["rate"]))
	} else {
		except_mean <- mean(pj_medians$Expr)
		except_sd <- sd(pj_medians$Expr)
		pj_medians$pval_diff <- pnorm(pj_medians$Expr, mean=except_mean, sd=except_sd)
		up_CI <- qnorm(0.95, mean=except_mean, sd=except_sd)
		low_CI <- qnorm(0.05, mean=except_mean, sd=except_sd)
	}
	mask <- pj_medians$log2_delta>0 
	mask[is.na(mask)] <- F
	pj_medians[mask ,"pval_diff"]<- 1 - pj_medians[mask ,"pval_diff"]
	pj_medians$diff_signif<-""
	for (i in rownames(pj_medians)){
		pj_medians[i,"diff_signif"]<- p2sym(pj_medians[i,"pval_diff"])
		}

	lab_col<-stats[levels(all_oncologic$project),"col"]
	lab_col[is.na(lab_col)]<-"black"
	stats<-cbind(pj_medians, stats[rownames(pj_medians),])
	stats$signif[is.na(stats$signif)]<-" "
	stats$Expr <- signif(stats$Expr, 3)
	stats$log2_delta <- signif(stats$log2_delta, 3)
	stats$pval_diff <- signif(stats$pval_diff, 3)
	stats$cox_coef <- signif(stats$cox_coef, 3)
	stats$pval <- signif(stats$pval, 3)
	p <- ggplot(all_oncologic, aes(x=project, y=all_oncologic$Expr, fill=project)) +  geom_hline(yintercept=c(up_CI, low_CI), linetype="dashed", color = "gray") + geom_boxplot() + geom_hline(yintercept=pj_medians$glob_medmean[1], linetype="dashed", color = "black") + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, colour = lab_col), legend.position="none")  + ylab(paste0(gene, " (RPKM)")) + labs(title=title)

	print(p + ylim(0, 3*max(stats$Expr)))
	return(stats)
}

# plotExprCancer_surv("ENSG00000165819", 0.25)
# plotExprCancer_surv("ENSG00000037897", 0.25, "METTL1")

plotExprCancerVsNormal_stats<-function(gene, title){	

	allCN <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))

	all_with_normaldata <- subset(allCN, !oncologic)
	all2<-allCN[as.vector(allCN$project) %in% levels(droplevels(all_with_normaldata$project)),]
	colnames(all2)[which(colnames(all2) == gene)]<-"Expr"
	
	CvNorStat<-c()
	pjctnames <- levels(droplevels(all2$project))
	for(pjct in pjctnames){
	#message(pjct)
		tmp_set_oncologic<-subset(all2, project==pjct & oncologic)
		tmp_set_normal<-subset(all2, project==pjct & !oncologic)
		median_a<-median(tmp_set_oncologic$Expr)
		median_b<-median(tmp_set_normal$Expr)
		pval_on<-wilcox.test(tmp_set_oncologic$Expr, tmp_set_normal$Expr)$p.value
		CvNorStat<-rbind(CvNorStat, c(nrow(tmp_set_normal), nrow(tmp_set_oncologic), median_b, median_a, log((median_a/median_b),2), pval_on, p2sym(pval_on)))
	}
	CvNorStat<-as.data.frame(CvNorStat, stringsAsFactors=F)
	CvNorStat$V1<-signif(as.numeric(CvNorStat$V1),3)
	CvNorStat$V2<-signif(as.numeric(CvNorStat$V2),3)
	CvNorStat$V3<-signif(as.numeric(CvNorStat$V3),3)
	CvNorStat$V4<-signif(as.numeric(CvNorStat$V4),3)
	CvNorStat$V5<-signif(as.numeric(CvNorStat$V5),3)
	CvNorStat$V6<-signif(as.numeric(CvNorStat$V6),3)
	colnames(CvNorStat)<-c("nNorm", "nCancer","ExprNorm","ExprCancer","log2Delta","pval"," ")
	rownames(CvNorStat)<-pjctnames
	
	CvNorStat$col<-"#000000"
	CvNorStat$log2Delta[is.na(CvNorStat$log2Delta)]<-0
	CvNorStat[CvNorStat$log2Delta < 0 & CvNorStat$pval < 0.05 , "col"]<-"#f8766d"
	CvNorStat[CvNorStat$log2Delta > 0 & CvNorStat$pval < 0.05 , "col"]<-"#00bfc4"

	all2<-droplevels(all2)
	
	p <- ggplot(all2, aes(x=project, y=Expr, fill=oncologic)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, colour = CvNorStat$col)) + ylab(paste0(gene, " (RPKM)")) + guides(fill=guide_legend(title="Tissue"))  + scale_fill_discrete(breaks=c(T, F), labels=c("Cancer", "Normal")) + labs(title=title)

	print(p + ylim(0, 3*max(c(max(CvNorStat$ExprNorm),max(CvNorStat$ExprCancer)))))

	CvNorStat$col<-NULL
	CvNorStat<-cbind(rownames(CvNorStat),CvNorStat)
	colnames(CvNorStat)[1]<-"TCGA_project_ID"
	return(CvNorStat)

}

#plotExprCancerVsNormal_stats("ENSG00000165819", "METTL3")

plotExprCancerVsMetastasis_stats<-function(gene, title){	

	allCM <- extract_all_curated(genes=gene, clinics=c("tumor_stage","days_to_death"))

	all_with_metastasis <- subset(allCM, metastasis)
	all3<-allCM[as.vector(allCM$project) %in% levels(droplevels(all_with_metastasis$project)),]
	all3<-subset(all3, oncologic)
	
	colnames(all3)[which(colnames(all3) == gene)]<-"Expr"
	
	CvMetStat<-c()
	pjctnames <- levels(droplevels(all3$project))
	for(pjct in pjctnames){
	#message(pjct)
		tmp_set_metastasis<-subset(all3, project==pjct & metastasis)
		tmp_set_primary<-subset(all3, project==pjct & !metastasis)
		median_a<-median(tmp_set_metastasis$Expr)
		median_b<-median(tmp_set_primary$Expr)
		pval_on<-wilcox.test(tmp_set_metastasis$Expr, tmp_set_primary$Expr)$p.value
		CvMetStat<-rbind(CvMetStat, c(nrow(tmp_set_primary), nrow(tmp_set_metastasis), median_b, median_a, log((median_a/median_b),2), pval_on, p2sym(pval_on)))
	}
	CvMetStat<-as.data.frame(CvMetStat, stringsAsFactors=F)
	CvMetStat$V1<-signif(as.numeric(CvMetStat$V1),3)
	CvMetStat$V2<-signif(as.numeric(CvMetStat$V2),3)
	CvMetStat$V3<-signif(as.numeric(CvMetStat$V3),3)
	CvMetStat$V4<-signif(as.numeric(CvMetStat$V4),3)
	CvMetStat$V5<-signif(as.numeric(CvMetStat$V5),3)
	CvMetStat$V6<-signif(as.numeric(CvMetStat$V6),3)
	colnames(CvMetStat)<-c("nPrim", "nMetas","ExprPrim","ExprMetas","log2Delta","pval"," ")
	rownames(CvMetStat)<-pjctnames
	
	CvMetStat$col<-"#000000"
	CvMetStat$log2Delta[is.na(CvMetStat$log2Delta)]<-0
	CvMetStat[CvMetStat$log2Delta < 0 & CvMetStat$pval < 0.05 , "col"]<-"#f8766d"
	CvMetStat[CvMetStat$log2Delta > 0 & CvMetStat$pval < 0.05 , "col"]<-"#00bfc4"

	all3<-droplevels(all3)
	
	p <- ggplot(all3, aes(x=project, y=Expr, fill=metastasis)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, colour = CvMetStat$col)) + ylab(paste0(gene, " (RPKM)")) + guides(fill=guide_legend(title="Tissue")) + labs(title=title) + scale_fill_discrete(breaks=c(F, T), labels=c("Primary", "Metas"))

	print(p + ylim(0, 3*max(c(max(CvMetStat$ExprPrim),max(CvMetStat$ExprMetas)))))

	CvMetStat$col<-NULL
	CvMetStat<-cbind(rownames(CvMetStat),CvMetStat)
	colnames(CvMetStat)[1]<-"TCGA_project_ID"
	return(CvMetStat)

}

#plotExprCancerVsMetastasis_stats("ENSG00000148584", "METTL3")
