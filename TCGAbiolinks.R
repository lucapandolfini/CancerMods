library("TCGAbiolinks")
library(SummarizedExperiment)
library(doMC)
library(foreach)
registerDoMC(14)
library(reshape2)

system("mkdir ~/DATA/TCGAbiolinks_dir")
setwd("~/DATA/TCGAbiolinks_dir")

all_projects <- TCGAbiolinks:::getGDCprojects()$project_id

exprs_projects<-c()
for (i in all_projects){
	message(i)
	data_cat<-TCGAbiolinks:::getProjectSummary(i)$data_categories[,3]
	if(sum(which(data_cat == "Transcriptome Profiling")) > 0 ){
		exprs_projects<-c(exprs_projects,i)
	}
}

mRNA_query<-list()
microRNA_query<-list()
microRNAiso_query<-list()
clinical_data<-list()
query1<-NULL
query2<-NULL
query3<-NULL

system("touch tcga.log")

for (i in exprs_projects){
	system('echo "#################################### TCGA Dataset ################################" >> tcga.log')
	system(paste0('echo ',i,' >> tcga.log'))
	system('echo "##################################################################################" >> tcga.log')
	system('echo "Get mRNA query..." >> tcga.log')
	tryCatch({ 
	query1 <- GDCquery(project = i,
           data.category = "Transcriptome Profiling",
           data.type = "Gene Expression Quantification", 
           workflow.type = "HTSeq - FPKM",
	   legacy = FALSE)
	}, error=function(e){system('echo "ERROR : mRNA dataset not available" >> tcga.log')})

	if(!is.null(query1)){
		mRNA_query[[i]] <- query1
		system('echo "Download mRNA datasets..." >> tcga.log')
		GDCdownload(query1, method = "client")
	}
		
fsystem('echo "Get microRNA query..." >> tcga.log')
	tryCatch({ 
	query2 <- GDCquery(project = i,
        data.category = "Transcriptome Profiling",
           data.type = "miRNA Expression Quantification", 
	   legacy = FALSE)
	}, error=function(e){system('echo "ERROR : microRNA dataset not available" >> tcga.log')})

	if(!is.null(query2)){
	 	microRNA_query[[i]] <- query2
		system('echo "Download microRNA datasets..." >> tcga.log')
		GDCdownload(query2, method = "client")
	}

	system('echo "Get microRNA isoform query..." >> tcga.log')
	tryCatch({ 
	query3 <- GDCquery(project = i,
        data.category = "Transcriptome Profiling",
           data.type = "Isoform Expression Quantification", 
	   legacy = FALSE)
	}, error=function(e){system('echo "ERROR : microRNA isoform dataset not available" >> tcga.log')})

	if(!is.null(query3)){
	 	microRNAiso_query[[i]] <- query3
		system('echo "Download microRNA isoform datasets..." >> tcga.log')
		GDCdownload(query3, method = "client")
	}


	clinical_data[[i]] <- GDCquery_clinic(project = i, type = "clinical")

	query1<-NULL
	query2<-NULL
	query3<-NULL

}

#save.image(file="queries.RData")

##########################################################################################################################

# Fixing the shit of TCGA misannotation.
# Manual tweak for the "Warning: There are more than one file for the same case"
# causing GDCprepare to fail for some projects

investigate<-microRNA_query[["TARGET-RT"]]
investigate$results[[1]] <- investigate$results[[1]][c(1:36,38:48,50),]
microRNA_query[["TARGET-RT"]]<-investigate

investigate<-microRNAiso_query[["TARGET-RT"]]
investigate$results[[1]] <- investigate$results[[1]][c(1:7,9:18,20:50),]
microRNAiso_query[["TARGET-RT"]]<-investigate

for (pj in names(mRNA_query)){
	mRNA_query[[pj]]$results[[1]]$case
}

mRNA_data<-foreach(i=names(mRNA_query)) %dopar% {
	return(GDCprepare(mRNA_query[[i]]))
}

#rm(i, query1, query2, data_cat, investigate)
#save.image(file="data.RData")

## Access mRNA_expr_list through project_id
mRNA_expr_data<-list()
for(n in 1:length(mRNA_data)){

	pj_id <- colData(mRNA_data[[n]])$project_id[1]
	mRNA_expr_data[[pj_id]] <- mRNA_data[[n]]

}
mRNA_data<-NULL

microRNA_expr_data<-list()
for(i in names(microRNA_query)){
	
	message(i)
	foo<-GDCprepare(microRNA_query[[i]])
	samples<-colnames(foo)
	sample_num<-c(1, grep("reads_per_million_miRNA_mapped_", samples))
	samples<-samples[sample_num]
	samples<-gsub("reads_per_million_miRNA_mapped_","",samples)
	foo<-foo[,sample_num]
	colnames(foo)<-samples

	microRNA_expr_data[[i]]<-foo
}


miID2Name<-read.table(file="hsa_miR_accessionTOname.txt", row.names=1, header=T)
microRNAiso_expr_data<-list()

microRNAiso_expr_data_tmp<-foreach(i=names(microRNAiso_query)) %dopar% {

	message(i)
	tmp<-GDCprepare(microRNAiso_query[[i]])

	tmp<-as.data.frame(tmp)
	tmp<-aggregate(data=tmp, reads_per_million_miRNA_mapped ~ miRNA_region + barcode, FUN=sum)
	tmp<-dcast(tmp, miRNA_region ~ barcode, value.var= "reads_per_million_miRNA_mapped")

	tmp<-tmp[grep("mature", tmp$miRNA_region),]
	tmp$miRNA_region<- gsub("mature,","", tmp$miRNA_region)

	tmp$miRNA_names<- as.vector(miID2Name[tmp$miRNA_region, "Name"])

	tmp<-tmp[!is.na(tmp$miRNA_names),]

	tmp$miRNA_region<-tmp$miRNA_names
	tmp$miRNA_names<-NULL
	colnames(tmp)[1]<-"miRNA_ID"

	return(list(i,tmp))
	tmp<-NULL

}

microRNAiso_expr_data<-list()
for(n in 1:length(microRNAiso_expr_data_tmp)){

	pj_id <- microRNAiso_expr_data_tmp[[n]][[1]]
	microRNAiso_expr_data[[pj_id]] <- microRNAiso_expr_data_tmp[[n]][[2]]

}

microRNAiso_expr_data_tmp<-NULL

#save.image("TCGA_expression_datasets.RData")

##########################################################################################################################
## Create a curate dataset with:
## - non duplicated patients
## - only projects with both mRNA and miRNA datasets available
## - 1:1 mapping between mRNA and miRNA data
##########################################################################################################################

miRid2mRNAid<-function(str){
	return(paste0(substr(str, 1,19),"R"))
}

TCGA_projects <- exprs_projects[grep("TCGA", exprs_projects)]
TARGET_projects <- exprs_projects[grep("TARGET", exprs_projects)]

curated_dataset<-list()

for(proj in TCGA_projects){

	message(proj)
	samples_mRNA<-colnames(mRNA_expr_data[[proj]])
	samples_microRNA<-colnames(microRNAiso_expr_data[[proj]])

	mirdata2 <- microRNAiso_expr_data[[proj]][,!(duplicated(miRid2mRNAid(samples_microRNA)))]
	colnames(mirdata2)<-miRid2mRNAid(colnames(mirdata2))
	expr_with_micro<-which(substr(samples_mRNA, 1,20) %in% miRid2mRNAid(samples_microRNA))
	mrnadata2 <- mRNA_expr_data[[proj]][,expr_with_micro]
	mrnadata2 <- mrnadata2[,!duplicated(substr(colnames(mrnadata2),1,20))]
	colnames(mrnadata2) <- substr(colnames(mrnadata2),1,20)
	
	rownames(mirdata2)<-mirdata2$miRNA_IDR
	mirdata2 <- mirdata2[,colnames(mrnadata2)]
	
	clinic<-colData(mRNA_expr_data[[proj]])
	clinic$InternalSampleID<-substr(rownames(clinic),1,20)
	clinic<-clinic[!duplicated(clinic$InternalSampleID),]
	rownames(clinic)<-clinic$InternalSampleID
	clinic<-clinic[colnames(mrnadata2),]
	clinic<-t(as.data.frame(clinic))
	
	mrnadata2<-assay(mrnadata2)

	data_list <- list(mrnadata2, mirdata2, clinic)
	names(data_list) <- c("mRNA","microRNA","Clinical")
	curated_dataset[[proj]] <- data_list

}

for(proj in c("TCGA-LAML",TARGET_projects[c(1,2,4)])){

	message(proj)
	samples_mRNA<-colnames(mRNA_expr_data[[proj]])
	samples_microRNA<-colnames(microRNAiso_expr_data[[proj]])

	mirdata2 <- microRNAiso_expr_data[[proj]][,!(duplicated(samples_microRNA))]
	expr_with_micro<-which(samples_mRNA %in% samples_microRNA)
	mrnadata2 <- mRNA_expr_data[[proj]][,expr_with_micro]
	mrnadata2 <- mrnadata2[,!duplicated(colnames(mrnadata2))]
	
	rownames(mirdata2)<-mirdata2$miRNA_ID
	mirdata2 <- mirdata2[,colnames(mrnadata2)]
	
	clinic<-colData(mRNA_expr_data[[proj]])
	clinic<-clinic[colnames(mrnadata2),]
	clinic<-t(as.data.frame(clinic))
	
	mrnadata2<-assay(mrnadata2)

	data_list <- list(mrnadata2, mirdata2, clinic)
	names(data_list) <- c("mRNA","microRNA","Clinical")
	curated_dataset[[proj]] <- data_list

}

#save(curated_dataset, file="TCGA_curated_dataset.RData")

#curated_dataset[["TCGA-OV"]][["mRNA"]][1:5,1:5]
#curated_dataset[["TCGA-OV"]][["microRNA"]][1:5,1:5]
#curated_dataset[["TCGA-OV"]][["Clinical"]][1:5,1:5]

#curated_dataset[["TARGET-NBL"]][["mRNA"]][1:5,1:5]
#curated_dataset[["TARGET-AML"]][["microRNA"]][1:5,1:5]
#curated_dataset[["TARGET-AML"]][["Clinical"]][1:5,1:5]

curated_mRNA_dataset<-list()

for(proj in exprs_projects){

	message(proj)
	samples_mRNA<-colnames(mRNA_expr_data[[proj]])

	mrnadata2 <- mRNA_expr_data[[proj]][,!(duplicated(samples_mRNA))]

	clinic<-colData(mRNA_expr_data[[proj]])
	clinic<-clinic[colnames(mrnadata2),]
	clinic<-t(as.data.frame(clinic))
	
	mrnadata2<-assay(mrnadata2)

	data_list <- list(mrnadata2, clinic)
	names(data_list) <- c("mRNA","Clinical")
	curated_mRNA_dataset[[proj]] <- data_list

}

#save(curated_mRNA_dataset, file="TCGA_curated_dataset.mRNAonly.RData")
