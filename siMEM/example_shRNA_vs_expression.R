# This R script calculates the association between shNRA of a specifice gene of interest and expression of other genes.
# In this example script, gene of interest is PRKDC (Entres id: 5591). So, here it'll identify the association between the PRKDC shRNA and expression of genes provided in the list. -ve difference means PRKDC is essential for cell lines having high expression concerend gene and vice versa.
# URL: http://neellab.github.io/simem/

# Required R packages
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("blme"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("locfit"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("preprocessCore"))
suppressPackageStartupMessages(library("reshape"))

# Required codes and datasets
source("../codes/data_format_lib.R")
source("../codes/model_lib.R")
source("../codes/simem_lib.R")
load("../data/breast_screens_with_weights.eset")
hp=read.delim("../data/hairpin_annotations.txt",header=T,as.is=T,check.names=F)
hpWeights=hp[,c("trcn_id","gene_id","weight")]
rnaseq=read.delim("../data/breast_rnaseq_qn.txt",header=T,as.is=T,check.names=F)

# Provide the list of genes e.g. CGC etc.
list=read.csv("../data/list",sep="\t",header=FALSE, stringsAsFactors=FALSE)
	for (genes in list$V1)
	{
	x <- genes
	y <- tolower(x)
	#ERROR HANDLING 
	possibleError <- tryCatch({
	print(paste("Start Loop ", x ,sep=""))

	yIndex=which(rnaseq$symbol==x)
	yExpr=unlist(rnaseq[yIndex,-grep("gene_id|symbol|ensembl_id",colnames(rnaseq))])
	yValues=cbind.data.frame(cell_line=names(yExpr),y=yExpr,stringsAsFactors=FALSE)
	
	# Set entrez id of gene of interest
	genesOfInterest=c(5591)
	
	results=simem(screens=breast_screens,geneIds = genesOfInterest,covariate="y",reagentWeights=hpWeights,annotationsPerCellLine=yValues,inverseVarianceWeights=TRUE,signalProbWeights=TRUE,analyzeReagents=TRUE,parallelNodes=8)
	
	# Writing output
	result <- results$gene
	result1 <- cbind(x,result)
        write.csv(result1,file=paste("result1_",x,sep=""))
	}
	,
	error=function(e) {
	e
	print(paste("Oops! --> Error in Loop ",x,sep = ""))
	})
	if(inherits(possibleError, "error")) next
	print(paste("  End Loop ",x,sep = ""))
	}

	#### Concatenating the result files ####
        files <- list.files(pattern = "result1_")
        out<-do.call("rbind", lapply(files, read.csv, header = TRUE))
	write.csv(out,file="siMEM_result")

        #### Removing the result files ####
        do.call(file.remove,list(list.files(pattern = "result1_")))
