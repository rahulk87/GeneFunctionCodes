# This R script calculates the association between expression of gene of interest and shRNA of other genes.
# In this example script, gene of interest is APC. So, here it'll identify the essential genes with high APC expression (-ve difference) and low APC expression (+ve difference)
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
options(width=100)
rnaseq=read.delim("../data/breast_rnaseq_qn.txt",header=T,as.is=T,check.names=F)

# Set gene of interest
geneIndex=which(rnaseq$symbol=="APC")

geneExpr=unlist(rnaseq[geneIndex,-grep("gene_id|symbol|ensembl_id",colnames(rnaseq))])
geneValues=cbind.data.frame(cell_line=names(geneExpr),gene=geneExpr,stringsAsFactors=FALSE)

# Set entrez gened IDs of query shRNAs
genesOfInterest=c(595,1019,1021)
results=simem(screens=breast_screens,geneIds = genesOfInterest,covariate="gene",reagentWeights=hpWeights,annotationsPerCellLine=geneValues,inverseVarianceWeights=TRUE,signalProbWeights=TRUE,analyzeReagents=TRUE,parallelNodes=8)

# Genome-wide analysis, it takes around 3-4 hours on HPC
#results=simem(screens=breast_screens,covariate="gene",reagentWeights=hpWeights,annotationsPerCellLine=geneValues,inverseVarianceWeights=TRUE,signalProbWeights=TRUE,analyzeReagents=TRUE,parallelNodes=8)

# Writing output
result <- results$gene
write.csv(result,file="result")
