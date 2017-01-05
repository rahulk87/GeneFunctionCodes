# This R script calculates the 2 class context specific gene essentiality.
# This example script identify the genes essential for the survival of ERBB2+ breast cell lines.
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
hp = read.delim("../data/hairpin_annotations.txt", header=T, as.is=T, check.names=F)
hpWeights = hp[,c("trcn_id", "gene_id", "weight")]
subtypes = read.delim("../data/cell_line_subtypes.txt", header=T, as.is=T, check.names=F)

# Set the cell lines classification
status = subtypes[,c("cell_line", "subtype_three_receptor")]
status$class = ifelse(status$subtype_three_receptor == "ERBB2", "ERBB2", "other")

# Set entrez gened IDs of query shRNAs
genesOfInterest=c(595,1019,1021)
results = simem(screens = breast_screens,geneIds = genesOfInterest,covariate = "class",reagentWeights = hpWeights, annotationsPerCellLine = status, inverseVarianceWeights = TRUE, signalProbWeights = TRUE,analyzeReagents = TRUE,covariateFactorOrder = c("other", "ERBB2"), parallelNodes = 8)

# Genome-wide analysis, it takes around 3-4 hours on HPC
#results = simem(screens = breast_screens,covariate = "erbb2",reagentWeights = hpWeights, annotationsPerCellLine = status, inverseVarianceWeights = TRUE, signalProbWeights = TRUE,analyzeReagents = TRUE,covariateFactorOrder = c("other", "ERBB2"), parallelNodes = 8)

# Writing output
result <- results$gene
write.csv(result,file="result")
