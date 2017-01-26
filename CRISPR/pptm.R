
setwd("./") 
source("./pptm_functions.R")


sample <- read.table(
	file="all.count.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)



sample_pptm <- get_pptm(
	sample
	)

write.table(
	sample_pptm,
	file="all.count.pptm.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)
