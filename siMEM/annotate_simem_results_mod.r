#
# Add zGARP values and group medians
# to siMEMv10 output
#

setwd("./")

i=1
list=read.csv("list1",sep="\t",header=FALSE, stringsAsFactors=FALSE)
        for (genes in list$V1)
        {
	i=i+1
        x <- genes
        y <- tolower(x)
        print (x)
	#### creating mutation file for each gene ####
        file<-read.table("../combined_exome_cnv_func_muts_131216_gt3_mod.txt",sep="\t",header=TRUE)
        cells <- as.data.frame(file[,1])
        gene_mut<-cbind(cells,file[,i])
        names(gene_mut) <- c("cell_line", "calls")
        write.table(gene_mut,file="gene_mut.txt",sep="\t",row.names=F,quote = FALSE)
        ##############################################

simem_results_file <- paste("result_",y,sep="")
cell_line_grouping_file <- "gene_mut.txt"
zgarp_file <- "../breast_zgarp_intercell_names_transposed.txt"
simem_p_threshold <- 0.05
diff <- 0

results <- read.table(
	file=simem_results_file,
	row.names=1,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)
results <- results[,1:5]
rownames(results) <- results$symbol


groups <- read.table(
	file=cell_line_grouping_file,
	header=TRUE,
	row.names=1,
	sep="\t"
	)


zgarp <- read.table(
	file=zgarp_file,
	header=TRUE,
	row.names=1,
	sep="\t"
	)


get_median <- function(x){
	return(median(x, na.rm=TRUE))
}


# define the RB1 defective and wild type groups
rb1_defective <- rownames(groups)[which(groups$calls == 1)]
l_rb1_defective <- length(rb1_defective)

rb1_wt <- rownames(groups)[which(groups$calls == 0)]
l_rb1_wt <- length(rb1_wt)


# calculate the group median zGARPs

rb1_defective_zgarp_meds <- apply(
	zgarp[rb1_defective,],
	2,
	get_median
	)

rb1_wt_zgarp_meds <- apply(
	zgarp[rb1_wt,],
	2,
	get_median
	)


# select significant siMEM results

results_sig <- results[which(
	results$pvalue_high <= simem_p_threshold & results$difference_high < diff
	),]

# add median and individual zGARP scores

#results_sig_with_meds <- cbind(
#	results_sig,
#	RB1_defective_median_zGARP=rb1_defective_zgarp_meds[rownames(results_sig)],
#	RB1_wild_type_median_zGARP=rb1_wt_zgarp_meds[rownames(results_sig)],
#	t(zgarp[rb1_defective,gsub("-", ".", rownames(results_sig))]),
#	t(zgarp[rb1_wt,gsub("-", ".", rownames(results_sig))])
#	)

#write.table(
#        results_sig_with_meds,
#        file=paste("result_sig_",y,".txt")
#        sep="\t",
#        col.names=TRUE,
#        row.names=FALSE,
#        quote=FALSE
#        )


#:1

# write the full set
#

results_with_meds <- cbind(
	results_sig,
	y,
	RB1_defective_median_zGARP=rb1_defective_zgarp_meds[rownames(results_sig)],
	RB1_wild_type_median_zGARP=rb1_wt_zgarp_meds[rownames(results_sig)],
	l_rb1_defective,
	l_rb1_wt
	)


write.table(
	results_with_meds,
	file=paste("result_medians_sig_",y,".txt",sep=""),
	sep="\t",
	col.names=TRUE,
	row.names=FALSE,
	quote=FALSE
	)
file.remove("gene_mut.txt")
}
