# This scipt run DESeq2 to run differential expression analysis of RNAseq count data.
# It needs design metrix (design_metadata.txt).
# 

require("DESeq2")

setwd("/Users/rkumar/mac_icr/B97_jenni_rnaseq")

#
# Create a table of read counts for b97_data
#
files <- c(
	"Sample_B97_0001.htseq.txt",
	"Sample_B97_0002.htseq.txt",
	"Sample_B97_0003.htseq.txt",
	"Sample_B97_0004.htseq.txt",
	"Sample_B97_0005.htseq.txt",
	"Sample_B97_0006.htseq.txt",
	"Sample_B97_0007.htseq.txt",
	"Sample_B97_0008.htseq.txt",
	"Sample_B97_0009.htseq.txt",
	"Sample_B97_0010.htseq.txt",
	"Sample_B97_0011.htseq.txt",
	"Sample_B97_0012.htseq.txt",
	"Sample_B97_0013.htseq.txt",
	"Sample_B97_0014.htseq.txt",
	"Sample_B97_0015.htseq.txt",
	"Sample_B97_0016.htseq.txt",
	"Sample_B97_0017.htseq.txt",
	"Sample_B97_0018.htseq.txt"
	)

data <- read.table(
	file=files[1],
	header=FALSE,
	sep="\t",
	row.names=1,
	stringsAsFactors=FALSE
	)

for(file in files[-1]){
	temp <- read.table(
		file=file,
		header=FALSE,
		sep="\t",
		row.names=1,
		stringsAsFactors=FALSE
		)
	data <- cbind(
		data,
		temp
		)
}
colnames(data) <- files

write.table(
	data,
	file="b97_data/combined_b97_rnaseq_files.txt",
	row.names=TRUE,
	col.names=TRUE,
	sep="\t"
	)


#
# read in the table of read counts for b97_data
#
data <- read.table(
	file="b97_data/combined_b97_rnaseq_files.txt",
	row.names=1,
	header=TRUE,
	sep="\t"
	)


# ============================================= #


design <- read.table(
	file="/Users/rkumar/mac_icr/B97_jenni_rnaseq/design_metadata.txt",
	header=TRUE,
	sep="\t",
	row.names=1
	)

# design cols
# cellline
# mutant
# treatment



#
# ======================================================= #
#



#
# A3B par Vs A3B DMSO
#
this_analysis <- "A3B_par_vs_A3B_DMSO"
this_data <- data[,c(1,2,3,7,8,9)]
this_design <- design[c(1,2,3,7,8,9),]
	dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ treatment
		)
	dds <- DESeq(dds)
	res <- results(dds)
	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
	plotDispEsts( dds )
	plotMA(res)
	dev.off()
	write.table(
		res,
		file=paste(
			this_analysis,
			"_full_results.txt",
			sep=""
			),
		row.names=TRUE,
		col.names=TRUE,
		sep="\t"
		)



#
# A3B_TP53_par_vs_A3B_TP53_DMSO
#
this_analysis <- "A3B_TP53_par_vs_A3B_TP53_DMSO"
this_data <- data[,c(4,5,6,13,14,15)]
this_design <- design[c(4,5,6,13,14,15),]
	dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ treatment
		)
	dds <- DESeq(dds)
	res <- results(dds)
	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
	plotDispEsts( dds )
	plotMA(res)
	dev.off()
	write.table(
		res,
		file=paste(
			this_analysis,
			"_full_results.txt",
			sep=""
			),
		row.names=TRUE,
		col.names=TRUE,
		sep="\t"
		)


#
# A3B_DMSO_vs_A3B_DOX
#
this_analysis <- "A3B_DMSO_vs_A3B_DOX"
this_data <- data[,c(7,8,9,10,11,12)]
this_design <- design[c(7,8,9,10,11,12),]
	dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ treatment
		)
	dds <- DESeq(dds)
	res <- results(dds)
	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
	plotDispEsts( dds )
	plotMA(res)
	dev.off()
	write.table(
		res,
		file=paste(
			this_analysis,
			"_full_results.txt",
			sep=""
			),
		row.names=TRUE,
		col.names=TRUE,
		sep="\t"
		)

#
# A3B_TP53_DMSO_vs_A3B_TP53_DOX
#
this_analysis <- "A3B_TP53_DMSO_vs_A3B_TP53_DOX"
this_data <- data[,c(13,14,15,16,17,18)]
this_design <- design[c(13,14,15,16,17,18),]
	dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ treatment
		)
	dds <- DESeq(dds)
	res <- results(dds)
	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
	plotDispEsts( dds )
	plotMA(res)
	dev.off()
	write.table(
		res,
		file=paste(
			this_analysis,
			"_full_results.txt",
			sep=""
			),
		row.names=TRUE,
		col.names=TRUE,
		sep="\t"
		)


#
# A3B_DMSO_vs_A3B_TP53_DMSO
#
this_analysis <- "A3B_DMSO_vs_A3B_TP53_DMSO"
this_data <- data[,c(7,8,9,13,14,15)]
this_design <- design[c(7,8,9,13,14,15),]
	dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ mutant
		)
	dds <- DESeq(dds)
	res <- results(dds)
	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
	plotDispEsts( dds )
	plotMA(res)
	dev.off()
	write.table(
		res,
		file=paste(
			this_analysis,
			"_full_results.txt",
			sep=""
			),
		row.names=TRUE,
		col.names=TRUE,
		sep="\t"
		)


#
# A3B_DOX_vs_A3B_TP53_DOX
#
this_analysis <- "A3B_DOX_vs_A3B_TP53_DOX"
this_data <- data[,c(10,11,12,16,17,18)]
this_design <- design[c(10,11,12,16,17,18),]
	dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ mutant
		)
	dds <- DESeq(dds)
	res <- results(dds)
	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
	plotDispEsts( dds )
	plotMA(res)
	dev.off()
	write.table(
		res,
		file=paste(
			this_analysis,
			"_full_results.txt",
			sep=""
			),
		row.names=TRUE,
		col.names=TRUE,
		sep="\t"
		)


#
# A3B_par_vs_A3B_DOX
#
this_analysis <- "A3B_par_vs_A3B_DOX"
this_data <- data[,c(1,2,3,10,11,12)]
this_design <- design[c(1,2,3,10,11,12),]
        dds = DESeqDataSetFromMatrix(
                countData=this_data,
                colData=this_design,
                design= ~ mutant
                )
        dds <- DESeq(dds)
        res <- results(dds)
        pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
        hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
        plotDispEsts( dds )
        plotMA(res)
        dev.off()
        write.table(
                res,
                file=paste(
                        this_analysis,
                        "_full_results.txt",
                        sep=""
                        ),
                row.names=TRUE,
                col.names=TRUE,
                sep="\t"
                )

#
# A3B_TP53_par_vs_A3B_TP53_DOX
#
this_analysis <- "A3B_TP53_par_vs_A3B_TP53_DOX"
this_data <- data[,c(4,5,6,16,17,18)]
this_design <- design[c(4,5,6,16,17,18),]
        dds = DESeqDataSetFromMatrix(
                countData=this_data,
                colData=this_design,
                design= ~ mutant
                )
        dds <- DESeq(dds)
        res <- results(dds)
        pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
        hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
        plotDispEsts( dds )
        plotMA(res)
        dev.off()
        write.table(
                res,
                file=paste(
                        this_analysis,
                        "_full_results.txt",
                        sep=""
                        ),
                row.names=TRUE,
                col.names=TRUE,
                sep="\t"
                )


#
# A3B_par_vs_A3B_TP53_par
#
this_analysis <- "A3B_par_vs_A3B_TP53_par"
this_data <- data[,c(1,2,3,4,5,6)]
this_design <- design[c(1,2,3,4,5,6),]
        dds = DESeqDataSetFromMatrix(
                countData=this_data,
                colData=this_design,
                design= ~ mutant
                )
        dds <- DESeq(dds)
        res <- results(dds)
        pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
        hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
        plotDispEsts( dds )
        plotMA(res)
        dev.off()
        write.table(
                res,
                file=paste(
                        this_analysis,
                        "_full_results.txt",
                        sep=""
                        ),
                row.names=TRUE,
                col.names=TRUE,
                sep="\t"
                )



#
# A3B_DMSO_vs_A3B_TP53_DOX
#
this_analysis <- "A3B_DMSO_vs_A3B_TP53_DOX"
this_data <- data[,c(7,8,9,16,17,18)]
this_design <- design[c(7,8,9,16,17,18),]
        dds = DESeqDataSetFromMatrix(
                countData=this_data,
                colData=this_design,
                design= ~ mutant
                )
        dds <- DESeq(dds)
        res <- results(dds)
        pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
        hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
        plotDispEsts( dds )
        plotMA(res)
        dev.off()
        write.table(
                res,
                file=paste(
                        this_analysis,
                        "_full_results.txt",
                        sep=""
                        ),
                row.names=TRUE,
                col.names=TRUE,
                sep="\t"
                )
