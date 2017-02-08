
E7<-read.csv("htseq_COV_E7.txt",sep="\t",header=FALSE,row.names=1)
A4<-read.csv("htseq_COV_A4.txt",sep="\t",header=FALSE,row.names=1)
D11<-read.csv("htseq_COV_D11.txt",sep="\t",header=FALSE,row.names=1)
Cas9<-read.csv("htseq_SUM149_CAS9.txt",sep="\t",header=FALSE,row.names=1)
TR1<-read.csv("htseq_SUM149_TR1.txt",sep="\t",header=FALSE,row.names=1)
TR2<-read.csv("htseq_SUM149_TR2.txt",sep="\t",header=FALSE,row.names=1)

require("DESeq2")
setwd("./")

#
# Create a table of read counts for b97_data
#

files <- c(

	"htseq_COV_E7.txt",

	"htseq_COV_A4.txt",

	"htseq_COV_D11.txt",

	"htseq_SUM149_CAS9.txt",

	"htseq_SUM149_TR1.txt",

	"htseq_SUM149_TR2.txt"

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

data_red <- data[rowSums(data > 5),]

write.table(

	data,

	file="combined_rnaseq_files.txt",

	row.names=TRUE,

	col.names=TRUE,

	sep="\t"

	)

write.table(

	data_red,

	file="combined_rnaseq_files_reduced.txt",

	row.names=TRUE,

	col.names=TRUE,

	sep="\t"

	)



#

# read in the table of read counts for b97_data

#

data <- read.table(

	file="combined_rnaseq_files.txt",

	row.names=1,

	header=TRUE,

	sep="\t"

	)





# ============================================= #





design <- read.table(

	file="design_metadata.txt",

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

# E7 vs A4

#

this_analysis <- "E7_A4"

this_data <- data[,c(1,2)]

this_design <- design[c(1,2),]

	dds = DESeqDataSetFromMatrix(

		countData=this_data,

		colData=this_design,

		design= ~ treatment

		)
	dds <- dds[rowSums(counts(dds)) > 5,]	

	dds <- DESeq(dds)

	res <- results(dds,contrast=c("treatment","A4","E7"))

	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)

	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")

	plotDispEsts( dds )

	plotMA(res)

	dev.off()
	
	res1<-as.matrix(res)

	A4_E7<-cbind(E7[rownames(res1),],A4[rownames(res1),],res1)
	
	write.table(

		A4_E7,

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

# E7 and D11

#

this_analysis <- "E7_D11"

this_data <- data[,c(1,3)]

this_design <- design[c(1,3),]

	dds = DESeqDataSetFromMatrix(

		countData=this_data,

		colData=this_design,

		design= ~ treatment

		)
	dds <- dds[rowSums(counts(dds)) > 5,]	

	dds <- DESeq(dds)

	res <- results(dds, contrast=c("treatment","D11","E7"))

	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)

	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")

	plotDispEsts( dds )

	plotMA(res)

	dev.off()
	
	res1<-as.matrix(res)

	D11_E7<-cbind(E7[rownames(res1),],D11[rownames(res1),],res1)

	write.table(

		D11_E7,

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

# Cas9 vs TR1

#

this_analysis <- "Cas9_TR1"

this_data <- data[,c(4,5)]

this_design <- design[c(4,5),]

	dds = DESeqDataSetFromMatrix(

		countData=this_data,

		colData=this_design,

		design= ~ treatment

		)
	dds <- dds[rowSums(counts(dds)) > 5,]	

	dds <- DESeq(dds)

	res <- results(dds, contrast=c("treatment","TR1","Cas9"))

	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)

	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")

	plotDispEsts( dds )

	plotMA(res)

	dev.off()
	
	res1<-as.matrix(res)

	TR1_Cas9<-cbind(Cas9[rownames(res1),],TR1[rownames(res1),],res1)

	write.table(

		TR1_Cas9,

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

# Cas9 vs TR2

#

this_analysis <- "Cas9_TR2"

this_data <- data[,c(4,6)]

this_design <- design[c(4,6),]

        dds = DESeqDataSetFromMatrix(

                countData=this_data,

                colData=this_design,

                design= ~ treatment

                )
        dds <- dds[rowSums(counts(dds)) > 5,]        

        dds <- DESeq(dds)

        res <- results(dds, contrast=c("treatment","TR2","Cas9"))

        pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)

        hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")

        plotDispEsts( dds )

        plotMA(res)

        dev.off()
        
        res1<-as.matrix(res)

	    TR2_Cas9<-cbind(Cas9[rownames(res1),],TR2[rownames(res1),],res1)

        write.table(

                TR2_Cas9,

                file=paste(

                        this_analysis,

                        "_full_results.txt",

                        sep=""

                        ),

                row.names=TRUE,

                col.names=TRUE,

                sep="\t"

                )

