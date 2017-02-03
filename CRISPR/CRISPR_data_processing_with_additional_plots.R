
# Z-normalisation 26th Jan 2017.

setwd("./") 
source("./CRISPR_screen_functions_110704.R")


# Sample 9 is HMLE T=O
T0 <- read.table(
	file="../../../shalign/D56/Sample_D56_0007/merged.fastq.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)


# Sample 10 is HMLE DMSO 1W
DMSO_1 <- read.table(
	file="../../../shalign/D56/Sample_D56_0001/merged.fastq.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

# Sample 11 is HMLE DMSO 2W
DMSO_2 <- read.table(
	file="../../../shalign/D56/Sample_D56_0004/merged.fastq.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

# Sample 12 is HMLE AZD 1W
Drug_1 <- read.table(
	file="../../../shalign/D56/Sample_D56_0002/merged.fastq.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)
	
# Sample 13 is HMLE AZD 2W
Drug_2 <- read.table(
	file="../../../shalign/D56/Sample_D56_0005/merged.fastq.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

## Reads count histogram
pdf("Histogram_read_counts_1000_all_guides.pdf")
#par(mfrow=c(3,2))
hist(T0$total.hits,xlim=c(0,10000),ylim=c(0,30000),xlab="Read counts",main="T0")
hist(DMSO_1$total.hits,xlim=c(0,12000),ylim=c(0,25000),xlab="Read counts",main="DMSO (1 Week)")
hist(DMSO_2$total.hits,xlim=c(0,10000),ylim=c(0,30000),xlab="Read counts",main="DMSO (2 Weeks)")
hist(Drug_1$total.hits,xlim=c(0,12000),ylim=c(0,30000),xlab="Read counts",main="Drug (1 Week)")
hist(Drug_2$total.hits,xlim=c(0,12000),ylim=c(0,25000),xlab="Read counts",main="Drug (2 Weeks)")
dev.off()

# Get rid of controls (Olfr and PLK1)
# as makes no sense to consider with 
# samples. We find that all controls
# with a common sequence collapse to
# a single id (regardless of the plate
# they were on)

T0_no_controls <- remove_controls(T0)
DMSO_1_no_controls <- remove_controls(DMSO_1)
DMSO_2_no_controls <- remove_controls(DMSO_2)
Drug_1_no_controls <- remove_controls(Drug_1)
Drug_2_no_controls <- remove_controls(Drug_2)

#
# find low abundance guides (unreliable measurements
# that cause chos in the analysis)
#
	
#HMLE3 screens	
low_abundance_guides <- get_low_abundance_guides(
	x=T0_no_controls,
	threshold=50, #default value
	)

#strip info after gene name so we can count genes
low_abundance_genes <- gsub(
	"_.+",
	"",
	low_abundance_guides,
	perl=TRUE
	)


# ========================== #
#
# Negative selection screens
#
# ========================== #



# samples are negative selection screens.
# 9 is t0,
# 10,11 are DMSO T=1 and T=2
# 12,13 are T=1 and T=2 plus AZD
# 14,15 are T=1 and T=2 plus 877 0.1
# 16,17 are T=1 and T=2 plus 877 0.5
# 18,19 are T=1 and T=2 plus 4070 0.1
# 20,21 are T=1 and T=2 plus 4070 0.5


T0_pptm <- get_pptm(
	T0_no_controls
	)
DMSO_1_pptm <- get_pptm(
	DMSO_1_no_controls
	)
DMSO_2_pptm <- get_pptm(
	DMSO_2_no_controls
	)
Drug_1_pptm <- get_pptm(
	Drug_1_no_controls
	)
Drug_2_pptm <- get_pptm(
	Drug_2_no_controls
	)

# Z-scores after removing low count
# rows the t0/no-drug sample


# viabilty effect at 1wk
znorm_10_9 <- znorm(
	x0=T0_pptm,
	x1=DMSO_1_pptm,
	exclude_guides=low_abundance_guides
	)

# viabilty effect at 2wks
znorm_11_9 <- znorm(
	x0=T0_pptm,
	x1=DMSO_2_pptm,
	exclude_guides=low_abundance_guides
	)


# AZD drug effect at 1wk
znorm_12_10 <- znorm(
	x0=DMSO_1_pptm,
	x1=Drug_1_pptm,
	exclude_guides=low_abundance_guides
	)

genes_names <- gsub(
        "_.+",
        "",
        znorm_12_10$shrna.id,
        perl=TRUE
        )

znorm_12_10_genes <- cbind(genes_names, znorm_12_10)
write.table(znorm_12_10_genes,sep="\t",quote=FALSE,file="x.txt")

##

#
# Histograms of Z-scores (as example)
# for two week drug and viability effect
#
pdf("Histograms_of_DE_Scores_excluding_common_low_abundance_2wk.pdf", width=5, height=5)


# viabilty effect at 2wks
hist(
	znorm_11_9$x1_x0_zscore,
	xlab="VE Z-score (2 weeks)",
	main="Guides > 50 counts"
	)

# AZD drug effect at 2wks
hist(
	znorm_13_11$x1_x0_zscore,
	xlab="DE AZD Z-score (2 weeks)",
	main="Guides > 50 counts"
	)

dev.off()



# Histograms of Z-scores (as example)
# for one week drug and viability effect
#
pdf("Histograms_of_DE_Scores_excluding_common_low_abundance_1wk.pdf", width=5, height=5)
#viabilty effect at 1wks
hist(
	znorm_10_9$x1_x0_zscore,
	xlab="VE Z-score (1 weeks)",
	main="Guides > 50 counts"
	)

# AZD drug effect at 1wks
hist(
	znorm_12_10$x1_x0_zscore,
	xlab="DE AZD Z-score (1 weeks)",
	main="Guides > 50 counts"
	)
	
dev.off()

#
# join up Z-scores with original data
# and write out
#


#######n comment out below

# viability at 1 wk
write.table(
	znorm_10_9,
	file="Sample10_9_filtered_viability_effect_Zscores_1W_260117.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# viabilty effect at 2wks
write.table(
	znorm_11_9,
	file="Sample11_9_filtered_viability_effect_Zscores_2W_260117.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# AZD drug effect at 1wk
write.table(
	znorm_12_10,
	file="Sample12_10_filtered_AZD_drug_effect_Zscores_1W_260117.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# AZD drug effect at 2wks
write.table(
	znorm_13_11,
	file="Sample13_11_filtered_AZD_drug_effect_Zscores_2W_260117.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)

############## stop commenting



#
# DE correction using VE
#

# vx

vx_1wk_de_corrected <- correct_viability(
	de_zscore_file="Sample12_10_filtered_AZD_drug_effect_Zscores_1W_260117.txt",
	ve_zscore_file="Sample10_9_filtered_viability_effect_Zscores_1W_260117.txt",
	plot_file="Sample12_10_filtered_AZD_drug_effect_Zscores_1W_260117_VE_correction_plots.pdf"
	)

write.table(
	vx_1wk_de_corrected,
	file="Sample12_10_filtered_AZD_drug_effect_Zscores_1W_260117_VE_correction.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)

vx_2wk_de_corrected <- correct_viability(
	de_zscore_file="Sample13_11_filtered_AZD_drug_effect_Zscores_2W_260117.txt",
	ve_zscore_file="Sample11_9_filtered_viability_effect_Zscores_2W_260117.txt",
	plot_file="Sample13_11_filtered_AZD_drug_effect_Zscores_2W_260117_VE_correction_plots.pdf"
	)

write.table(
	vx_2wk_de_corrected,
	file="Sample13_11_filtered_AZD_drug_effect_Zscores_2W_260117_VE_correction.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)

## RAD9A & RAD9B 1 week
genes_names <- gsub(
 	"_.+",
 	"",
 	vx_1wk_de_corrected$shrna.id,
 	perl=TRUE
 	)

vx_1wk_de_corrected_genes <- cbind(genes_names, vx_1wk_de_corrected)
	
## Plot
png("DE_1weeks.png",height=2800,width=2800,res=400)
with(vx_1wk_de_corrected_genes,plot(x1_log,corrected_de,col="grey",xlab="log2 read counts",ylab="Corrected drug Effect Z score (1 Week)",main="1000 screen (1 Week)"))
abline(h=3,col=4,lty=2)
abline(h=2,col=3,lty=2)
abline(h=-3,col=4,lty=2)
abline(h=-2,col=3,lty=2)
with(subset(vx_1wk_de_corrected_genes, corrected_de > 3), points(x1_log, corrected_de, col="black"))
with(subset(vx_1wk_de_corrected_genes, corrected_de < -3), points(x1_log, corrected_de, col="black"))
with(subset(vx_1wk_de_corrected_genes, genes_names=="RAD9A"), points(x1_log, corrected_de, col="red",pch=16, bg=24))
with(subset(vx_1wk_de_corrected_genes, genes_names=="RAD9B"), points(x1_log, corrected_de, col="green",pch=16, bg=24))
legend("topright",legend=c("RAD9A gRNAs", "RAD9B gRNAs"),cex=1.2,col=c("red","green"),pch=c(16,16),pt.bg=c("red","green"))
legend("bottomright",legend=c("Corrected DE Z cutoff -2/+2", "Corrected DE Z cutoff -3/+3"),cex=1.2,col=c("green","blue"),lty=c(2,2))
dev.off()

## RAD9A & RAD9B 2 weeks
genes_names <- gsub(
 	"_.+",
 	"",
 	vx_2wk_de_corrected$shrna.id,
 	perl=TRUE
 	)

vx_2wk_de_corrected_genes <- cbind(genes_names, vx_2wk_de_corrected)
	
## Plot
png("DE_2weeks.png",height=2800,width=2800,res=400)
with(vx_2wk_de_corrected_genes,plot(x1_log,corrected_de,col="grey",xlab="log2 read counts",ylab="Corrected drug Effect Z score (2 Weeks)",main="1000 screen (2 Weeks)"))
abline(h=3,col=4,lty=2)
abline(h=2,col=3,lty=2)
abline(h=-3,col=4,lty=2)
abline(h=-2,col=3,lty=2)
with(subset(vx_2wk_de_corrected_genes, corrected_de > 3), points(x1_log, corrected_de, col="black"))
with(subset(vx_2wk_de_corrected_genes, corrected_de < -3), points(x1_log, corrected_de, col="black"))
with(subset(vx_2wk_de_corrected_genes, genes_names=="RAD9A"), points(x1_log, corrected_de, col="red",pch=16, bg=24))
with(subset(vx_2wk_de_corrected_genes, genes_names=="RAD9B"), points(x1_log, corrected_de, col="green",pch=16, bg=24))
legend("topright",legend=c("RAD9A gRNAs", "RAD9B gRNAs"),cex=1.2,col=c("red","green"),pch=c(16,16),pt.bg=c("red","green"))
legend("bottomright",legend=c("Corrected DE Z cutoff -2/+2", "Corrected DE Z cutoff -3/+3"),cex=1.2,col=c("green","blue"),lty=c(2,2))
dev.off()


# Read count intensity explains Z-score well. The most extreme Z-scores
# all have low read count (≤100?)
#

# 1week AZD
pdf(
	"Sample12_10_filtered_AZD_drug_effect_Zscores_1W_260117_Intensity_vs_Z-score_plot.pdf",
	width=5,
	height=5
	)
plot(
	vx_1wk_de_corrected$x1_x0_zscore,
	log2(vx_1wk_de_corrected$total.hits),
	xlab="drug effect Z-score",
	ylab="log2 read count (DMSO)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="1 week DMSO",
	xlim=c(-15,10)
	)
plot(
	vx_1wk_de_corrected$x1_x0_zscore,
	log2(vx_1wk_de_corrected$total.hits.1),
	xlab="drug effect Z-score",
	ylab="log2 read count (AZD)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="1 week AZD",
	xlim=c(-15,10)
	)
plot(
	vx_1wk_de_corrected$corrected_de,
	log2(vx_1wk_de_corrected$total.hits),
	xlab="corrected drug effect Z-score",
	ylab="log2 read count (DMSO)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="1 week DMSO",
	xlim=c(-15,10)
	)
plot(
	vx_1wk_de_corrected$corrected_de,
	log2(vx_1wk_de_corrected$total.hits.1),
	xlab="corrected drug effect Z-score",
	ylab="log2 read count (AZD)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="1 week AZD",
	xlim=c(-15,10)
	)
dev.off()

#2 week AZD

pdf(
	"Sample13_11_filtered_AZD_drug_effect_Zscores_2W_260117_Intensity_vs_Z-score_plot.pdf",
	width=5,
	height=5
	)
plot(
	vx_2wk_de_corrected$x1_x0_zscore,
	log2(vx_2wk_de_corrected$total.hits),
	xlab="drug effect Z-score",
	ylab="log2 read count (DMSO)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="2 week DMSO",
	xlim=c(-15,10)
	)
plot(
	vx_2wk_de_corrected$x1_x0_zscore,
	log2(vx_2wk_de_corrected$total.hits.1),
	xlab="drug effect Z-score",
	ylab="log2 read count (AZD)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="2 week AZD",
	xlim=c(-15,10)
	)
plot(
	vx_2wk_de_corrected$corrected_de,
	log2(vx_2wk_de_corrected$total.hits),
	xlab="corrected drug effect Z-score",
	ylab="log2 read count (DMSO)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="2 week DMSO",
	xlim=c(-15,10)
	)
plot(
	vx_2wk_de_corrected$corrected_de,
	log2(vx_2wk_de_corrected$total.hits.1),
	xlab="corrected drug effect Z-score",
	ylab="log2 read count (AZD)",
	pch=19,
	col=rgb(0,0,0,0.25),
	main="2 week AZD",
	xlim=c(-15,10)
	)
dev.off()


