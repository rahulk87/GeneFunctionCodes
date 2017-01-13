setwd("/Users/rkumar/mac_icr/work/crispr/D48/shalign") 
source("CRISPR_screen_functions_110704.R")

# Sample TO
sampleT0 <- read.table(
	file="T0.results",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)


# Sample T1 (Drug treated)
sampleT1 <- read.table(
	file="T1.results",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

# Get rid of controls (Olfr and PLK1)
# as makes no sense to consider with 
# samples. We find that all controls
# with a common sequence collapse to
# a single id (regardless of the plate
# they were on)

sampleT0_no_controls <- remove_controls(sampleT0)
sampleT1_no_controls <- remove_controls(sampleT1)

#
# find low abundance guides (unreliable measurements
# that cause chos in the analysis)
#
	
low_abundance_guides <- get_low_abundance_guides(
	x=sampleT0_no_controls,
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
# Positive selection screens
#
# ========================== #

sampleT0_pptm <- get_pptm(
	sampleT0_no_controls
	)
sampleT1_pptm <- get_pptm(
	sampleT1_no_controls
	)

# Z-scores after removing low count

# viabilty effect
znorm_T0_T1 <- znorm(
	x0=sampleT0_pptm,
	x1=sampleT1_pptm,
	exclude_guides=low_abundance_guides
	)

#
# Histograms of Z-scores (as example)
# for two week drug and viability effect
#
pdf("Histograms_of_Z_Scores_excluding_common_low_abundance.pdf", width=5, height=5)
#viabilty effect at 1wks
hist(
	znorm_T0_T1$x1_x0_zscore,
	xlab="Z-score",
	main="Guides > 50 counts"
	)
plot(
        sampleT0$total.hits,
        sampleT1$total.hits,
        xlab="T0 (Total count)",
        ylab="T1 (Total count)",
        log="xy"
        )
plot(
        sampleT0_pptm$pptm_psuedo,
        sampleT1_pptm$pptm_psuedo,
        xlab="T0 (PPTM count)",
        ylab="T1 (PPTM count)",
	log="xy"
        )
plot(
	znorm_T0_T1$x0_log,
	znorm_T0_T1$x1_log,
	xlab="T0 (log PPTM count)",
	ylab="T1 (log PPTM count)"
	)

dev.off()

#
# join up Z-scores with original data
# and write out
#

#######n comment out below

# viability at 1 wk
write.table(
	znorm_T0_T1,
	file="Sample_T0_T1_Zscores.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)

