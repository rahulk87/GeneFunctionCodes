## This script calculate copy number using PureCN package.
## Input files: Mutect1 vcf file and GATK coverage file.

## ----load-purecn, echo=FALSE, message=FALSE------------------------------
library(PureCN)

## ----example_files, message=FALSE, warning=FALSE, results='hide'---------
gatk.normal.file <- "BTBC183_germline_cov.txt"
#gatk.tumor.file <- "BTBC183_pdx_cov.txt" 
gatk.tumor.file <- "BTBC183_primary_cov.txt" 
#vcf.file <- "BTBC183.germline.pdx.mutect.vcf"
vcf.file <- "BTBC183.germline.primary.mutect.vcf"

## ----gccorrect-----------------------------------------------------------
gc.gene.file <- "/Users/rkumar/mac_icr/work/Becky_PDXs/BXset/exome/common_files/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod_GC_gene.bed" 
#coverage <- correctCoverageBias(gatk.normal.file, gc.gene.file, "BTBC183_germline_cov_gc.txt")
#coverage1 <- correctCoverageBias(gatk.tumor.file, gc.gene.file, "BTBC183_pdx_cov_gc.txt")
coverage1 <- correctCoverageBias(gatk.tumor.file, gc.gene.file, "BTBC183_primary_cov_gc.txt")

## ----Reloading GC corrected files files----------------------------------
gatk.normal.file.gc <- "BTBC183_germline_cov_gc.txt"
#gatk.tumor.file.gc <- "BTBC183_pdx_cov_gc.txt"
gatk.tumor.file.gc <- "BTBC183_primary_cov_gc.txt"

## ----runpurecn-----------------------------------------------------------
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file.gc, gatk.tumor.file=gatk.tumor.file.gc, 
    vcf.file=vcf.file, sampleid='BTBC183.germline.primary',
    gc.gene.file=gc.gene.file, 
    #args.filterVcf=list(snp.blacklist=snp.blacklist, stats.file=mutect.stats.file), 
    #args.segmentation=list(exon.weight.file=exon.weight.file), 
    post.optimize=FALSE, plot.cnv=FALSE)

## ----createoutput--------------------------------------------------------
file.rds <- 'BTBC183.germline.primary.rds'
saveRDS(ret, file=file.rds)
pdf('BTBC183.germline.primary.pdf', width=10, height=12)
plotAbs(ret, type='all')
