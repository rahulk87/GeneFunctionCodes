## This script align the raw fastq reads to reference genome usin BWA and call variants using GATK
## You may need to change the file path, software module version etc.
# NOTE: Some commands can change if you are using the different version of software.

#!/bin/bash
#BSUB -u rahul.kumar\@icr.ac.uk
#BSUB -J BTBC14_P
#BSUB -e /scratch/DBC/GENFUNC/NGS_Projects/BTBC14_Steve_WES/BTBC14_P/gatk/run.err
#BSUB -o /scratch/DBC/GENFUNC/NGS_Projects/BTBC14_Steve_WES/BTBC14_P/gatk/run.out
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -P DBCDOBZAK
#BSUB -W 96:00
#BSUB -q normal
#

# load modules for bwa, picard etc
module load bwa/0.7.9a
module load samtools/1.2
module load java/1.6.0u45
module load picard/1.130


# Required paths
path1="/skaro/DBC/GENFUNC/NGS_projects/BTBC14_Steve_WES/BTBC14_P" # raw data 
path2="/scratch/DBC/GENFUNC/NGS_Projects/BTBC14_Steve_WES/BTBC14_P/gatk" # wroking directory
path3="/scratch/DBC/GENFUNC/NGS_Projects/genomes" # genomes

# Copy and unzip the fastq files
cp $path1/*_1.fq.gz $path2/BTBC14_P.R1.fq.gz
gunzip $path2/BTBC14_P.R1.fq.gz

cp $path1/*_2.fq.gz $path2/BTBC14_P.R2.fq.gz
gunzip $path2/BTBC14_P.R2.fq.gz

# run bwa mem
/apps/bwa/0.7.9a/bwa mem \
-t 8 \
$path3/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
$path2/BTBC14_P.R1.fq \
$path2/BTBC14_P.R2.fq \
> $path2/BTBC14_P.aln.sam

# convert sams to bams
/apps/samtools/1.2/bin/samtools view -bS \
$path2/BTBC14_P.aln.sam \
> $path2/BTBC14_P.aln.bam

# remove the fasqt files and sam file to save space
rm $path2/BTBC14_P.R1.fq 
rm $path2/BTBC14_P.R2.fq 
rm $path2/BTBC14_P.aln.sam

# sort and index bams
/apps/samtools/1.2/bin/samtools sort -m 16000000000 \
$path2/BTBC14_P.aln.bam \
$path2/BTBC14_P.aln.sorted

/apps/samtools/1.2/bin/samtools index \
$path2/BTBC14_P.aln.sorted.bam

# fix read group information
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar AddOrReplaceReadGroups \
I=$path2/BTBC14_P.aln.sorted.bam \
O=$path2/BTBC14_P.aln.sorted.RGfixed.bam \
RGID=BTBC14_P \
RGLB=BTBC14_P \
RGPL=illumina \
RGPU=NNNNNNNN \
RGSM=BTBC14_P \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
TMP_DIR=./ \
VALIDATION_STRINGENCY=LENIENT

# remove duplicates
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=true \
INPUT=$path2/BTBC14_P.aln.sorted.RGfixed.bam \
OUTPUT=$path2/BTBC14_P.aln.sorted.RGfixed.rmdup.bam \
METRICS_FILE=$path2/BTBC14_P.metrics \
ASSUME_SORTED=true \
CREATE_INDEX=true \
TMP_DIR=./ \
VALIDATION_STRINGENCY=LENIENT


# get WGS metrics
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar CollectWgsMetrics \
INPUT=$path2/BTBC14_P.aln.sorted.RGfixed.rmdup.bam \
OUTPUT=$path2/BTBC14_P.aln.sorted.RGfixed.rmdup.bam.WGSmetrics \
REFERENCE_SEQUENCE=$path3/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
TMP_DIR=./ \
MAX_RECORDS_IN_RAM=5000000 \
VALIDATION_STRINGENCY=LENIENT \
INCLUDE_BQ_HISTOGRAM=true

## Reorder bam file
java -Xmx64G -jar /apps/picard-tools/1.130/picard.jar \
ReorderSam \
I=$path2/BTBC14_P.aln.sorted.RGfixed.rmdup.bam \
O=$path2/BTBC14_P.aln.sorted.RGfixed.rmdup_reordered.bam \
R=$path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa 

/apps/samtools/1.2/bin/samtools index \
$path2/BTBC14_P.aln.sorted.RGfixed.rmdup_reordered.bam \

## run realigner target creator
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC14_P.aln.sorted.RGfixed.rmdup_reordered.bam \
-known $path3/1000G_indels_for_realignment.b37.vcf \
-o $path2/BTBC14_P.forIndelRealigner.intervals \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC \
-filterMBQ  \
-filterNoBases

## run Indel aligner
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC14_P.aln.sorted.RGfixed.rmdup_reordered.bam \
-targetIntervals $path2/BTBC14_P.forIndelRealigner.intervals \
-known $path3/1000G_indels_for_realignment.b37.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path2/BTBC14_P.forIndelRealigner.bam \

## run picard Sortsam
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/picard-tools/1.130/picard.jar SortSam \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM= 10000000 \
SORT_ORDER=coordinate \
I=$path2/BTBC14_P.forIndelRealigner.bam \
O=$path2/BTBC14_P.forIndelRealigner.sorted.bam \
CREATE_INDEX=true \
TMP_DIR= $path2

## run BaseRecalibrator
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC14_P.forIndelRealigner.sorted.bam \
-o $path2/BTBC14_P.forIndelRealigner.sorted.bam.table \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC -filterMBQ  -filterNoBases \
-knownSites $path3/dbsnp_132.b37.vcf \
-knownSites $path3/1000G_indels_for_realignment.b37.vcf

## run PrintReads
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T PrintReads \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC14_P.forIndelRealigner.sorted.bam \
-o $path2/BTBC14_P.recaliberated.bam \
-BQSR $path2/BTBC14_P.forIndelRealigner.sorted.bam.table \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC -filterMBQ  -filterNoBases

## run HaplotypeCaller
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
--dbsnp $path3/dbsnp_132.b37.vcf \
-I $path2/BTBC14_P.recaliberated.bam \
-o $path2/BTBC14_P.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-nt 1 \
-nct 1 \
-dt NONE \
-A Coverage \
-filterRNC -filterMBQ  -filterNoBases \
-stand_call_conf 30.0 \
-stand_emit_conf 10.0

## run select SNPs
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T SelectVariants \
-selectType SNP \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-V $path2/BTBC14_P.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path2/BTBC14_P_SNPs.vcf

## run VariantRecalibrator
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-input $path2/BTBC14_P_SNPs.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $path3/hapmap_3.3.b37.vcf \
-resource:omni,known=false,training=true,truth=false,prior=12.0 $path3/1000G_omni2.5.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $path3/dbsnp_132.b37.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 $path3/1000G_phase3_v4_20130502.sites.vcf \
-an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP \
-recalFile BTBC14_P.recal \
-tranchesFile BTBC14_P.tranches \
-rscriptFile BTBC14_P.plots.R

# Making model plots
perl -pi -e 's/,legend=FALSE//g' BTBC14_P.plots.R
R CMD BATCH BTBC14_P.plots.R

## run ApplyRecalibration
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-input $path2/BTBC14_P_SNPs.vcf \
--ts_filter_level 99.0 \
-tranchesFile BTBC14_P.tranches \
-recalFile BTBC14_P.recal \
-mode SNP \
-o $path2/BTBC14_P_SNPs_VQSR.vcf

## run select INDELs
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T SelectVariants \
-selectType INDEL \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-V $path2/BTBC14_P.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path2/BTBC14_P_INDELs.vcf

## run VariantRecalibrator for INDELs
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-input $path2/BTBC14_P_INDELs.vcf \
-resource:mills,known=false,training=true,truth=true,prior=12.0 $path3/Mills_and_1000G_gold_standard.indels.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $path3/dbsnp_132.b37.vcf \
--maxGaussians 4 \
-an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode INDEL \
-recalFile BTBC14_P_INDELs.recal \
-tranchesFile BTBC14_P_INDELs.tranches \
-rscriptFile BTBC14_P_INDELs.plots.R

## run ApplyRecalibration for INDELs
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-input $path2/BTBC14_P_INDELs.vcf \
--ts_filter_level 99.0 \
-tranchesFile BTBC14_P_INDELs.tranches \
-recalFile BTBC14_P_INDELs.recal \
-mode INDEL \
-o $path2/BTBC14_P_INDELs_VQSR.vcf

## run DepthOfCoverage
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC14_P.recaliberated.bam \
-o $path2/BTBC14_P.DepthOfCoverage \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-ct 4 -ct 6 -ct 10 \
-pt readgroup
