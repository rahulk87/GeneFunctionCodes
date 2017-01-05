#!/bin/bash
#BSUB -u rahul.kumar\@icr.ac.uk
#BSUB -J bwa_BTBC183_pdx
#BSUB -e /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/pdx/bwa.err
#BSUB -o /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/pdx/bwa.out
#BSUB -n 8
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
path1="/scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/raw_data/BTBC183/pdx" # raw data 
path2="/scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/pdx" # wroking directory
path3="/scratch/DBC/GENFUNC/NGS_Projects/genomes" # genomes

# zcat the R1 files together and the R2 files together
cp $path1/*_1.fq.gz $path2/BTBC183_pdx.R1.fq.gz
gunzip $path2/BTBC183_pdx.R1.fq.gz

cp $path1/*_2.fq.gz $path2/BTBC183_pdx.R2.fq.gz
gunzip $path2/BTBC183_pdx.R2.fq.gz

# run bwa mem
/apps/bwa/0.7.9a/bwa mem \
-t 8 \
$path3/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
$path2/BTBC183_pdx.R1.fq \
$path2/BTBC183_pdx.R2.fq \
> $path2/BTBC183_pdx.aln.sam

# convert sams to bams
/apps/samtools/1.2/bin/samtools view -bS \
$path2/BTBC183_pdx.aln.sam \
> $path2/BTBC183_pdx.aln.bam

# remove the fasqt files and sam file to save space
rm $path2/BTBC183_pdx.R1.fq 
rm $path2/BTBC183_pdx.R2.fq 
rm $path2/BTBC183_pdx.aln.sam
comment

# sort and index bams
/apps/samtools/1.2/bin/samtools sort -m 16000000000 \
$path2/BTBC183_pdx.aln.bam \
$path2/BTBC183_pdx.aln.sorted

/apps/samtools/1.2/bin/samtools index \
$path2/BTBC183_pdx.aln.sorted.bam

# fix read group information
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar AddOrReplaceReadGroups \
I=$path2/BTBC183_pdx.aln.sorted.bam \
O=$path2/BTBC183_pdx.aln.sorted.RGfixed.bam \
RGID=BTBC183_pdx \
RGLB=BTBC183_pdx \
RGPL=illumina \
RGPU=NNNNNNNN \
RGSM=BTBC183_pdx \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
TMP_DIR=./ \
VALIDATION_STRINGENCY=LENIENT

# remove duplicates
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=true \
INPUT=$path2/BTBC183_pdx.aln.sorted.RGfixed.bam \
OUTPUT=$path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup.bam \
METRICS_FILE=$path2/BTBC183_pdx.metrics \
ASSUME_SORTED=true \
CREATE_INDEX=true \
TMP_DIR=./ \
VALIDATION_STRINGENCY=LENIENT


# get WGS metrics
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar CollectWgsMetrics \
INPUT=$path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup.bam \
OUTPUT=$path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup.bam.WGSmetrics \
REFERENCE_SEQUENCE=$path3/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
TMP_DIR=./ \
MAX_RECORDS_IN_RAM=5000000 \
VALIDATION_STRINGENCY=LENIENT \
INCLUDE_BQ_HISTOGRAM=true

## Reorder bam file
java -Xmx64G -jar /apps/picard-tools/1.130/picard.jar \
ReorderSam \
I=$path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup.bam \
O=$path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup_reordered.bam \
R=$path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa 

/apps/samtools/1.2/bin/samtools index \
$path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup_reordered.bam \

## run realigner target creator
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup_reordered.bam \
-known $path3/1000G_indels_for_realignment.b37.vcf \
-o $path2/BTBC183_pdx.forIndelRealigner.intervals \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC \
-filterMBQ  \
-filterNoBases

## run Indel aligner
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC183_pdx.aln.sorted.RGfixed.rmdup_reordered.bam \
-targetIntervals $path2/BTBC183_pdx.forIndelRealigner.intervals \
-known $path3/1000G_indels_for_realignment.b37.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path2/BTBC183_pdx.forIndelRealigner.bam \

## run picard Sortsam
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/picard-tools/1.130/picard.jar SortSam \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM= 10000000 \
SORT_ORDER=coordinate \
I=$path2/BTBC183_pdx.forIndelRealigner.bam \
O=$path2/BTBC183_pdx.forIndelRealigner.sorted.bam \
CREATE_INDEX=true \
TMP_DIR= $path2

## run BaseRecalibrator
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC183_pdx.forIndelRealigner.sorted.bam \
-o $path2/BTBC183_pdx.forIndelRealigner.sorted.bam.table \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC -filterMBQ  -filterNoBases \
-knownSites $path3/dbsnp_132.b37.vcf \
-knownSites $path3/1000G_indels_for_realignment.b37.vcf

## run PrintReads
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T PrintReads \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC183_pdx.forIndelRealigner.sorted.bam \
-o $path2/BTBC183_pdx.recaliberated.bam \
-BQSR $path2/BTBC183_pdx.forIndelRealigner.sorted.bam.table \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC -filterMBQ  -filterNoBases

## run DepthOfCoverage
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path2/BTBC183_pdx.recaliberated.bam \
-o $path2/BTBC183_pdx.DepthOfCoverage \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-ct 4 -ct 6 -ct 10 \
-pt readgroup
