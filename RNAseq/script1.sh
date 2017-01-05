## This script align RNAseq reads to the refernce genome using tophat2 and prepare the assembly using cufflinks
## You may need to change the file path, software module version etc.
 
#!/bin/bash
#BSUB -u rahul.kumar\@icr.ac.uk
#BSUB -J tophat_BTBC183_primary
#BSUB -e /scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/analysis/BTBC183/primary/tophat.err
#BSUB -o /scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/analysis/BTBC183/primary/tophat.out
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -P DBCDOBZAK
#BSUB -W 96:00
#BSUB -q normal

module load boost/1.60.0
module load samtools/0.1.19
module load bowtie2/2.2.6 
module load bamtools/2.4.0
module load cufflinks/2.2.1
module load tophat/2.1.0
module load java/sun7
module load picard-tools/1.141
module load python/2.7.11 intel/comp/17.0.1 intel/mkl/11.3u3

path1="/scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/analysis/BTBC183/primary"
path2="/scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/raw_data/BTBC183/primary"
path3="/scratch/DBC/GENFUNC/rkumar/NGS_projects/genomes"

tophat2 \
-p 12 \
-g 1 \
--no-coverage-search \
--library-type fr-firststrand \
-G $path3/Homo_sapiens.GRCh37.82.chr.gtf \
--tmp-dir $path1/tmp \
-o $path1 \
$path3/ens61_human_chromosomes_and_MT.fa \
$path2/*_1.fq.gz \
$path2/*_2.fq.gz

/apps/java/sun7/1.7.0u80/bin/java -Xmx60g -XX:ParallelGCThreads=1 -jar /apps/picard-tools/1.141/picard.jar ReorderSam \
TMP_DIR=$path1 \
I=$path1/accepted_hits.bam \
O=$path1/accepted_hits.reorder.bam \
REFERENCE=$path3/ens61_human_chromosomes_and_MT.fa \
CREATE_INDEX=true 

/apps/java/sun7/1.7.0u80/bin/java -Xmx60g -XX:ParallelGCThreads=1 -jar /apps/picard-tools/1.141/picard.jar AddOrReplaceReadGroups \
TMP_DIR=$path1 \
I=$path1/accepted_hits.reorder.bam \
O=$path1/accepted_hits.reorder.sort.bam \
RGID=BTBC183_primary \
RGLB=BTBC183_primary \
RGPL=ILLUMINA \
RGPU= BTBC183_primary \
RGSM=BTBC183_primary \
SORT_ORDER=coordinate \
CREATE_INDEX=true

## guided assembly
/apps/cufflinks/2.2.1/cufflinks \
-p 12 \
-G $path3/Homo_sapiens.GRCh37.82.chr.gtf \
--library-type fr-firststrand \
-o $path1/guided \
$path1/accepted_hits.reorder.sort.bam

## de novo assembly
/apps/cufflinks/2.2.1/cufflinks \
-p 12 \
--library-type fr-firststrand \
-o $path1/denovo \
$path1/accepted_hits.reorder.sort.bam

/apps/java/sun7/1.7.0u80/bin/java -Xmx60g -XX:ParallelGCThreads=1 -jar /apps/picard-tools/1.141/picard.jar CollectRnaSeqMetrics \
TMP_DIR=$path1 \
REF_FLAT=$path3/refFlat.txt \
STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
INPUT=$path1/accepted_hits.reorder.sort.bam \
OUTPUT=$path1/accepted_hits.reorder.sort.RNAmetrics \
CHART_OUTPUT=$path1/BTBC183_primary.rnaseq_coverage.pdf \
ASSUME_SORTED=true

awk '$3=="transcript"' $path1/transcripts.gtf | sort -n -k 1 -n -k 4,5 > $path1/BTBC183_primary.transcripts.gtf

awk '$3=="exon"' $path1/transcripts.gtf | sort -n -k 1 -n -k 4,5 > $path1/BTBC183_primary.exon.gtf

htseq-count\
 /scratch/ngs/3.analysis/Project_C30/Sample_C30_0013/Sample_C30_0013.sort.bam \
/scratch/ngs/resources/GRCh37/genomes/Homo_sapiens.GRCh37.61.gtf \
-f bam -r name -s reverse -i gene_name > /scratch/ngs/3.analysis/Project_C30/Sample_C30_0013/Sample_C30_0013.htseq.txt
