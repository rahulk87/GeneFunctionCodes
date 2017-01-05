## This script run cuffmerge and cuffdiff for the differential expression analysis
## You may need to change the file path, software module version etc.

#!/bin/bash
#BSUB -u rahul.kumar\@icr.ac.uk
#BSUB -J cuffmerge_BTBC183
#BSUB -e /scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/analysis/BTBC183/compare/cuffmerge.err
#BSUB -o /scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/analysis/BTBC183/compare/cuffmerge.out
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

path1="/scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/analysis/BTBC183/compare"
path2="/scratch/DBC/GENFUNC/rkumar/NGS_projects/Becky_BTBC/analysis/BTBC183/"
path3="/scratch/DBC/GENFUNC/rkumar/NGS_projects/genomes"

## Cuffmerge
cuffmerge \
-g $path3/Homo_sapiens.GRCh37.82.chr.gtf \
-s $path3/ens61_human_chromosomes_and_MT.fa \
-p 8 \
assemblies.txt

## Cuffdiff
cuffdiff \
-o diff_out \
-b $path3/ens61_human_chromosomes_and_MT.fa \
-p 8 \
-L BTBC183_primary,BTBC183_pdx \
--library-type fr-firststrand \
-u merged_asm/merged.gtf \
$path2/primary/accepted_hits.reorder.sort.bam \
$path2/pdx/accepted_hits.reorder.sort.bam

## Cuffcompare
cuffcompare \
-o compare \
-s $path3/ens61_human_chromosomes_and_MT.fa \
-r $path3/Homo_sapiens.GRCh37.82.chr.gtf  \
$path2/primary/guided/transcripts.gtf \
$path2/pdx/guided/transcripts.gtf 
