## This scipt call somatic variants using Mutect1 and Mutect2
## You may need to change the file path, software module version etc.

#!/bin/bash
#BSUB -u rahul.kumar@icr.ac.uk
#BSUB -J mutect_PARP1
#BSUB -e /scratch/DBC/GENFUNC/NGS_Projects/SUM149_Inger_WES/SUM149_PARP1/mutect/mutect.err
#BSUB -o /scratch/DBC/GENFUNC/NGS_Projects/SUM149_Inger_WES/SUM149_PARP1/mutect/mutect.out
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -P DBCDOBZAK
#BSUB -W 96:00
#BSUB -q normal
#

path1="/scratch/DBC/GENFUNC/NGS_Projects/SUM149_Inger_WES/SUM149_PARENTAL/gatk"
path2="/scratch/DBC/GENFUNC/NGS_Projects/SUM149_Inger_WES/SUM149_PARP1"
path3="/scratch/DBC/GENFUNC/NGS_Projects/genomes" # genomes

## run Mutect PARENTAL vs PARP1
/apps/java/jdk1.7.0_80/bin/java -Xmx32g -jar /apps/mutect/1.1.7/mutect-1.1.7.jar \
--analysis_type MuTect \
--reference_sequence $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
--cosmic $path3/b37_cosmic_v54_120711.vcf \
--dbsnp $path3/dbsnp_132.b37.vcf \
--input_file:normal $path1/SUM149_PARENTAL.recaliberated.bam \
--input_file:tumor $path2/gatk/SUM149_PARP1.recaliberated.bam \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
--out $path2/mutect/PARENTAL.PARP1.mutect.txt \
--vcf $path2/mutect/PARENTAL.PARP1.mutect.vcf \
--coverage_file $path2/mutect/PARENTAL.PARP1.coverage.wig.txt

## run Mutect2 PARENTAL vs PARP1
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.5-0/GenomeAnalysisTK.jar \
-T MuTect2 \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I:normal $path1/SUM149_PARENTAL.recaliberated.bam \
-I:tumor $path2/gatk/SUM149_PARP1.recaliberated.bam \
--dbsnp $path3/dbsnp_132.b37.vcf \
--cosmic $path3/b37_cosmic_v54_120711.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path2/mutect/PARENTAL.PARP1.mutect2.vcf
