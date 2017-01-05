## This scipt call somatic variants using Mutect1 and Mutect2

#!/bin/bash
#BSUB -u rahul.kumar@icr.ac.uk
#BSUB -J BTBC183_mutect
#BSUB -e /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/mutect_calling/mutect.err
#BSUB -o /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/mutect_calling/mutect.out
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -P DBCDOBZAK
#BSUB -W 96:00
#BSUB -q normal
#

path1="/scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/mutect_calling"
path2="/scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183"
path3="/scratch/DBC/GENFUNC/NGS_Projects/genomes" # genomes

## run Mutect germline vs primary
/apps/java/jdk1.7.0_80/bin/java -Xmx32g -jar /apps/mutect/1.1.7/mutect-1.1.7.jar \
--analysis_type MuTect \
--reference_sequence $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
--cosmic $path3/b37_cosmic_v54_120711.vcf \
--dbsnp $path3/dbsnp_132.b37.vcf \
--input_file:normal $path2/germline/BTBC183_germline.recaliberated.bam \
--input_file:tumor $path2/primary/BTBC183_primary.recaliberated.bam \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
--out $path1/BTBC183.germline.primary.mutect.txt \
--vcf $path1/BTBC183.germline.primary.mutect.vcf \
--coverage_file $path1/BTBC183.germline.primary.mutect.coverage.wig.txt

## run Mutect germline vs pdx
/apps/java/jdk1.7.0_80/bin/java -Xmx32g -jar /apps/mutect/1.1.7/mutect-1.1.7.jar \
--analysis_type MuTect \
--reference_sequence $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
--cosmic $path3/b37_cosmic_v54_120711.vcf \
--dbsnp $path3/dbsnp_132.b37.vcf \
--input_file:normal $path2/germline/BTBC183_germline.recaliberated.bam \
--input_file:tumor $path2/pdx/BTBC183_pdx.recaliberated.bam \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
--out $path1/BTBC183.germline.pdx.mutect.txt \
--vcf $path1/BTBC183.germline.pdx.mutect.vcf \
--coverage_file $path1/BTBC183.germline.pdx.mutect.coverage.wig.txt

## run Mutect2 germline vs primary
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.5-0/GenomeAnalysisTK.jar \
-T MuTect2 \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I:normal $path2/germline/BTBC183_germline.recaliberated.bam \
-I:tumor $path2/primary/BTBC183_primary.recaliberated.bam \
--dbsnp $path3/dbsnp_132.b37.vcf \
--cosmic $path3/b37_cosmic_v54_120711.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path1/BTBC183.germline.primary.mutect2.vcf

## run Mutect2 germline vs pdx
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.5-0/GenomeAnalysisTK.jar \
-T MuTect2 \
-R $path3/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I:normal $path2/germline/BTBC183_germline.recaliberated.bam \
-I:tumor $path2/pdx/BTBC183_pdx.recaliberated.bam \
--dbsnp $path3/dbsnp_132.b37.vcf \
--cosmic $path3/b37_cosmic_v54_120711.vcf \
-L $path3/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path1/BTBC183.germline.pdx.mutect2.vcf
