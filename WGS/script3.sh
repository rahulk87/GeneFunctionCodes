## This shell script uses Mutect2 output and run SNPEff to find out protein coding changes.

#!/bin/bash
#BSUB -u rahul.kumar@icr.ac.uk
#BSUB -J processing_BTBC183
#BSUB -e /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/mutect_calling/proc.err
#BSUB -o /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky_BTBC/Exome/analysis/BTBC183/mutect_calling/proc.out
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -P DBCDOBZAK
#BSUB -W 96:00
#BSUB -q normal

####
#Primary
####

#snpEff
/apps/java/jdk1.7.0_80/bin/java -Xmx60g -jar ~/snpEff/snpEff.jar GRCh37.75  BTBC183.germline.primary.mutect2.vcf > BTBC183.germline.primary.mutect2.snpeff.vcf

grep -w 'PASS'  BTBC183.germline.primary.mutect2.snpeff.vcf > BTBC183.germline.primary.mutect2.snpeff.pass.vcf

perl /scratch/DBC/GENFUNC/NGS_Projects/scripts/get_allele_freqs_mod.pl --vcf BTBC183.germline.primary.mutect2.snpeff.pass.vcf --columns 9,10 --out BTBC183.germline.primary.mutect2.snpeff.pass.vaf.vcf

cut -f1-5 BTBC183.germline.primary.mutect2.snpeff.pass.vaf.vcf >a
cut -f10,12-17 BTBC183.germline.primary.mutect2.snpeff.pass.vaf.vcf >b
cut -f8 BTBC183.germline.primary.mutect2.snpeff.pass.vaf.vcf | cut -d '|' -f2,3,4,5,10,11 |perl -pi -e 's/\|/\t/g' >c
paste a b c |grep -v '#' >BTBC183.germline.primary.mutect2.snpeff.pass.vaf.filtered.vcf
rm a b c

## snpeff output annotation with CGC status
perl /scratch/DBC/GENFUNC/NGS_Projects/scripts/annotate_snpeff_output_mutect2.pl BTBC183.germline.primary.mutect2.snpeff.pass.vaf.filtered.vcf > BTBC183.germline.primary.mutect2.snpeff.pass.vaf.filtered.ann.vcf


####
#PDX
####

#/apps/java/jdk1.7.0_80/bin/java -Xmx60g -jar ~/snpEff/snpEff.jar GRCh37.75  BTBC183.germline.pdx.mutect2.vcf > BTBC183.germline.pdx.mutect2.snpeff.vcf

grep -w 'PASS'  BTBC183.germline.pdx.mutect2.snpeff.vcf > BTBC183.germline.pdx.mutect2.snpeff.pass.vcf

perl /scratch/DBC/GENFUNC/NGS_Projects/scripts/get_allele_freqs_mod.pl --vcf BTBC183.germline.pdx.mutect2.snpeff.pass.vcf --columns 9,10 --out BTBC183.germline.pdx.mutect2.snpeff.pass.vaf.vcf

cut -f1-5 BTBC183.germline.pdx.mutect2.snpeff.pass.vaf.vcf >a
cut -f10,12-17 BTBC183.germline.pdx.mutect2.snpeff.pass.vaf.vcf >b
cut -f8 BTBC183.germline.pdx.mutect2.snpeff.pass.vaf.vcf | cut -d '|' -f2,3,4,5,10,11 |perl -pi -e 's/\|/\t/g' >c
paste a b c |grep -v '#' >BTBC183.germline.pdx.mutect2.snpeff.pass.vaf.filtered.vcf
rm a b c

## snpeff output annotation with CGC status
perl /scratch/DBC/GENFUNC/NGS_Projects/scripts/annotate_snpeff_output_mutect2.pl BTBC183.germline.pdx.mutect2.snpeff.pass.vaf.filtered.vcf > BTBC183.germline.pdx.mutect2.snpeff.pass.vaf.filtered.ann.vcf
