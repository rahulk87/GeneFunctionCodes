#!/bin/bash
#BSUB -u rkumar\@icr.ac.uk
#BSUB -J crispr_D48
#BSUB -e /scratch/DBC/GENFUNC/CRISPR_SCREENS/analysis_results/D48/crispr.err
#BSUB -o /scratch/DBC/GENFUNC/CRISPR_SCREENS/analysis_results/D48/crispr.out
#BSUB -n 8
#BSUB -R "span[ptile=16]"
#BSUB -P DBCDOBZAK
#BSUB -W 96:00
#BSUB -q normal


perl crispr-align.pl \
--fastq_file_path /scratch/DBC/GENFUNC/CRISPR_SCREENS/raw_data/D48/ \
--library_file_path /scratch/DBC/GENFUNC/CRISPR_SCREENS/crispr_libraries/Human_genome_library_guides_for_shalign.txt \
--project_name D48 \
--resources_file /scratch/DBC/GENFUNC/CRISPR_SCREENS/analysis_results/D48/resources_current.txt \
--email rkumar\@icr.ac.uk
