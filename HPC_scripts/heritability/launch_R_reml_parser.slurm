#!/bin/bash
#SBATCH --job-name=maize_parse_LDAK_fungal_h2
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/LDAK/maize_parse_LDAK_fungal_h2%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/LDAK/maize_parse_LDAK_fungal_h2%j.err



ml anaconda
conda activate GWAS
ml r/4.5.0

Rscript \
 /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/heritability/parse_LDAK_reml.R \
 /scratch/cfc5873/maize_ldak_h2_results \
 /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/heritability/maize_fungal_h2_summary.csv
