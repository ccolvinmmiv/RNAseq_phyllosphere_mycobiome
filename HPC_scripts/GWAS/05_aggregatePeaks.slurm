#!/usr/bin/bash
#SBATCH --job-name=aggregate_maize_688_peaks
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=04:00:00
#SBATCH --partition=open
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/GWAS/aggregate_peaks/aggregate_maize_688_peaks_eQTL%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/GWAS/aggregate_peaks/aggregate_maize_688_peaks_eQTL%j.err

ml r/4.5.0
ml anaconda
conda activate GWAS

R CMD BATCH -/scratch/cfc5873/mycobiome_eQTL_outputs/maize/significant_SNPs/called_peaks/called_peaks*.csv -/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/all_maize_688_genos_eQTL_peaks.csv /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/aggregatePeaks2.R /storage/work/cfc5873/Mycobiome/logs/aggregate_peaks/maize_aggregatePeaks2.Rout
