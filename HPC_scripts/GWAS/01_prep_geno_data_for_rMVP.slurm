#!/usr/bin/bash
#SBATCH --job-name=maize_prep_for_rMVP
#SBATCH --partition=open
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/prep_geno_data_rMVP/maize_prep_for_rMVP%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/prep_geno_data_rMVP/maize_prep_for_rMVP%j.err


ml r/4.5.0
ml anaconda
conda activate GWAS


Rscript /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/prep_geno_data_for_rMVP.R \
 /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb.recode.vcf \
 /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb.mvp \
 /storage/work/cfc5873/Mycobiome/logs/prep_geno_data_rMVP/Maize_WiDiv_GeneticMarkers_filtered_688_genos_prepGenotypic.Rout
