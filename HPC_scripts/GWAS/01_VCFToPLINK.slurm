#!/usr/bin/bash
#SBATCH --job-name=convertVCFToPLINK_maize_688
#SBATCH --partition=open
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/VCF_to_PLINK/convertVCFToPLINK_maize_688%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/VCF_to_PLINK/convertVCFToPLINK_maize_688%j.err


ml anaconda
conda activate GWAS


plink2 \
 --vcf /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb.recode.vcf \
 --make-bed \
 --out /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb \
 --const-fid 
 
