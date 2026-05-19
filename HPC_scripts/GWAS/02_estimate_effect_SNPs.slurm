#!/usr/bin/bash
#SBATCH --job-name=Est_Effect_SNPs_maize_688
#SBATCH --partition=open
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/Est_Effect_SNPs/Est_Effect_SNPs_maize_688%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/Est_Effect_SNPs/Est_Effect_SNPs_maize_688%j.err

ml anaconda
ml jdk/1.8.0
conda activate GWAS


java \
 -Xmx30g -jar /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/gec.jar \
 --effect-number \
 --plink-binary /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb \
 --genome \
 --out /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/maize_688_thin_1kb_gec
