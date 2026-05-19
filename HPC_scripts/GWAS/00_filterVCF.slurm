#!/usr/bin/bash
#SBATCH --job-name=filterVCF
#SBATCH --partition=open
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/filter_VCF/maize_filter_VCF%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/filter_VCF/maize_filter_VCF%j.err

ml anaconda
conda activate GWAS

#vcftools \
 #--gzvcf /scratch/cfc5873/WiDiv.vcf.gz  \
 #--out /scratch/cfc5873/Maize_WiDiv_GeneticMarkers_filtered_688_genos \
 #--keep /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/maize_final_filtered_688_GWAS_genotypes_list.txt \
 #--recode

#python \
 #/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/HomoMAFfiltering2.py \
 #/scratch/cfc5873/Maize_WiDiv_GeneticMarkers_filtered_688_genos.recode.vcf \
 #/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.recode.HomoMAF.vcf \
 #0.05 0.1

#Maize Marker set very high density (~10 million SNPs) after filtering above, so pruning down to SNP every 3 kb
#gzip /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.recode.HomoMAF.vcf

vcftools \
  --gzvcf /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.recode.HomoMAF.vcf.gz \
  --thin 1000 \
  --recode \
  --out /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb


echo "Final SNP count:"
grep -vc "^#" /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb.recode.vcf
