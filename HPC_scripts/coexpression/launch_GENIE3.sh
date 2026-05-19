#!/bin/bash
#SBATCH --job-name=GENIE3_maize
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --partition=open
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/coexpression/maize_GENIE3_%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/coexpression/maize_GENIE3_%j.err

ml anaconda
conda activate GENIE
ml r/4.5.0

# Run the R script


Rscript /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/GENIE3_all_genes.R   
