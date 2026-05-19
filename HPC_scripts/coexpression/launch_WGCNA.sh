#!/bin/bash
#SBATCH --job-name=maize_WGCNA
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=open
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/coexpression/maize_WGCNA%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/coexpression/maize_WGCNA%j.err

ml anaconda
conda activate WGCNA
ml r/4.5.0

Rscript /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/WGCNA.R


