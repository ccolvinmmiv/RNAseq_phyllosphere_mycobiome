#!/bin/bash
#SBATCH --job-name=orthofinder_maize_sorghum_soybean
#SBATCH --account=one_sc_default
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/orthofinder/orthofinder_%j.out
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/orthofinder/orthofinder_%j.err

ml anaconda
conda activate orthofinder3

# Input directory containing protein FASTA files
INPUT_DIR="/storage/work/cfc5873/Mycobiome/protein_files/raw_fasta_files/clean_headers"

# Output directory
OUT_DIR="/storage/work/cfc5873/Mycobiome/protein_files/orthofinder_outputs"

echo "Starting OrthoFinder run"

orthofinder \
    -f $INPUT_DIR \
    -o $OUT_DIR \
    -t $SLURM_CPUS_PER_TASK \
    -a $SLURM_CPUS_PER_TASK

echo "OrthoFinder completed"
