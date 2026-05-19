#!/bin/bash
#SBATCH --job-name=ne_maize_download
#SBATCH --partition=open
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/ne_maize_download%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/ne_maize_download%j.err


module purge
module load anaconda

# mkdir -p batches
# split -l 200 ne_maize_file_list.txt batches/batch_

# Set scratch folder
SCRATCH_DIR=/scratch/cfc5873/ne_2020_maize_raw_reads
mkdir -p $SCRATCH_DIR
cd $SCRATCH_DIR



# File containing URLs for this batch (replace below when submitting)
BATCH_FILE=$1

# Download in parallel using 4 threads

cat $BATCH_FILE | xargs -n 1 -P $SLURM_CPUS_PER_TASK wget -c -t 5 --retry-connrefused --waitretry=10



# run with this: for f in /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/batches/batch_*; do sbatch /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/ne_maize_download_template.sh "$f"; done
