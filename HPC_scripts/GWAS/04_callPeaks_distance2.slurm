#!/bin/bash
#SBATCH --array=0-9
#SBATCH --job-name=maize_688_call_peaks_eQTL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=01:00:00
#SBATCH --partition=open
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/GWAS/peak_calling/maize_688_call_peaks_eQTL%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/GWAS/peak_calling/maize_688_call_peaks_eQTL%j.err


ml anaconda
conda activate GWAS

INPUT_DIR="/scratch/cfc5873/mycobiome_eQTL_outputs/maize/significant_SNPs"
OUTPUT_DIR="/scratch/cfc5873/mycobiome_eQTL_outputs/maize/significant_SNPs/called_peaks"
mkdir -p "$OUTPUT_DIR"

# CONFIG
N_CHUNKS=10   # must match --array=0-(N_CHUNKS-1)

# Gather files
mapfile -t files < <(ls ${INPUT_DIR}/*_signals.csv)
N_FILES=${#files[@]}

if [[ $N_FILES -eq 0 ]]; then
  echo "No input files found"
  exit 1
fi

# Compute chunk boundaries
CHUNK_SIZE=$(( (N_FILES + N_CHUNKS - 1) / N_CHUNKS ))

START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
END=$(( START + CHUNK_SIZE - 1 ))

if [[ $START -ge $N_FILES ]]; then
  echo "No files assigned to task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

if [[ $END -ge $N_FILES ]]; then
  END=$(( N_FILES - 1 ))
fi

echo "Task $SLURM_ARRAY_TASK_ID processing files $START to $END"

# Process files in this chunk
for (( i=$START; i<=$END; i++ )); do
  trait_file="${files[$i]}"
  trait=$(basename "$trait_file" _signals.csv)

  echo "Processing $trait"

  python /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/call_peaks1.py \
    "$trait_file" \
    "${OUTPUT_DIR}/called_peaks_${trait}.csv" \
    "$trait" \
    1000000 \
    3
done