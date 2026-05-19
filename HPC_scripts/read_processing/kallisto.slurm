#!/bin/bash
#SBATCH --array=1-25
#SBATCH --job-name=maize_kallisto
#SBATCH --time=48:00:00
#SBATCH --partition=open
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/kallisto/maize_kallisto_%A_%a.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/kallisto/maize_kallisto_%A_%a.err

ml anaconda
conda activate GWAS  #kallisto version 0.46.0

TRIM_DIR="/scratch/cfc5873/ne_2020_maize_trimmed_reads"
SAMPLE_LIST="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/maize_raw_samplesheet.txt"

# block size
BLOCK=30

# Unified output directory
BASE_OUT="/scratch/cfc5873/kallisto_out/ne_2020_maize"

# 
mkdir -p "$BASE_OUT"

# start & end indices for this array task
start=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * BLOCK ))

# loop through the samples for this task
for i in $(seq $start $end); do

    sample=$(sed -n "${i}p" "$SAMPLE_LIST")
    [[ -z "$sample" ]] && continue

    fq1="${TRIM_DIR}/${sample}_1_paired.fastq.gz"
    fq2="${TRIM_DIR}/${sample}_2_paired.fastq.gz"

    outdir="${BASE_OUT}/${sample}"
    mkdir -p "$outdir"


    echo "Running kallisto on sample $sample"

    kallisto quant --threads=${SLURM_CPUS_PER_TASK} \
        -i /storage/work/cfc5873/Mycobiome/kallisto_index/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.55.transcript_primaryTranscriptOnly.idx \
        -o "$outdir" \
        "$fq1" "$fq2"

done