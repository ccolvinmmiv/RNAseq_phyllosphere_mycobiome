#!/bin/bash
#SBATCH --job-name=maize_trimmomatic
#SBATCH --array=1-50   # open only allows 200 jobs at a time per user
#SBATCH --partition=open
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/trimmomatic/maize_trimmomatic%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/trimmomatic/maize_trimmomatic%j.err




#First run this: 
#   ls /scratch/cfc5873/ne_2020_maize_raw_reads/*_1.fastq.gz | xargs -n1 basename | sed 's/_1\.fastq\.gz//' > maize_raw_samplesheet.txt

ml jdk/1.8.0
RAW_DIR="/scratch/cfc5873/ne_2020_maize_raw_reads"
OUT_DIR="/scratch/cfc5873/ne_2020_maize_trimmed_reads"
TRIMJAR="/storage/work/cfc5873/Mycobiome/scripts/trimmomatic-0.39.jar"
ADAPTERS="/storage/work/cfc5873/Mycobiome/scripts/TruSeq3-PE-2.fa"
SAMPLESHEET="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/maize_raw_samplesheet.txt"

mkdir -p "$OUT_DIR"

# Total samples
total_samples=$(wc -l < "$SAMPLESHEET")
# Jobs in array
num_jobs=50

# Calculate chunk for this job
samples_per_job=$(( (total_samples + num_jobs - 1) / num_jobs ))
start=$(( (SLURM_ARRAY_TASK_ID - 1) * samples_per_job + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * samples_per_job ))
if [ "$end" -gt "$total_samples" ]; then end=$total_samples; fi

# Loop over all samples assigned to this job
for i in $(seq $start $end); do
    sample=$(sed -n "${i}p" "$SAMPLESHEET")
    echo "Processing sample: $sample"

    # Define input and output files
    r1="${RAW_DIR}/${sample}_1.fastq.gz"
    r2="${RAW_DIR}/${sample}_2.fastq.gz"

    out1="${OUT_DIR}/${sample}_1_paired.fastq.gz"
    out1u="${OUT_DIR}/${sample}_1_unpaired.fastq.gz"
    out2="${OUT_DIR}/${sample}_2_paired.fastq.gz"
    out2u="${OUT_DIR}/${sample}_2_unpaired.fastq.gz"

  # Run Trimmomatic
    java -jar "$TRIMJAR" PE -phred33 \
      "$r1" "$r2" \
      "$out1" "$out1u" \
      "$out2" "$out2u" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35

    echo "Completed sample: $sample"

done