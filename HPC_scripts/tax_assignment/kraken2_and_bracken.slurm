#!/bin/bash
#SBATCH --job-name=maize_kraken2_bracken
#SBATCH --array=10-25   # open only allows 200 jobs at a time per user
#SBATCH --mem=400G
#SBATCH --account=sic3_cr_default     #Credit account
#SBATCH --partition=himem
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=36:00:00
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/kraken/maize_kraken_bracken%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/kraken/maize_kraken_bracken%j.err



# Load environment. GWAS Env contains both kraken2 and bracken, but we will use local version of kraken2 since it is more up to date
ml anaconda
conda activate GWAS

# Paths
KRAKENDB="/scratch/cfc5873/kraken2_db"
TRIM_DIR="/scratch/cfc5873/ne_2020_maize_trimmed_reads"
SAMPLE_LIST="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/maize_raw_samplesheet.txt"

# Use local kraken2 (Version 2.1) instead of the cluster version
KRAKEN2_BIN="/storage/home/cfc5873/kraken2/kraken2"

# Kraken outputs
KRAKEN_OUT="/scratch/cfc5873/kraken_out/ne_2020_maize"
mkdir -p "$KRAKEN_OUT"

# Bracken outputs
BRACKEN_OUT="/scratch/cfc5873/bracken_out/ne_2020_maize"
mkdir -p "$BRACKEN_OUT"

# Determine workload for this array task
total_samples=$(wc -l < "$SAMPLE_LIST")
num_jobs=25

samples_per_job=$(( (total_samples + num_jobs - 1) / num_jobs ))
start=$(( (SLURM_ARRAY_TASK_ID - 1) * samples_per_job + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * samples_per_job ))
if [ "$end" -gt "$total_samples" ]; then end=$total_samples; fi

echo "Task $SLURM_ARRAY_TASK_ID processing samples $start to $end"

# Run each assigned sample
for i in $(seq $start $end); do
    sample=$(sed -n "${i}p" "$SAMPLE_LIST")
    echo "Processing $sample"

    READ1="${TRIM_DIR}/${sample}_1_paired.fastq.gz"
    READ2="${TRIM_DIR}/${sample}_2_paired.fastq.gz"

    REPORT="${KRAKEN_OUT}/${sample}.report"
    OUTPUT="${KRAKEN_OUT}/${sample}.kraken"

   
    #KRAKEN2 (using local version 2.1)
    "$KRAKEN2_BIN" \
        --db "$KRAKENDB" \
        --threads $SLURM_CPUS_PER_TASK \
        --paired \
        --report "$REPORT" \
        --output "$OUTPUT" \
        "$READ1" "$READ2"


    #BRACKEN
    # Levels to compute: Class Order Family Genus Species
    for LEVEL in C O F G S; do
        bracken \
            -d "$KRAKENDB" \
            -i "$REPORT" \
            -o "${BRACKEN_OUT}/${sample}.bracken.${LEVEL}" \
            -r 150 \
            -l $LEVEL
    done

    echo "Completed $sample"
done