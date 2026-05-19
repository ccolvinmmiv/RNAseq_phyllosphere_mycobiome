#!/bin/bash
#SBATCH --job-name=maize_LDAK_fungal_h2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/LDAK/maize_LDAK_fungal_h2%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/LDAK/maize_LDAK_fungal_h2%j.err

ml anaconda  
conda activate GWAS

BFILE=/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb
PHENO=/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/heritability/inputs/maize_LDAK_phenotypes.txt
OUTDIR=/scratch/cfc5873/maize_ldak_h2_results
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p ${OUTDIR}


# Kinship matrix
ldak6 --calc-kins-direct ${OUTDIR}/maize_data_kins \
      --bfile ${BFILE} \
      --ignore-weights YES \
      --power -0.25

# REML (multi-phenotype)
# Read trait names into bash array (assumes header row exists)
read -a traits <<< $(head -n1 $PHENO | tr "\t" " ")

# Traits start at column 3 (after FID IID)
for ((i=3; i<=${#traits[@]}; i++)); do
    trait_name=${traits[i-1]}
    echo "Running REML for $trait_name (mpheno $i)"
    ldak6 --reml $OUTDIR/maize_fungal_${trait_name}_reml \
          --grm ${OUTDIR}/maize_data_kins \
          --pheno $PHENO \
          --mpheno $i \
          --constrain YES
done

echo "LDAK REML completed successfully."
