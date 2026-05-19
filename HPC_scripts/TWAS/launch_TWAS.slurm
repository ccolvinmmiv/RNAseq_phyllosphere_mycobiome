#!/usr/bin/bash
#SBATCH --array=2
#SBATCH --job-name=TWAS_maize_688_genos
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --partition=open
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/TWAS/maize_TWAS%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/TWAS/maize_TWAS%j.err


ml r/4.5.0
ml anaconda
conda activate GWAS

USE_OPENMP=1

# List of subset files for TWAS (full paths)
pheno_file_list="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/TWAS/input_files/pheno_file_path_list.txt"

p=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${pheno_file_list} |  awk '{print $1}'`
expression="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/TWAS/input_files/maize_688_counts_tpm.csv"
gene_location_file="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/TWAS/input_files/maize_gene_positions.csv"
  

echo $p
echo $g
echo $gene_location_file

cd /scratch/cfc5873
mkdir -p maize_all_genos_mycobiome_TWAS_${SLURM_ARRAY_TASK_ID}
cd maize_all_genos_mycobiome_TWAS_${SLURM_ARRAY_TASK_ID}

cp ${p} .
cp ${gene_location_file} .
cp ${expression} .
cp /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/TWAS/GAPIT_TWAS.R .

ls -lht

pushd /scratch/cfc5873/maize_all_genos_mycobiome_TWAS_${SLURM_ARRAY_TASK_ID}
R CMD BATCH -${p} -${expression} -${gene_location_file} GAPIT_TWAS.R maize_GAPIT_TWAS_${SLURM_ARRAY_TASK_ID}.Rout


cp TWAS.CMLM_*.csv /scratch/cfc5873/mycobiome_TWAS_outputs/all_genes/maize

