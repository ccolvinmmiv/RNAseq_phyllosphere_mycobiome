#!/usr/bin/bash
#SBATCH --array=1
#SBATCH --job-name=mlm_GWAS_maize_688_eQTL
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5G
#SBATCH --time=2-00:00:00
#SBATCH --partition=open
#SBATCH --output=/storage/work/cfc5873/Mycobiome/logs/GWAS/maize_688_eQTL%j.log
#SBATCH --error=/storage/work/cfc5873/Mycobiome/logs/GWAS/maize_688_eQTL%j.err


ml r/4.5.0
ml anaconda
conda activate GWAS

USE_OPENMP=1

#pheno_file_list="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/phenotype_files/pheno_subset_paths.txt"

#p=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${pheno_file_list} |  awk '{print $1}'`
p="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/phenotype_files/eQTL_exp_data/maize_log_CPM_TMM_TWAS_genes.csv"
g="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/genotype_files/Maize_WiDiv_GeneticMarkers_filtered_688_genos.HomoMAF.thin1kb.mvp"

total_markers=1006790               	
effective_markers=480364.94    

echo $p
echo $g
echo $total_markers
echo $effective_markers

cd /scratch/cfc5873
mkdir -p maize_all_genos_eQTL_${SLURM_ARRAY_TASK_ID}
cd maize_all_genos_eQTL_${SLURM_ARRAY_TASK_ID}

cp ${p} .
cp ${g}.* .
cp /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/GWAS/MLMGWAS_trycatch.R .

ls -lht

pushd /scratch/cfc5873/maize_all_genos_eQTL_${SLURM_ARRAY_TASK_ID}
R CMD BATCH -${p} -${g} -${total_markers} -${effective_markers} MLMGWAS_trycatch.R maize_all_genos_eQTL_MLM_24_CPU_${SLURM_ARRAY_TASK_ID}.Rout


cp *_signals.csv /scratch/cfc5873/mycobiome_eQTL_outputs/maize/significant_SNPs
cp *.MLM.csv /scratch/cfc5873/mycobiome_eQTL_outputs/maize/all_SNPs



