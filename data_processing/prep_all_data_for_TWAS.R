library("data.table")
library(tidyverse)


setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")

#################### Sorghum Prep ############################


#load counts/gene expression data
sorghum_raw_tpm_data <- read_csv("expression_data/raw_expression_data/ne_2021_sorgh_merged_gene_tpms.csv")

sorgh_used_samples_names <- read_csv("genotype_files/ne_sorgh_2021_representative_samples_by_TPMs.csv")
sorgh_used_genos_736 <- read_table("genotype_files/sorgh_final_filtered_736_GWAS_genotypes_list.txt", col_names = "GENOTYPE")


sorgh_genos_to_keep_for_TWAS <- sorgh_used_samples_names %>%
  filter(GENOTYPE %in% sorgh_used_genos_736$GENOTYPE) %>%
  pull(SAMPLE)

sorghum_raw_tpm_data_correct_samples <- sorghum_raw_tpm_data %>% #Need to subset to the same 736 genotypes/samples used for GWAS/TWAS
  dplyr::select(TranscriptID, any_of(sorgh_genos_to_keep_for_TWAS))

counts <- sorghum_raw_tpm_data_correct_samples %>%
  rename_with(.fn = ~str_replace(.x,"4..._", "")) %>%
  pivot_longer(!TranscriptID, names_to = "GENOTYPE", values_to = "TPM") %>%
  pivot_wider(names_from = TranscriptID, values_from = TPM)

write_csv(counts, "data/TWAS_gene_TPM_files/sorghum_736_counts_tpm.csv")
sorgh_counts <- read_csv("data/TWAS_gene_TPM_files/sorghum_736_counts_tpm.csv")


# Load/compute gene info (position)
GFF <- fread("data/gene_data/Sbicolor_730_v5.1.gene.gff3", skip = "#", header = FALSE)
colnames(GFF) <- c(
  "seqid", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)
GFF[, gene_id := sub(".*ID=([^;]+).*", "\\1", attributes)]
GFF <- GFF %>%
  filter(type == "gene") %>%
  mutate(
    pos = floor((start + end) / 2), #using midpoint of each gene as the position
    chr = as.numeric(str_remove(seqid, "Chr")),
    gene = str_remove(gene_id, ".v5.1")
  ) %>%
  dplyr::select(gene, chr, pos)

write_csv(GFF, "data/gene_data/sorghum_gene_positions.csv")











#################### Soybean Prep ############################

#load counts/gene expression data
soybean_raw_tpm_data <- read_csv("expression_data/raw_expression_data/soybean_merged_gene_tpms.csv")

#Getting Soybean Genotype Names
soybean_run_metadata <- read_csv("genotype_files/soybean_run_metadata.csv")
soybean_sample_metadata <- read_csv("genotype_files/soybean_sample_metadata.csv")

soybean_run_sample_combined_metadata <- soybean_run_metadata %>%
  left_join(soybean_sample_metadata, join_by(`Run title` == `Sample name`)) %>%
  select(Accession.x, Cultivar) %>%
  rename(RunID = Accession.x,
         GENOTYPE_REAL = Cultivar)

#Need to average accross two genotypes that each have two samples (goes from 622 total samples/620 genotypes to 620 total samples and 620 unique genotypes)
dup_ids <- soybean_run_sample_combined_metadata %>%
  count(GENOTYPE_REAL) %>%
  filter(n > 1) %>%
  pull(GENOTYPE_REAL)

expr_wide <- soybean_raw_tpm_data %>%
  pivot_longer(cols = -TranscriptID, names_to = "RunID", values_to = "expression_count") %>%
  pivot_wider(id_cols = RunID, names_from = TranscriptID, values_from = expression_count) %>%
  left_join(soybean_run_sample_combined_metadata, join_by(RunID)) %>%
  rename(Sample = GENOTYPE_REAL) %>%
  select(Sample, where(is.numeric))

expr_single <- expr_wide %>%
  filter(!Sample %in% dup_ids)

expr_dup <- expr_wide %>%
  filter(Sample %in% dup_ids)

expr_dup_collapsed <- expr_dup %>%
  group_by(Sample) %>%
  summarise(across(everything(), ~mean(.x, na.rm = FALSE)), .groups = "drop")

final_expr <- bind_rows(expr_single, expr_dup_collapsed)

soybean_counts <- final_expr %>%
  rename(GENOTYPE = Sample)

write_csv(soybean_counts, "data/TWAS_gene_TPM_files/soybean_620_counts_tpm.csv")


# Load/compute gene info (position)
GFF <- fread("data/gene_data/Gmax_275_Wm82.a2.v1.gene.gff3", skip = "#", header = FALSE)
colnames(GFF) <- c(
  "seqid", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)
GFF[, gene_id := sub(".*ID=([^;]+).*", "\\1", attributes)]
GFF <- GFF %>%
  filter(type == "gene") %>%
  mutate(
    pos = floor((start + end) / 2), #using midpoint of each gene as the position
    chr = as.numeric(str_remove(seqid, "Chr")),
    gene = str_remove(gene_id, ".Wm82.a2.v1")
  ) %>%
  dplyr::select(gene, chr, pos)

write_csv(GFF, "data/gene_data/soybean_gene_positions.csv")







#################### Maize Prep ############################

#load counts/gene expression data
maize_raw_tpm_data <- read_csv("expression_data/raw_expression_data/ne_2020_maize_merged_gene_tpms.csv")

final_widiv_key <- read_csv("genotype_files/final_widiv_genotypes_fix_key.csv") %>%
  select(!correct) %>%
  mutate(wgs = str_replace(wgs, "DK83IBI", "DK83IBI3"))

GWAS_688_Genos_wgs_names <- read_table("genotype_files/maize_final_filtered_688_GWAS_genotypes_list.txt", col_names = "WGS")

genos_to_keep_for_coexp <- final_widiv_key %>%
  filter(wgs %in% GWAS_688_Genos_wgs_names$WGS) %>%
  pull(orig)

maize_raw_count_data_correct_samples <- maize_raw_tpm_data %>% #Need to subset to the same 688 genotypes/samples used for GWAS/TWAS
  rename_with(toupper) %>%
  select(TRANSCRIPTID, any_of(genos_to_keep_for_coexp)) %>%
  rename(TranscriptID = TRANSCRIPTID)

maize_counts <- maize_raw_count_data_correct_samples %>%
  pivot_longer(!TranscriptID, names_to = "GENOTYPE", values_to = "TPM") %>%
  pivot_wider(names_from = TranscriptID, values_from = TPM) %>%
  rename_with(~ str_remove(.x, "_T\\d+$")) 
# ^ Remove transcript indicator for downstream (Soybean and Sorghum use a different convention and are removed in main TWAS script)
# Doing this now for maize allows use of same TWAS script for everything without adjustment

write_csv(maize_counts, "data/TWAS_gene_TPM_files/maize_688_counts_tpm.csv")


# Load/compute gene info (position)
GFF <- fread("data/gene_data/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.55.gene.gff3", skip = "#", header = FALSE)
colnames(GFF) <- c(
  "seqid", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)
GFF[, gene_id := sub(".*ID=([^;]+).*", "\\1", attributes)]
GFF <- GFF %>%
  filter(type == "gene") %>%
  mutate(
    pos = floor((start + end) / 2), #using midpoint of each gene as the position
    chr = as.numeric(str_remove(seqid, "Chr")),
    gene = str_remove(gene_id, ".Zm_B73_REFERENCE_NAM_5.0.55")
  ) %>%
  dplyr::select(gene, chr, pos)

write_csv(GFF, "data/gene_data/maize_gene_positions.csv")













