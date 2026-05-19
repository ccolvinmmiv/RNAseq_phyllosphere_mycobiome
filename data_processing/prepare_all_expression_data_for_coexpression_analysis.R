library(tidyverse)
library(tximport)
library(edgeR)
library(fuzzyjoin)


setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")


################ Normalizing/preparing maize expression data for coexpression ################
maize_raw_count_data <- read_csv("expression_data/raw_expression_data/ne_2020_maize_merged_gene_raw_counts.csv")

final_widiv_key <- read_csv("genotype_files/final_widiv_genotypes_fix_key.csv") %>%
  select(!correct) %>%
  mutate(wgs = str_replace(wgs, "DK83IBI", "DK83IBI3"))

GWAS_688_Genos_wgs_names <- read_table("genotype_files/maize_final_filtered_688_GWAS_genotypes_list.txt", col_names = "WGS")

genos_to_keep_for_coexp <- final_widiv_key %>%
  filter(wgs %in% GWAS_688_Genos_wgs_names$WGS) %>%
  pull(orig)

maize_raw_count_data_correct_samples <- maize_raw_count_data %>% #Need to subset to the same 688 genotypes/samples used for GWAS/TWAS
  rename_with(toupper) %>%
  select(TRANSCRIPTID, any_of(genos_to_keep_for_coexp)) %>%
  rename(TranscriptID = TRANSCRIPTID)

flipped_maize_raw_count_data <- maize_raw_count_data_correct_samples %>%
  pivot_longer(cols = -TranscriptID, names_to = "Sample", values_to = "expression_count") %>%
  pivot_wider(id_cols = Sample, names_from = TranscriptID, values_from = expression_count)


maize_counts_t <- flipped_maize_raw_count_data %>%
  column_to_rownames("Sample") %>%  # set sample IDs as rownames
  t() %>%                              # transpose: now rows = transcripts
  as.data.frame() %>%
  rownames_to_column("TranscriptID")

# Aggregate to gene level (Only does anything when using all transcripts not primary only)
maize_counts_gene <- maize_counts_t %>%
  mutate(GeneID = sub("_T\\d+$", "", TranscriptID)) %>%
  group_by(GeneID) %>%
  summarise(across(-TranscriptID, sum)) %>%
  column_to_rownames("GeneID")

# Flip back to original format
maize_counts_gene <- t(maize_counts_gene) %>% 
  as.data.frame() %>%
  rownames_to_column("ID")

# Pull genotype vector
genotypes_maize <- maize_counts_gene$ID

# Expression matrix only
maize_counts <- maize_counts_gene %>%
  select(-ID) %>%
  as.data.frame()

# Make sure the row names are genotype IDs
rownames(maize_counts) <- genotypes_maize

# Create DGEList for edgeR
maize_counts_edger <- t(maize_counts)  # genes x samples
maize_dge <- DGEList(counts = maize_counts_edger)

# Remove genes with zero counts in >50% of samples
maize_keep <- rowSums(maize_dge$counts > 0) > (ncol(maize_dge$counts)/2)
maize_dge <- maize_dge[maize_keep, , keep.lib.sizes = FALSE]

# Normalization (TMM)
maize_dge <- calcNormFactors(maize_dge, method = "TMM")

# TMM logCPM
maize_logCPM_TMM <- edgeR::cpm(maize_dge, log = TRUE, prior.count = 1) %>%
  t() %>% as.data.frame()  # back to samples x genes



# Convert all matrices to tibbles with ID column
maize_logCPM_TMM_tb <- maize_logCPM_TMM %>%
  rownames_to_column("ID") %>%
  as_tibble()

# Save full normalized matrix for future reference
maize_logCPM_TMM_tb_to_save <- maize_logCPM_TMM_tb %>%
  rename(GENOTYPE = ID)

write_csv(maize_logCPM_TMM_tb_to_save, "expression_data/normalized_counts/maize_logCPM_TMM_full.csv")

### Flipping/prepping data for WGCNA co expression

maize_WGCNA_expr <- read_csv("expression_data/normalized_counts/maize_logCPM_TMM_full.csv")

maize_WGCNA_mat <- as.data.frame(maize_WGCNA_expr)
rownames(maize_WGCNA_mat) <- maize_WGCNA_mat$GENOTYPE
maize_WGCNA_mat$GENOTYPE <- NULL

saveRDS(maize_WGCNA_mat, "expression_data/for_coexpression_analysis/maize_WGCNA_ready.rds")

#Prep Expression matrix for GENIE3

maize_GENIE_mat <- as.data.frame(maize_WGCNA_expr)
rownames(maize_GENIE_mat) <- maize_GENIE_mat$GENOTYPE
maize_GENIE_mat$GENOTYPE <- NULL 
#maize_GENIE_mat <- maize_GENIE_mat %>% #If interested in particular WGCNA module, you can do filtering here before running GENIE3 to make it run faster
#  select(any_of(all_maize_brown_genes))

saveRDS(maize_GENIE_mat, "expression_data/for_coexpression_analysis/maize_GENIE3_exp_mat.rds")









################ Normalizing/preparing sorghum expression data for coexpression ################
sorghum_raw_count_data <- read_csv("expression_data/raw_expression_data/ne_2021_sorgh_merged_gene_raw_counts.csv")


sorgh_used_samples_names <- read_csv("genotype_files/ne_sorgh_2021_representative_samples_by_TPMs.csv")
sorgh_used_genos_736 <- read_table("genotype_files/sorgh_final_filtered_736_GWAS_genotypes_list.txt", col_names = "GENOTYPE")


sorgh_genos_to_keep_for_coexp <- sorgh_used_samples_names %>%
  filter(GENOTYPE %in% sorgh_used_genos_736$GENOTYPE) %>%
  pull(SAMPLE)

sorghum_raw_count_data_correct_samples <- sorghum_raw_count_data %>% #Need to subset to the same 736 genotypes/samples used for GWAS/TWAS
  select(TranscriptID, any_of(sorgh_genos_to_keep_for_coexp))





flipped_sorghum_raw_count_data <- sorghum_raw_count_data_correct_samples %>%
  pivot_longer(cols = -TranscriptID, names_to = "Sample", values_to = "expression_count") %>%
  pivot_wider(id_cols = Sample, names_from = TranscriptID, values_from = expression_count)


sorghum_counts_t <- flipped_sorghum_raw_count_data %>%
  column_to_rownames("Sample") %>%  # set sample IDs as rownames
  t() %>%                              # transpose: now rows = transcripts
  as.data.frame() %>%
  rownames_to_column("TranscriptID")

# Aggregate to gene level (Only does anything when using all transcripts not primary only)
sorghum_counts_gene <- sorghum_counts_t %>%
  mutate(GeneID = sub("\\.[0-9]+$", "", TranscriptID)) %>%
  group_by(GeneID) %>%
  summarise(across(-TranscriptID, sum)) %>%
  column_to_rownames("GeneID")

# Flip back to original format
sorghum_counts_gene <- t(sorghum_counts_gene) %>% 
  as.data.frame() %>%
  rownames_to_column("ID")

# Pull genotype vector
genotypes_sorghum <- sorghum_counts_gene$ID

# Expression matrix only
sorghum_counts <- sorghum_counts_gene %>%
  select(-ID) %>%
  as.data.frame()

# Make sure the row names are genotype IDs
rownames(sorghum_counts) <- genotypes_sorghum

# Create DGEList for edgeR
sorghum_counts_edger <- t(sorghum_counts)  # genes x samples
sorghum_dge <- DGEList(counts = sorghum_counts_edger)

# Remove genes with zero counts in >50% of samples
sorghum_keep <- rowSums(sorghum_dge$counts > 0) > (ncol(sorghum_dge$counts)/2)
sorghum_dge <- sorghum_dge[sorghum_keep, , keep.lib.sizes = FALSE]

# Normalization (TMM)
sorghum_dge <- calcNormFactors(sorghum_dge, method = "TMM")

# TMM logCPM
sorghum_logCPM_TMM <- edgeR::cpm(sorghum_dge, log = TRUE, prior.count = 1) %>%
  t() %>% as.data.frame()  # back to samples x genes



# Convert all matrices to tibbles with ID column
sorghum_logCPM_TMM_tb <- sorghum_logCPM_TMM %>%
  rownames_to_column("ID") %>%
  as_tibble()

# Save full normalized matrix for future reference
sorghum_logCPM_TMM_tb_to_save <- sorghum_logCPM_TMM_tb %>%
  rename(GENOTYPE = ID)

write_csv(sorghum_logCPM_TMM_tb_to_save, "expression_data/normalized_counts/sorghum_logCPM_TMM_full.csv")

### Flipping/prepping data for WGCNA co expression

sorghum_WGCNA_expr <- read_csv("expression_data/normalized_counts/sorghum_logCPM_TMM_full.csv")

sorghum_WGCNA_mat <- as.data.frame(sorghum_WGCNA_expr)
rownames(sorghum_WGCNA_mat) <- sorghum_WGCNA_mat$GENOTYPE
sorghum_WGCNA_mat$GENOTYPE <- NULL

saveRDS(sorghum_WGCNA_mat, "expression_data/for_coexpression_analysis/sorghum_WGCNA_ready.rds")

#Prep Expression matrix for GENIE3

sorghum_GENIE_mat <- as.data.frame(sorghum_WGCNA_expr)
rownames(sorghum_GENIE_mat) <- sorghum_GENIE_mat$GENOTYPE
sorghum_GENIE_mat$GENOTYPE <- NULL 
#sorghum_GENIE_mat <- sorghum_GENIE_mat %>% #If interested in particular WGCNA module, you can do filtering here before running GENIE3 to make it run faster
#  select(any_of(all_sorghum_brown_genes))

saveRDS(sorghum_GENIE_mat, "expression_data/for_coexpression_analysis/sorghum_GENIE3_exp_mat.rds")


#plantTFDB sorghum TF list comes in V3 naming format, have to use name conversion file from phytozome to get into V5 names

sorgh_V3_TF_list <- read_table("expression_data/for_coexpression_analysis/Sbi_TF_list.txt")
sorgh_V3_to_V5_sheet <- read_table("genotype_files/Sbicolor_730_v5.1.locus_transcript_name_map.txt")

sorgh_V3_V5_TF_list <- sorgh_V3_TF_list %>%
  left_join(sorgh_V3_to_V5_sheet, join_by(Gene_ID == oldlocusName)) %>%
  rename(V3_gene_ID = Gene_ID) %>%
  rename(V5_gene_ID = new_locusName)
write_csv(sorgh_V3_V5_TF_list, "expression_data/for_coexpression_analysis/Sobic_V5_TF_list.csv")







################ Normalizing/preparing soybean expression data for coexpression ################
soybean_raw_count_data <- read_csv("expression_data/raw_expression_data/soybean_merged_gene_raw_counts.csv")

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

expr_wide <- soybean_raw_count_data %>%
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

flipped_soybean_raw_count_data <- final_expr

soybean_counts_t <- flipped_soybean_raw_count_data %>%
  column_to_rownames("Sample") %>%  # set sample IDs as rownames
  t() %>%                              # transpose: now rows = transcripts
  as.data.frame() %>%
  rownames_to_column("TranscriptID")

# Aggregate to gene level (Only does anything when using all transcripts not primary only)
soybean_counts_gene <- soybean_counts_t %>%
  mutate(GeneID = sub("\\.[0-9]+$", "", TranscriptID)) %>%
  group_by(GeneID) %>%
  summarise(across(-TranscriptID, sum)) %>%
  column_to_rownames("GeneID")

# Flip back to original format
soybean_counts_gene <- t(soybean_counts_gene) %>% 
  as.data.frame() %>%
  rownames_to_column("ID")

# Pull genotype vector
genotypes_soybean <- soybean_counts_gene$ID

# Expression matrix only
soybean_counts <- soybean_counts_gene %>%
  select(-ID) %>%
  as.data.frame()

# Make sure the row names are genotype IDs
rownames(soybean_counts) <- genotypes_soybean

# Create DGEList for edgeR
soybean_counts_edger <- t(soybean_counts)  # genes x samples
soybean_dge <- DGEList(counts = soybean_counts_edger)

# Remove genes with zero counts in >50% of samples
soybean_keep <- rowSums(soybean_dge$counts > 0) > (ncol(soybean_dge$counts)/2)
soybean_dge <- soybean_dge[soybean_keep, , keep.lib.sizes = FALSE]

# Normalization (TMM)
soybean_dge <- calcNormFactors(soybean_dge, method = "TMM")

# TMM logCPM
soybean_logCPM_TMM <- edgeR::cpm(soybean_dge, log = TRUE, prior.count = 1) %>%
  t() %>% as.data.frame()  # back to samples x genes



# Convert all matrices to tibbles with ID column
soybean_logCPM_TMM_tb <- soybean_logCPM_TMM %>%
  rownames_to_column("ID") %>%
  as_tibble()

# Save full normalized matrix for future reference
soybean_logCPM_TMM_tb_to_save <- soybean_logCPM_TMM_tb %>%
  rename(GENOTYPE = ID)

write_csv(soybean_logCPM_TMM_tb_to_save, "expression_data/normalized_counts/soybean_logCPM_TMM_full.csv")

### Flipping/prepping data for WGCNA co expression

soybean_WGCNA_expr <- read_csv("expression_data/normalized_counts/soybean_logCPM_TMM_full.csv")

soybean_WGCNA_mat <- as.data.frame(soybean_WGCNA_expr)
rownames(soybean_WGCNA_mat) <- soybean_WGCNA_mat$GENOTYPE
soybean_WGCNA_mat$GENOTYPE <- NULL

saveRDS(soybean_WGCNA_mat, "expression_data/for_coexpression_analysis/soybean_WGCNA_ready.rds")

#Prep Expression matrix for GENIE3

soybean_GENIE_mat <- as.data.frame(soybean_WGCNA_expr)
rownames(soybean_GENIE_mat) <- soybean_GENIE_mat$GENOTYPE
soybean_GENIE_mat$GENOTYPE <- NULL 
#soybean_GENIE_mat <- soybean_GENIE_mat %>% #If interested in particular WGCNA module, you can do filtering here before running GENIE3 to make it run faster
#  select(any_of(all_soybean_brown_genes))

saveRDS(soybean_GENIE_mat, "expression_data/for_coexpression_analysis/soybean_GENIE3_exp_mat.rds")


soybean_tf_table <- read.table("expression_data/for_coexpression_analysis/Gma_TF_list.txt", header = TRUE)







