library(tidyverse)

setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")


#######################      Sorghum      ######################
sorgh_v5_tf_list <- read_csv("expression_data/for_coexpression_analysis/Sobic_V5_TF_list.csv") %>%
  pull(V5_gene_ID)

sorghum_TWAS_sig_genes <- read_csv("TWAS_results/sig_genes/sorghum_736_TWAS_all_sig_hits.csv") #%>%
  pull(SNP)
sorghum_TWAS_sig_genes <- sorghum_TWAS_sig_genes %>%
  filter(SNP %in% sorgh_v5_tf_list)
length(unique(sorghum_TWAS_sig_genes))

sorghum_log_CPM_TMM_all_genes <- read_csv("expression_data/normalized_counts/sorghum_logCPM_TMM_full.csv")

sorghum_log_CPM_TMM_TWAS_genes <- sorghum_log_CPM_TMM_all_genes %>%
  select(GENOTYPE, any_of(sorghum_TWAS_sig_genes)) %>%
  mutate(GENOTYPE = str_remove(GENOTYPE, "^4\\d{3}_"))

write_csv(sorghum_log_CPM_TMM_TWAS_genes, "eQTL/input_exp_data/sorghum_log_CPM_TMM_TWAS_genes.csv")





#######################      Maize      ######################

final_widiv_key <- read_csv("genotype_files/final_widiv_genotypes_fix_key.csv") %>%
  select(!correct) %>%
  mutate(wgs = str_replace(wgs, "DK83IBI", "DK83IBI3"))

maize_TWAS_sig_genes <- read_csv("TWAS_results/sig_genes/maize_688_TWAS_all_sig_hits.csv") %>%
  pull(SNP)

length(unique(maize_TWAS_sig_genes))

maize_log_CPM_TMM_all_genes <- read_csv("expression_data/normalized_counts/maize_logCPM_TMM_full.csv")

maize_log_CPM_TMM_TWAS_genes <- maize_log_CPM_TMM_all_genes %>%
  select(GENOTYPE, any_of(maize_TWAS_sig_genes)) %>%
  left_join(final_widiv_key, join_by(GENOTYPE==orig)) %>%
  select(!GENOTYPE) %>%
  relocate(wgs) %>%
  rename(GENOTYPE = wgs)

write_csv(maize_log_CPM_TMM_TWAS_genes, "eQTL/input_exp_data/maize_log_CPM_TMM_TWAS_genes.csv")






#######################      Soybean      ######################


soybean_TWAS_sig_genes <- read_csv("TWAS_results/sig_genes/soybean_620_TWAS_all_sig_hits.csv") %>%
  pull(SNP)

length(unique(soybean_TWAS_sig_genes))

soybean_log_CPM_TMM_all_genes <- read_csv("expression_data/normalized_counts/soybean_logCPM_TMM_full.csv")

soybean_log_CPM_TMM_TWAS_genes <- soybean_log_CPM_TMM_all_genes %>%
  select(GENOTYPE, any_of(soybean_TWAS_sig_genes))

write_csv(soybean_log_CPM_TMM_TWAS_genes, "eQTL/input_exp_data/soybean_log_CPM_TMM_TWAS_genes.csv")



