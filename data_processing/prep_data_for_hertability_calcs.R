library(tidyverse)











##########################     Sorghum     ##############################


sorghum_GWAS_pheno_file_full <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv")

sorghum_file_for_LDAK <- sorghum_GWAS_pheno_file_full %>%
  mutate(across(.cols = -GENOTYPE, .fns = ~ . * 10000)) %>%
  mutate(FID = 0,
         IID = GENOTYPE) %>%
  relocate(FID, IID) %>%
  select(!GENOTYPE) %>%
  rename_with(.fn = ~str_replace_all(.x, " ", ".")) #Spaces in colnames break LDAK because columns are already space delineated

write.table(sorghum_file_for_LDAK,
            "data/LDAK/sorghum_LDAK_phenotypes.txt",
            sep = " ",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            na = "NA")


##########################     Maize     ##############################


maize_GWAS_pheno_file_full <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2020_maize_COMBINED_ALL_LEVELS.csv")

maize_file_for_LDAK <- maize_GWAS_pheno_file_full %>%
  mutate(across(.cols = -GENOTYPE, .fns = ~ . * 10000)) %>%
  mutate(FID = 0,
         IID = GENOTYPE) %>%
  relocate(FID, IID) %>%
  select(!GENOTYPE) %>%
  rename_with(.fn = ~str_replace_all(.x, " ", ".")) #Spaces in colnames break LDAK because columns are already space delineated

write.table(maize_file_for_LDAK,
            "data/LDAK/maize_LDAK_phenotypes.txt",
            sep = " ",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            na = "NA")




##########################     Soybean     ##############################


soybean_GWAS_pheno_file_full <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_soybean_COMBINED_ALL_LEVELS.csv")

soybean_file_for_LDAK <- soybean_GWAS_pheno_file_full %>%
  mutate(across(.cols = -GENOTYPE, .fns = ~ . * 10000)) %>%
  mutate(FID = 0,
         IID = GENOTYPE) %>%
  relocate(FID, IID) %>%
  select(!GENOTYPE) %>%
  rename_with(.fn = ~str_replace_all(.x, " ", ".")) #Spaces in colnames break LDAK because columns are already space delineated

write.table(soybean_file_for_LDAK,
            "data/LDAK/soybean_LDAK_phenotypes.txt",
            sep = " ",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            na = "NA")






