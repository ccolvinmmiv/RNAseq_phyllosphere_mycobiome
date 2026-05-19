#Processing raw bracken outputs to prepare for co-occurance networks, GWAS/TWAS, etc.
#Requires cleaning genotype names, removing non-fungal taxa, winsorization for data normalization, etc.

library(tidyverse)

setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")

#Taxonomic levels to process (Class, Order, Family, Genus, Species)
TL_list <- c("C", "O", "F", "G", "S")
TL_list <- c("C")


#Read full NCBI taxID conversion table and pull fungal/insect IDs
tax_ID_list <- read_csv("data/raw_bracken_outputs/full_NCBI_taxID_conversion_table.csv")
fungal_only_tax_ID_table <- tax_ID_list %>%
  filter(kingdom == "Fungi")
fungal_only_tax_ID_list <- fungal_only_tax_ID_table %>%
  pull(taxid)


##################################################################################

#Original sorghum RNAseq genotype file from Mangal et al. (2025)
#Convert VCF format to long format with SNP info for each genotype
mangal_2025_genos_vcf <- read_table("genotype_files/SbDiv_RNAseq_GeneticMarkers_Mangal2025.vcf", skip = 23, col_names = TRUE) %>%
  pivot_longer(
    cols = 10:ncol(.),
    names_to = "GENOTYPE",
    values_to = "GENOTYPE_CODE"
  ) %>%
  mutate(
    snp_type = case_when(
      GENOTYPE_CODE == "0|0" ~ REF,
      GENOTYPE_CODE == "1|1" ~ ALT,
      GENOTYPE_CODE %in% c("0|1", "1|0") ~ paste0(REF, ALT),
      TRUE ~ NA_character_
    )
  ) %>%
  select(ID, GENOTYPE, GENOTYPE_CODE, snp_type)

#List of all genotypes from Mangal et al. (2025) for filtering
mangal_genos_list <- mangal_2025_genos_vcf %>% #This has 738 genotypes
  distinct(GENOTYPE) %>%
  pull(GENOTYPE)



for(TL in TL_list) 
{
  #Starting by getting whitelist of genotypes with atleast 2500 reads at class level
  if(TL == "C") { 
    
    print("Starting by getting whitelist of samples with atleast 2500 reads at class level")
    
    #Read class-level Bracken output and clean genotype names
    combined_bracken_X <- read_tsv(str_c("data/raw_bracken_outputs/combined_ne_sorgh_2021_bracken_C.tsv")) %>%
      rename_with(.fn = ~str_remove(.x, str_c(".bracken.", TL)))
     
    #Keep only fungal taxa
    filtered_combined_bracken_X <- filter(combined_bracken_X, taxonomy_id %in% fungal_only_tax_ID_list)
  
    #Select raw abundance columns and clean column names
    pre_otu_table_X <- select(filtered_combined_bracken_X, name, contains("_num")) %>%
      rename_with(.fn = ~str_remove(.x, "_num"))
  
 
    pre_otu_table_2_X <- pre_otu_table_X %>%
      rename("NAME" = name)
    pre_otu_table_2_X_names <- colnames(pre_otu_table_2_X) 
    colnames(pre_otu_table_2_X) <- pre_otu_table_2_X_names
  
    #Pivot to long and back to wide to have GENOTYPE as columns
    longer_pre_otu_table_2_X <- pivot_longer(pre_otu_table_2_X, !NAME, names_to = "GENOTYPE", values_to = "READS" )
    pre_otu_table_3_X <- pivot_wider(longer_pre_otu_table_2_X, id_cols = NAME, values_from = READS, names_from = GENOTYPE) 
    
    #Compute total reads per genotype and filter genotypes passing 2500 read threshold
    passing_genotypes_at_class_level <- pre_otu_table_3_X %>%
      pivot_longer(!NAME, names_to = "GENOTYPE", values_to = "NUMBER_READS")%>%
      group_by(GENOTYPE)%>% 
      mutate(TOTAL_READS = sum(NUMBER_READS, na.rm = TRUE)) %>% 
      filter(TOTAL_READS > 2500) %>% # Require minimum 2500 fungal reads at this level (Class) to keep
      pull(GENOTYPE) 
    
    passing_genotypes_at_class_level <- unique(passing_genotypes_at_class_level)
  
    
    
    
    # For repeated genotypes, finding most "average" sample (which will be kept for GWAS/TWAS)
    expr_tpm <- read_csv("expression_data/raw_expression_data/ne_2021_sorgh_merged_gene_tpms.csv") %>%
      select(!TranscriptID)
    expr_tpm <- expr_tpm[rowSums(expr_tpm) > 0, ]
    # Build sample map (plotID + genotype)
    sample_map <- tibble(SAMPLE = colnames(expr_tpm)) %>%
      mutate(
        plotID   = str_extract(SAMPLE, "^4\\d{3}"),
        GENOTYPE = str_remove(SAMPLE, "^4\\d{3}_")
      ) 
    
    # Log-transform TPMs 
    expr_log <- log2(expr_tpm + 1)
  
    # PCA across samples
    expr_pca <- prcomp(t(expr_log), scale. = TRUE)
    
    pcs <- as_tibble(expr_pca$x[, 1:5], rownames = "SAMPLE") %>%
      left_join(sample_map, by = "SAMPLE")
    
    # Select genotype-representative ("centroid") sample
    rep_samples <- pcs %>%
      group_by(GENOTYPE) %>%
      mutate(
        dist_to_centroid = sqrt(rowSums(
          (across(starts_with("PC")) -
             colMeans(across(starts_with("PC")), na.rm = TRUE))^2
        ))
      ) %>%
      slice_min(dist_to_centroid, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(SAMPLE, plotID, GENOTYPE)
    
    # Freeze for reuse 
    write_csv(rep_samples, "genotype_files/ne_sorgh_2021_representative_samples_by_TPMs.csv")
    
    }
  
  
  #Now actually processing files
  
  
  print(TL)
  
  #Read Bracken outputs and keep only fungal taxa
  combined_bracken_X <- read_tsv(str_c("data/raw_bracken_outputs/combined_ne_sorgh_2021_bracken_", TL, ".tsv")) %>%
    rename_with(.fn = ~str_remove(.x, str_c(".bracken.", TL)))
  filtered_combined_bracken_X <- filter(combined_bracken_X, taxonomy_id %in% fungal_only_tax_ID_list)
  
  #Prepare OTU table
  pre_otu_table_X <- select(filtered_combined_bracken_X, name, contains("_num")) %>%
    rename_with(.fn = ~str_remove(.x, "_num"))
  
  
  pre_otu_table_2_X <- pre_otu_table_X %>%
    rename("NAME" = name)
  pre_otu_table_2_X_names <- colnames(pre_otu_table_2_X) 
  colnames(pre_otu_table_2_X) <- pre_otu_table_2_X_names
  
  longer_pre_otu_table_2_X <- pivot_longer(pre_otu_table_2_X, !NAME, names_to = "GENOTYPE", values_to = "READS" )
  
  
  pre_otu_table_3_X <- pivot_wider(longer_pre_otu_table_2_X, id_cols = NAME, values_from = READS, names_from = GENOTYPE) 
  
  write_csv(pre_otu_table_3_X, str_c("data/bracken_intermediates/pre_otu_table_3_ne_2021_sorgh_", TL, ".csv"), col_names = TRUE)
  pre_otu_table_3_X <- read_csv(str_c("data/bracken_intermediates/pre_otu_table_3_ne_2021_sorgh_", TL, ".csv"))
  
  #Compute relative abundance, filter to class-level passing genotypes, and mark taxa with reads >24
  total_and_rel_abund_table_X <- pre_otu_table_3_X %>%
    pivot_longer(!NAME, names_to = "GENOTYPE", values_to = "NUMBER_READS")%>%
    group_by(GENOTYPE)%>% 
    mutate(TOTAL_READS = sum(NUMBER_READS, na.rm = TRUE)) %>% 
    filter(GENOTYPE %in% passing_genotypes_at_class_level) %>% # All genotypes required minimum 2500 fungal reads at class level to keep
    rowwise() %>%
    mutate(REL_ABUNDANCE = NUMBER_READS/TOTAL_READS) %>% 
    mutate(reads_gt_24 = ifelse(NUMBER_READS > 24, 1, 0))

  #Filter taxa present in at least 1/3 of genotypes
  gt_24_taxa_sums_X <- total_and_rel_abund_table_X %>%
    group_by(NAME) %>%
    summarise(sum_gt_24 = sum(reads_gt_24)) %>%
    ungroup()
  
  number_geno_samples <- length(unique(total_and_rel_abund_table_X$GENOTYPE))

  filtered_taxa_X <- gt_24_taxa_sums_X %>%
    filter(sum_gt_24 >= (number_geno_samples/3)) %>%   
    pull(NAME)
  
  filtered_total_and_rel_abund_table_X <- total_and_rel_abund_table_X %>%
    filter(NAME %in% filtered_taxa_X)
  
  write_csv(filtered_total_and_rel_abund_table_X, str_c("data/bracken_intermediates/filtered_total_and_rel_abund_table_NE_2021_sorgh_", TL, ".csv"), col_names = TRUE)
  
  #Pivot back to wide for GWAS / phenotype files
  filtered_wider_rel_abund_table_X <- select(filtered_total_and_rel_abund_table_X, NAME, GENOTYPE, REL_ABUNDANCE) %>%
    pivot_wider(names_from = NAME, values_from = REL_ABUNDANCE) 
  write_csv(filtered_wider_rel_abund_table_X, str_c("data/bracken_intermediates/filtered_wider_rel_abund_table_NE_2021_sorgh_", TL, ".csv"), col_names = TRUE)
  

  pre_pheno_X <- read_csv(str_c("data/bracken_intermediates/filtered_wider_rel_abund_table_ne_2021_sorgh_", TL, ".csv")) %>%
    rename_with(.fn = ~str_replace(.x, "-", "_")) %>%
    rename_with(.fn = ~str_replace(.x, "-", "_")) 
  
  write_csv(pre_pheno_X, str_c("data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_", TL, ".csv"), col_names = TRUE)
  
  sample_X <- read_csv(str_c("data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_", TL, ".csv"))
  
  #Winsorize abundance values to reduce effects of extreme outliers
  sample_X_winsor <- sample_X %>%
    mutate(across(.cols = -GENOTYPE, .fns = winsorize))
  
  write_csv(sample_X_winsor, str_c("data/bracken_intermediates/full_winsor_pheno_file_ne_2021_sorgh_", TL, ".csv"), col_names = TRUE)

  most_rep_samples_to_keep <- read_csv("genotype_files/ne_sorgh_2021_representative_samples_by_TPMs.csv") %>%
    pull(SAMPLE)
  
  gwas_pheno <- sample_X_winsor %>%
    filter(GENOTYPE %in% most_rep_samples_to_keep) %>%
    mutate(GENOTYPE = str_remove(GENOTYPE, "^4\\d{3}_")) %>%
    filter(GENOTYPE %in% mangal_genos_list) %>%
    #group_by(GENOTYPE) %>%
    #summarise(across(everything(), ~mean(.x, na.rm = FALSE)), .groups = "drop") %>%
    mutate(across(-GENOTYPE, ~na_if(.x, 0)))

  write_csv(gwas_pheno, str_c("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_", TL, ".csv"), col_names = TRUE)
  
}







# Combining raw files (All samples that meet depth threshold at class level, no normalization, all 0's still 0's) (For SpiecEasi or graphing)
full_raw_phenos_C <- read_csv("data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_C.csv")
full_raw_phenos_O <- read_csv("data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_O.csv")
full_raw_phenos_F <- read_csv("data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_F.csv")
full_raw_phenos_G <- read_csv("data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_G.csv")
full_raw_phenos_S <- read_csv("data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_S.csv")


combined_all_levels_raw_phenos <- full_raw_phenos_C %>%
  left_join(full_raw_phenos_O, join_by(GENOTYPE)) %>%
  left_join(full_raw_phenos_F, join_by(GENOTYPE)) %>%
  left_join(full_raw_phenos_G, join_by(GENOTYPE)) %>%
  left_join(full_raw_phenos_S, join_by(GENOTYPE)) 
write_csv(combined_all_levels_raw_phenos, 
          "data/non-normalized_pheno_files/full_ne_2021_sorgh_ALL_LEVELS_raw_phenos.csv", 
          col_names = TRUE)


  
# Combining GWAS files (one sample per genotype, selected most representitative expression sample across genotypes that had multiple replicates, 0's > NA's)
GWAS_phenos_C <- read_csv("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_C.csv")
GWAS_phenos_O <- read_csv("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_O.csv")
GWAS_phenos_F <- read_csv("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_F.csv")
GWAS_phenos_G <- read_csv("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_G.csv")
GWAS_phenos_S <- read_csv("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_S.csv")


combined_all_levels_GWAS_phenos <- GWAS_phenos_C %>%
  left_join(GWAS_phenos_O, join_by(GENOTYPE)) %>%
  left_join(GWAS_phenos_F, join_by(GENOTYPE)) %>%
  left_join(GWAS_phenos_G, join_by(GENOTYPE)) %>%
  left_join(GWAS_phenos_S, join_by(GENOTYPE)) 

write_csv(combined_all_levels_GWAS_phenos, 
          "data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv", 
          col_names = TRUE)
final_gwas_genos_list <- combined_all_levels_GWAS_phenos %>%
  pull(GENOTYPE)
writeLines(final_gwas_genos_list, "genotype_files/sorgh_final_filtered_736_GWAS_genotypes_list.txt")

combined_all_levels_GWAS_phenos <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv")
subset_data_frame_function(combined_all_levels_GWAS_phenos, 2:length(colnames(combined_all_levels_GWAS_phenos)), "data/processed_normalized_phenotype_files_for_GWAS/subsets/sorghum_subsets/sorghum_fungal_winsor_GWAS_subset_", 200)

#For TWAS (Subsets of 50)
combined_all_levels_GWAS_phenos <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv")
subset_data_frame_function(combined_all_levels_GWAS_phenos, 2:length(colnames(combined_all_levels_GWAS_phenos)), "data/TWAS_subsets/sorghum/sorghum_fungal_winsor_TWAS_subset_", 50)















#Custom winsorization function to help deal with non-normal data distribution without dropping samples
#Custom winsorize function (at 0.01 and 0.95)
winsorize <- function(x, lower = 0.01, upper = 0.95) {
  q <- quantile(x, probs = c(lower, upper), na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}


# data: Data frame to split into separate files by col
# cols: cols to split into separate files; tidy-select
# out: path to write file subsets to as csv files
# subsetSize: number of columns of cols per out file
subset_data_frame_function <- function(data, cols, out, subsetSize)
{
  metadata <- dplyr::select(data, !cols)
  data_split <- dplyr::select(data, cols)
  cols_to_split <- ncol(data_split)
  n_subsets <- ceiling(cols_to_split/subsetSize)
  
  for(i in 1:n_subsets)
  {
    outfile <- paste0(out, i, '.csv')
    end_index <- i*subsetSize
    start_index <- end_index - (subsetSize - 1)
    if(end_index > cols_to_split){end_index <- cols_to_split}
    subset <- data_split[, start_index:end_index]
    subset <- bind_cols(metadata, subset)
    write.csv(subset, outfile, quote = FALSE, row.names = FALSE)
  }
}



