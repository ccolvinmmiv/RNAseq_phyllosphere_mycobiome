library(tidyverse)


setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")




# Taxonomic levels
TL_map <- c(
  C = "Class",
  O = "Order",
  `F` = "Family",
  G = "Genus",
  S = "Species"
)

TL_list <- c("C", "O", "F", "G", "S")

# Function to process one taxonomic level file
process_abundance_file <- function(file_path, tl_code, host_species) {
  
  df <- read_csv(file_path)
  
  # Convert wide sample matrix to taxa rows
  df_long <- df %>%
    pivot_longer(
      cols = -GENOTYPE,
      names_to = "Taxon",
      values_to = "Relative_abundance"
    )
  
  # Rebuild wide matrix with taxa as rows
  df_wide <- df_long %>%
    pivot_wider(
      names_from = GENOTYPE,
      values_from = Relative_abundance
    )
  
  # Add metadata columns
  df_wide <- df_wide %>%
    mutate(
      taxon_id = str_extract(Taxon, "[0-9]+$"),
      Taxon = str_remove(Taxon, "_[0-9]+$"),
      Taxonomic_level = TL_map[[tl_code]],
      Host_species = host_species
    ) %>%
    relocate(
      Host_species,
      Taxonomic_level,
      Taxon,
      taxon_id
    )
  
  return(df_wide)
}







# Maize
maize_combined <- map_dfr(TL_list, function(tl) {
  
  process_abundance_file(
    file_path = paste0(
      "data/bracken_intermediates/full_raw_pheno_file_ne_2020_maize_",
      tl,
      ".csv"
    ),
    tl_code = tl,
    host_species = "Maize"
  )
  
})

write_csv(
  maize_combined,
  "data/Supplementary_Dataset_S1_Maize_fungal_relative_abundance_matrix.csv"
)



# Sorghum
sorghum_combined <- map_dfr(TL_list, function(tl) {
  
  process_abundance_file(
    file_path = paste0(
      "data/bracken_intermediates/full_raw_pheno_file_ne_2021_sorgh_",
      tl,
      ".csv"
    ),
    tl_code = tl,
    host_species = "Sorghum"
  )
  
})

write_csv(
  sorghum_combined,
  "data/Supplementary_Dataset_S2_Sorghum_fungal_relative_abundance_matrix.csv"
)



# Sorghum
soybean_combined <- map_dfr(TL_list, function(tl) {
  
  process_abundance_file(
    file_path = paste0(
      "data/bracken_intermediates/full_raw_pheno_file_soybean_",
      tl,
      ".csv"
    ),
    tl_code = tl,
    host_species = "Soybean"
  )
  
})

write_csv(
  soybean_combined,
  "data/Supplementary_Dataset_S3_Soybean_fungal_relative_abundance_matrix.csv"
)









