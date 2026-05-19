source("http://zzlab.net/GAPIT/gapit_functions.txt") #install GAPIT

#install.packages("data.table") # if not installed, then install before running
library("data.table") # load library
library(tidyverse)



#Defining arguments from slurm script
args <- commandArgs(trailingOnly = FALSE)
pheno_file <- str_remove(args[length(args)-2], fixed('-'))
counts_file <- str_remove(args[length(args)-1], fixed('-'))
gene_location_file <- str_remove(args[length(args)], fixed('-'))

phe <- read.csv(pheno_file)
counts <- read.csv(counts_file)
myGM <- read.csv(gene_location_file)


# Starting TWAS file prep
tpm_data <- counts
transcript_ids <- colnames(tpm_data)[-1]
gene_ids <- sub("\\.\\d+$", "", transcript_ids)

colnames(tpm_data)[-1] <- gene_ids

tpm_gene_level <- tpm_data %>%
  pivot_longer(-GENOTYPE, names_to = "Gene", values_to = "TPM") %>%
  group_by(GENOTYPE, Gene) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  pivot_wider(names_from = Gene, values_from = TPM)

filtered_counts <- counts %>%
  pivot_longer(-GENOTYPE, names_to = "gene", values_to = "TPM") %>%
  group_by(gene) %>%
  filter(mean(TPM < 0.1, na.rm = TRUE) <= 0.5) %>%
  pivot_wider(names_from = gene, values_from = TPM) %>%
  #column_to_rownames("GENOTYPE") %>%
  rename_with(~ str_remove(.x, "\\.\\d+$")) %>%
  dplyr::select(sort(names(.))) 
genes_to_keep <- colnames(filtered_counts)

# Load/compute gene info (position)
myGM <- myGM %>%
  dplyr::filter(gene %in% genes_to_keep) %>%
  arrange(gene) %>%
  mutate(chr = as.numeric(chr)) %>%
  filter(!is.na(chr))

chr_genes_only <- myGM$gene

final_counts_df <- filtered_counts %>%
  dplyr::select(GENOTYPE, dplyr::any_of(chr_genes_only)) %>%
  dplyr::select(sort(names(.))) 

#use quantile method to handle outliers and transform data to 0-2
Quantile_df <- final_counts_df %>%
  mutate(across(-1, ~ {
    min_x <- quantile(.x, 0.05, na.rm = TRUE)
    max_x <- quantile(.x, 0.95, na.rm = TRUE)
    out <- 2 * (.x - min_x) / (max_x - min_x)
    out[out > 2] <- 2
    out[out < 0] <- 0
    out
  }))


Quantile.t <- as.data.frame(Quantile_df[ , -1])

myGD <- Quantile_df %>%
  dplyr::rename(taxa = GENOTYPE) 

final_GD_rownames <- myGD$taxa

rownames(myGD) <- final_GD_rownames
######


# Actually running TWAS
# Loop through phenotype columns (skip the first column (genotype))
for (trait_col in 2:ncol(phe)) {
  trait_name <- colnames(phe)[trait_col]
  cat("Running TWAS for trait:", trait_name, "\n")
  
  myY <- phe[, c(1, trait_col)]  # first column is ID, second is current trait
  rownames(myY) <- myY$GENOTYPE
  
  # Run GAPIT TWAS
  myGAPIT <- GAPIT(
    Y = as.data.frame(myY),
    GD = as.data.frame(myGD),
    GM = as.data.frame(myGM),
    PCA.total = 3,
    model = "CMLM",
    SNP.MAF = 0,
    file.output = FALSE
  )
  
  # Collect TWAS results
  values <- as.data.frame(myGAPIT$GWAS)
  values$FDR <- p.adjust(values$P.value, method = "BH")
  
  # Save to CSV
  out_file <- paste0("TWAS.CMLM_", trait_name, ".csv")
  write_csv(values, out_file)
}

