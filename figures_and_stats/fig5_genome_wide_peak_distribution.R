library(GenomicRanges)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggnewscale)
library(viridis)
library(Hmisc) 
library(ineq)

setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")

sorgh_eQTL_peaks_df <- read_csv("eQTL/all_sorgh_736_genos_eQTL_peaks.csv") %>%
  mutate(CHROM = str_split(top_SNP, "_", simplify = TRUE)[, 1]) %>%
  mutate(CHROM = as.numeric(str_remove(CHROM, "Chr"))) %>%
  mutate(peak_type = "eQTL")
sorgh_GWAS_peaks_df <- read_csv("processed_outputs/GWAS/all_sorgh_736_genos_peaks.csv") %>%
  mutate(CHROM = str_split(top_SNP, "_", simplify = TRUE)[, 1]) %>%
  mutate(CHROM = as.numeric(str_remove(CHROM, "Chr"))) %>%
  mutate(peak_type = "GWAS")

sorgh_combined_peaks <- bind_rows(sorgh_eQTL_peaks_df, sorgh_GWAS_peaks_df)

maize_eQTL_peaks_df <- read_csv("eQTL/all_maize_688_genos_eQTL_peaks.csv") %>%
  mutate(CHROM = str_split(top_SNP, "_", simplify = TRUE)[, 1]) %>%
  mutate(CHROM = as.numeric(str_remove(CHROM, "chr"))) %>%
  mutate(peak_type = "eQTL")
maize_GWAS_peaks_df <- read_csv("processed_outputs/GWAS/all_maize_688_genos_peaks.csv") %>%
  mutate(CHROM = str_split(top_SNP, "_", simplify = TRUE)[, 1]) %>%
  mutate(CHROM = as.numeric(str_remove(CHROM, "chr"))) %>%
  mutate(peak_type = "GWAS")

maize_combined_peaks <- bind_rows(maize_eQTL_peaks_df, maize_GWAS_peaks_df)

##################################################


maize_combined_plot_df <- return_Hotspot_df_Continuous(maize_combined_peaks, species = "maize")
sorghum_combined_plot_df <- return_Hotspot_df_Continuous(sorgh_combined_peaks, species = "sorghum")



######################################################
#print supplemental tables
all_maize_peaks_print_format <- maize_combined_peaks %>%
  mutate(species = "maize")
all_sorgh_peaks_print_format <- sorgh_combined_peaks %>%
  mutate(species = "sorghum")

combined_all_peaks_print_format <- bind_rows(all_maize_peaks_print_format, all_sorgh_peaks_print_format) %>%
  select(species, peak_type, trait, CHROM, POS, top_SNP, top_Pvalue, pStart, pStop, pLength, num_SNPs) %>%
  mutate(trait = str_remove(trait, "called_peaks_"),
         trait = str_remove(trait, ".MLM.csv"))

eQTL_peaks_print <- combined_all_peaks_print_format %>%
  filter(peak_type == "eQTL")
GWAS_peaks_print <- combined_all_peaks_print_format %>%
  filter(peak_type == "GWAS")

write_csv(eQTL_peaks_print, "processed_outputs/eQTL/combined_eQTL_peaks_supplemental.csv")
write_csv(GWAS_peaks_print, "processed_outputs/GWAS/combined_GWAS_peaks_supplemental.csv")

# Chr04: 57.99-59.96 Mb peaks (GWAS + eQTL hotspot)
sorgh_chr_4_hotspot_peaks <- combined_all_peaks_print_format %>%
  filter(species == "sorghum", CHROM == "4", pStart >= 57988801, pStop <= 59955000)
write_csv(sorgh_chr_4_hotspot_peaks, "processed_outputs/sorgh_chr4_hotspot_GWAS_eQTL_peaks_supplemental.csv")

# Chr09: 61.62-63.28 Mb peaks (GWAS hotspot)
sorgh_chr_9_hotspot_peaks <- combined_all_peaks_print_format %>%
  filter(species == "sorghum", CHROM == "9", pStart >= 61620801, pStop <= 63277600)
write_csv(sorgh_chr_9_hotspot_peaks, "processed_outputs/sorgh_chr9_hotspot_GWAS_eQTL_peaks_supplemental.csv")

# Chr10: 58.88-60.74 Mb peaks (eQTL hotspot)
sorgh_chr_10_hotspot_peaks <- combined_all_peaks_print_format %>%
  filter(species == "sorghum", CHROM == "10", pStart >= 58884801, pStop <= 60744700)
write_csv(sorgh_chr_10_hotspot_peaks, "processed_outputs/sorgh_chr10_hotspot_GWAS_eQTL_peaks_supplemental.csv")

######################### Test overlaps ###########################





maize_res <- run_proximity_enrichment(maize_GWAS_peaks_df, maize_eQTL_peaks_df,
                                      species = "maize")
#10k permutations
#============================
#  Proximity enrichment: maize 
#Observed overlaps: 63 
#Mean permuted: 9.8563 
#Fold enrichment: 6.39 
#Permutation p-value: 9.999e-05 
#============================


sorgh_res <- run_proximity_enrichment(sorgh_GWAS_peaks_df, sorgh_eQTL_peaks_df,
                                      species = "sorghum")
#10k permutations
#============================
#  Proximity enrichment: sorghum 
#Observed overlaps: 153 
#Mean permuted: 46.3499 
#Fold enrichment: 3.3 
#Permutation p-value: 9.999e-05 
#============================
###############################################################################################################





#### Fig 5 Panel A 


# Expect fungal_h2_df to exist with:
# taxon | species | h2 | mean_rel_abundance

sorghum_h2_LDAK_out <- read_csv("processed_outputs/LDAK/sorghum_fungal_h2_summary.csv") %>%
  mutate(Species = "Sorghum") %>%
  mutate(Trait = str_replace_all(Trait, "\\.", " "))

sorghum_rel_abundances <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv") %>%
  select(-GENOTYPE)

# Clean column names to just the NCBI taxid (after the last underscore)
sorghum_taxon_ids <- names(sorghum_rel_abundances) %>%
  sapply(function(x) str_split(x, "_")[[1]] %>% tail(1)) %>%
  unname()

colnames(sorghum_rel_abundances) <- sorghum_taxon_ids

# Calculate mean relative abundance per taxon (ignoring NAs)
sorghum_mean_rel_abundance <- colMeans(sorghum_rel_abundances, na.rm = TRUE)

# Create dataframe for merging with h2 data
sorghum_mean_ab_df <- tibble(
  taxon_id = names(sorghum_mean_rel_abundance),
  mean_rel_abundance = sorghum_mean_rel_abundance
)

#
sorghum_fungal_h2_df <- sorghum_h2_LDAK_out %>%
  mutate(taxon_id = str_split(Trait, "_") %>% sapply(tail, 1)) %>% # get taxid for matching
  left_join(sorghum_mean_ab_df, by = "taxon_id") %>%
  filter(Converged == TRUE)

##############

maize_rel_abundances <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2020_maize_COMBINED_ALL_LEVELS.csv") %>%
  select(-GENOTYPE)
  
maize_h2_LDAK_out <- read_csv("processed_outputs/LDAK/maize_fungal_h2_summary.csv") %>%
  mutate(Species = "Maize") %>%
  mutate(Trait = str_replace_all(Trait, "\\.", " "))

# Clean column names to just the NCBI taxid (after the last underscore)
maize_taxon_ids <- names(maize_rel_abundances) %>%
  sapply(function(x) str_split(x, "_")[[1]] %>% tail(1)) %>%
  unname()

colnames(maize_rel_abundances) <- maize_taxon_ids

# Calculate mean relative abundance per taxon (ignoring NAs)
maize_mean_rel_abundance <- colMeans(maize_rel_abundances, na.rm = TRUE)

# Create dataframe for merging with h2 data
maize_mean_ab_df <- tibble(
  taxon_id = names(maize_mean_rel_abundance),
  mean_rel_abundance = maize_mean_rel_abundance
)

# Example merge with fungal_h2_df
maize_fungal_h2_df <- maize_h2_LDAK_out %>%
  mutate(taxon_id = str_split(Trait, "_") %>% sapply(tail, 1)) %>% # get taxid for matching
  left_join(maize_mean_ab_df, by = "taxon_id") %>%
  filter(Converged == TRUE)


combined_h2_rel_abundance <- bind_rows(maize_fungal_h2_df, sorghum_fungal_h2_df) 

h2_final <- combined_h2_rel_abundance %>%
  select(
    Species,
    Trait,
    taxon_id,
    h2,
    h2_SE,
    N,
    mean_rel_abundance
  )

write_csv(h2_final, "processed_outputs/LDAK/combined_h2_dist_table.csv")

## Checking h2 correlations across species
shared_taxa <- intersect(maize_fungal_h2_df$taxon_id,
                         sorghum_fungal_h2_df$taxon_id)

maize_shared <- maize_fungal_h2_df %>%
  filter(taxon_id %in% shared_taxa) %>%
  select(taxon_id, h2) %>%
  mutate(h2_maize = h2) %>%
  select(taxon_id, h2_maize) 

sorgh_shared <- sorghum_fungal_h2_df %>%
  filter(taxon_id %in% shared_taxa) %>%
  select(taxon_id, h2) %>%
  mutate(h2_sorghum = h2) %>%
  select(taxon_id, h2_sorghum) 

shared_h2_df <- inner_join(maize_shared, sorgh_shared, by = "taxon_id")

# Add GWAS presence

maize_peak_binary <- maize_gwas_counts %>%
  mutate(has_peak_maize = n > 0) %>%
  select(trait_ID, has_peak_maize)

sorgh_peak_binary <- sorgh_gwas_counts %>%
  mutate(has_peak_sorghum = n > 0) %>%
  select(trait_ID, has_peak_sorghum)

shared_h2_df <- shared_h2_df %>%
  left_join(maize_peak_binary, by = c("taxon_id" = "trait_ID")) %>%
  left_join(sorgh_peak_binary, by = c("taxon_id" = "trait_ID")) %>%
  mutate(
    has_peak_maize = replace_na(has_peak_maize, FALSE),
    has_peak_sorghum = replace_na(has_peak_sorghum, FALSE)
  )
# Correlation of h2 across species

cor_test <- cor.test(shared_h2_df$h2_maize,
                     shared_h2_df$h2_sorghum,
                     method = "spearman")

cor_test

# Define heritable threshold 
h2_thresh <- 0.1

shared_h2_df <- shared_h2_df %>%
  mutate(
    heritable_maize = h2_maize > h2_thresh,
    heritable_sorghum = h2_sorghum > h2_thresh
  )

# Contingency table

heritability_table <- table(shared_h2_df$heritable_maize,
                            shared_h2_df$heritable_sorghum)

heritability_table

fisher_herit <- fisher.test(heritability_table)

fisher_herit

# GWAS overlap consistency

gwas_table <- table(shared_h2_df$has_peak_maize,
                    shared_h2_df$has_peak_sorghum)

gwas_table

fisher_gwas <- fisher.test(gwas_table)

fisher_gwas


#####################
#Plot panel A

# First, calculate binned h2 and proportion per species
combined_h2_rel_abundance <- combined_h2_rel_abundance %>%
  filter(!is.na(h2) & !is.na(mean_rel_abundance)) %>%
  group_by(Species) %>%
  mutate(h2_bin = cut(h2, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) %>%
  group_by(Species, h2_bin) %>%
  mutate(bin_count = n()) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(prop_in_bin = bin_count / sum(bin_count)) %>%
  ungroup()

bin_width <- 0.05

binned_df <- combined_h2_rel_abundance %>%
  filter(!is.na(h2), !is.na(mean_rel_abundance)) %>%
  mutate(h2_bin = floor(h2 / bin_width) * bin_width) %>%
  group_by(Species, h2_bin) %>%
  summarise(
    taxa_in_bin = n(),
    mean_abundance_bin = mean(mean_rel_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Species) %>%
  mutate(prop_taxa = taxa_in_bin / sum(taxa_in_bin)) %>%
  ungroup()



methods <- c("spearman", "pearson", "kendall")

#Checking correlation Maize
for (m in methods) {
  res <- cor.test(maize_fungal_h2_df$mean_rel_abundance, maize_fungal_h2_df$h2, method = m, exact = FALSE)
  cat("\nMethod:", m, "\n")
  cat("Correlation:", res$estimate, "\n")
  cat("p-value:", res$p.value, "\n")
}


#Checking sorghum
for (m in methods) {
  res <- cor.test(sorghum_fungal_h2_df$mean_rel_abundance, sorghum_fungal_h2_df$h2, method = m, exact = FALSE)
  cat("\nMethod:", m, "\n")
  cat("Correlation:", res$estimate, "\n")
  cat("p-value:", res$p.value, "\n")
}
#Plotting

max_val <- max(binned_df$mean_abundance_bin, na.rm = TRUE)

panelA_plot_maize <- ggplot(filter(binned_df, Species == "Maize"),
                            aes(h2_bin, prop_taxa, fill = mean_abundance_bin)) +
  geom_col(color = "black", width = bin_width) +
  scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0, 1)) +
  annotate("text", x = 0.12, y = 0.22, label = "Maize", size = 4, fontface = "bold") +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Mean rel. abundance",
    limits = c(0.0, max_val),
    breaks = c(0.01, 0.03, 0.05),
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(4, "cm"),   # smaller legend bar
      barheight = unit(0.25, "cm")
    )
  ) +
  labs(#x = expression("SNP heritability ("*h^2*")"),
    x = "",
       y = "Proportion of taxa") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(size = 9, color = "black"),
    legend.position = c(0.7, 0.7),
    legend.direction = "horizontal",
    legend.title = element_text(size = 9, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    legend.background = element_blank()
  )

panelA_plot_sorghum <- ggplot(filter(binned_df, Species == "Sorghum"),
                              aes(h2_bin, prop_taxa, fill = mean_abundance_bin)) +
  geom_col(color = "black", width = bin_width) +
  scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0, 1)) +
  annotate("text", x = 0.18, y = 0.23, label = "Sorghum", size = 4, fontface = "bold") +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Mean rel. abundance",
    limits = c(0.0, max_val),
    breaks = c(0.01, 0.03, 0.05),
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(2, "cm"),   # smaller legend bar
      barheight = unit(0.25, "cm")
    )
  ) +
  labs(x = expression("SNP heritability ("*h^2*")"),
       y = "Proportion of taxa") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(size = 9, color = "black"),
    #legend.position = c(0.7, 0.7),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_text(size = 9, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    legend.background = element_blank()
  )
panelA_combined <- (panelA_plot_maize / panelA_plot_sorghum) +
  plot_layout(heights = c(1, 1))#, guides = "collect") &
  #theme(legend.position = c(0.5, 1.02))

panelA_combined

########################
#### Fig 5 Panel B #####
########################

total_maize_taxa <- length(unique(maize_taxon_ids))
total_sorghum_taxa <- length(unique(sorghum_taxon_ids))

# Count GWAS peaks per taxon
maize_gwas_counts <- maize_GWAS_peaks_df %>%
  mutate(trait_ID = sapply(trait, function(x) str_split(x, "_")[[1]] %>% tail(1))) %>%
  mutate(trait_ID = str_remove_all(trait_ID, ".MLM.csv")) %>%
  count(trait_ID) %>%
  mutate(species = "Maize")

sorgh_gwas_counts <- sorgh_GWAS_peaks_df %>%
  mutate(trait_ID = sapply(trait, function(x) str_split(x, "_")[[1]] %>% tail(1))) %>%
  mutate(trait_ID = str_remove_all(trait_ID, ".MLM.csv")) %>%
  count(trait_ID) %>%
  mutate(species = "Sorghum")


# Raw GWAS counts for MS
maize_total_peaks <- nrow(maize_GWAS_peaks_df)
sorgh_total_peaks <- nrow(sorgh_GWAS_peaks_df)

maize_taxa_with_peaks <- sum(maize_gwas_counts$n > 0)
sorgh_taxa_with_peaks <- sum(sorgh_gwas_counts$n > 0)

cat("\nGWAS RAW COUNTS\n")
cat("Maize: ", maize_taxa_with_peaks, " taxa with peaks; ", maize_total_peaks, " total peaks\n")
cat("Sorghum: ", sorgh_taxa_with_peaks, " taxa with peaks; ", sorgh_total_peaks, " total peaks\n")

# Add taxa with 0 peaks
maize_zero <- tibble(
  trait_ID = setdiff(unique(maize_fungal_h2_df$taxon_id), maize_gwas_counts$trait_ID),
  n = 0,
  species = "Maize"
)

sorghum_zero <- tibble(
  trait_ID = setdiff(unique(sorghum_fungal_h2_df$taxon_id), sorgh_gwas_counts$trait_ID),
  n = 0,
  species = "Sorghum"
)

# Combine counts including zeros
gwas_counts <- bind_rows(maize_gwas_counts, sorgh_gwas_counts,
                         maize_zero, sorghum_zero) %>%
  mutate(category = case_when(
    n == 0 ~ "0",
    n == 1 ~ "1",
    n >= 2 ~ "2+"
  ))

# Merge h2 with GWAS detectability

maize_detectability_df <- maize_fungal_h2_df %>%
  select(taxon_id, h2) %>%
  left_join(maize_gwas_counts, by = c("taxon_id" = "trait_ID")) %>%
  mutate(n = replace_na(n, 0),
         has_peak = n > 0,
         species = "Maize")

sorgh_detectability_df <- sorghum_fungal_h2_df %>%
  select(taxon_id, h2) %>%
  left_join(sorgh_gwas_counts, by = c("taxon_id" = "trait_ID")) %>%
  mutate(n = replace_na(n, 0),
         has_peak = n > 0,
         species = "Sorghum")

detectability_df <- bind_rows(maize_detectability_df, sorgh_detectability_df)

# Logistic regression: does h2 predict GWAS detection?

detectability_model <- glm(has_peak ~ h2 * species,
                           data = detectability_df,
                           family = binomial)

summary(detectability_model)


# Compute proportions per species
gwas_prop <- gwas_counts %>%
  group_by(species) %>%
  count(category) %>%
  mutate(prop_taxa = n / sum(n))
gwas_prop <- gwas_prop %>%
  mutate(category = factor(category, levels = c("2+", "1", "0")))

# Plot stacked bar
panelB_gwas_stacked <- ggplot(gwas_prop,
                              aes(x = species, y = prop_taxa, fill = category)) +
  geom_col(color = "black") +
  scale_fill_manual(
    values = c(
      "0" = "lightcyan1",  
      "1" = "cornflowerblue",  
      "2+" = "blue"  
    ),
    name = "GWAS peaks\nper taxon"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.5), nrow = 3, title.position = "left")) +
  labs(x = NULL, y = "Proportion of taxa") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 9),
    legend.title = element_text(color = "black", size = 9),
    #legend.text = element_text(color = "black", size = 9),
    legend.position = "top",
    legend.title.position = "left",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9, color = "black", margin = margin(b = 1)),
    legend.key.height = unit(0.4, "lines"),
    legend.key.width = unit(0.4, "lines"),
    legend.background = element_rect(fill = NA, color = NA),
    plot.title = element_blank(),
    legend.box.spacing = unit(1, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.spacing.y = unit(1, "pt")
  )

panelB_gwas_stacked <- panelB_gwas_stacked +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.spacing.x = unit(6, "pt"),       
    legend.box.spacing = unit(4, "pt")
  )

panelB_gwas_stacked


# Test if the proportion of 0/1/2+ differs between species
# Build contingency table from counts (not proportions)
gwas_counts_table <- gwas_counts %>%
  group_by(species, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
  column_to_rownames("species") %>%
  as.matrix()

# Chi-squared test
chi_test <- chisq.test(gwas_counts_table)
chi_test





########################
#### Fig 5 Panel C #####
########################


# Function to flag eQTLs near GWAS peaks
flag_eqtl_near_gwas <- function(gwas_df, eqtl_df, radius = 250000) {
  gwas_gr <- GRanges(seqnames = gwas_df$CHROM,
                     ranges = IRanges(gwas_df$POS, gwas_df$POS))
  eqtl_gr <- GRanges(seqnames = eqtl_df$CHROM,
                     ranges = IRanges(eqtl_df$POS, eqtl_df$POS))
  gwas_expanded <- resize(gwas_gr, width = 2*radius + 1, fix = "center")
  overlaps <- countOverlaps(eqtl_gr, gwas_expanded) > 0
  eqtl_df$near_gwas <- overlaps
  eqtl_df
}

# Flag eQTLs near GWAS peaks for maize and sorghum
maize_eqtl_flagged <- flag_eqtl_near_gwas(maize_GWAS_peaks_df, maize_eQTL_peaks_df) %>%
  mutate(species = "Maize")

sorgh_eqtl_flagged <- flag_eqtl_near_gwas(sorgh_GWAS_peaks_df, sorgh_eQTL_peaks_df) %>%
  mutate(species = "Sorghum")

# Combine species
eqtl_combined <- bind_rows(maize_eqtl_flagged, sorgh_eqtl_flagged)

# Count eQTL categories per species
eqtl_counts <- eqtl_combined %>%
  mutate(trait = str_remove(trait, "called_peaks_"),
                      str_remove(trait, ".MLM.csv")) %>%
  group_by(species, gene_ID = trait) %>%
  summarise(near_gwas = any(near_gwas), .groups = "drop") %>%
  mutate(category = case_when(
    near_gwas ~ "eQTL near GWAS",
    !near_gwas ~ "eQTL not near GWAS"
  ))

sorgh_eQTL_genes_tested <- read_csv("TWAS_results/sig_genes/sorghum_736_TWAS_all_sig_hits.csv") %>%
  pull(SNP)
maize_eQTL_genes_tested <- read_csv("TWAS_results/sig_genes/maize_688_TWAS_all_sig_hits.csv") %>%
  pull(SNP)

# Add genes with no eQTL detected
maize_no_eqtl <- setdiff(maize_eQTL_genes_tested, eqtl_counts$gene_ID[eqtl_counts$species == "Maize"])
sorgh_no_eqtl <- setdiff(sorgh_eQTL_genes_tested, eqtl_counts$gene_ID[eqtl_counts$species == "Sorghum"])

eqtl_zero <- tibble(
  gene_ID = c(maize_no_eqtl, sorgh_no_eqtl),
  category = "No eQTL",
  species = c(rep("Maize", length(maize_no_eqtl)),
              rep("Sorghum", length(sorgh_no_eqtl)))
)

# Combine all categories
eqtl_summary <- bind_rows(eqtl_counts %>% select(gene_ID, category, species),
                          eqtl_zero) %>%
  group_by(species, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(prop_genes = n / sum(n)) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c("eQTL near GWAS", "eQTL not near GWAS", "No eQTL")))

# Plot stacked bar
panelC_eqtl_stacked <- ggplot(eqtl_summary,
                              aes(x = species, y = prop_genes, fill = category)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c(
    "No eQTL" = "#fee0d2",             # light pink
    "eQTL not near GWAS" = "orange",  # medium red
    "eQTL near GWAS" = "darkorange3"       # dark red
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Proportion of genes", fill = NULL) +
  guides(fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.5), nrow = 3, title.position = "top")) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 9),
    legend.title = element_text(color = "black", size = 9),
    legend.position = "top",
    legend.title.position = "",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9, color = "black", margin = margin(b = 1)),
    legend.key.height = unit(0.4, "lines"),
    legend.key.width = unit(0.4, "lines"),
    legend.background = element_rect(fill = NA, color = NA),
    plot.title = element_blank(),
    legend.box.spacing = unit(1, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.spacing.y = unit(1, "pt")
  )
panelC_eqtl_stacked <- panelC_eqtl_stacked +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.spacing.x = unit(6, "pt"),
    legend.box.spacing = unit(4, "pt")
  )
panelC_eqtl_stacked





















##### Plotting   #######
maize_res_plot <- make_enrichment_plot(maize_res)
maize_res_plot
#ggsave("graphs/test_maize_res_plot.png", plot = maize_res_plot, dpi = 300, units = "in", height = 2.5, width = 2.5, bg = "white") 

sorgh_res_plot <- make_enrichment_plot(sorgh_res)
sorgh_res_plot



maize_combined_plot <- plot_MultiType_Hotspot_Continuous(maize_combined_peaks, species = "maize")
sorghum_combined_plot <- plot_MultiType_Hotspot_Continuous(sorgh_combined_peaks, species = "sorghum")


#Supplemental scaled hotspot plots
maize_combined_plot_scaled <- plot_MultiType_Hotspot_scaled(maize_combined_peaks, species = "maize")
sorghum_combined_plot_scaled <- plot_MultiType_Hotspot_scaled(sorgh_combined_peaks, species = "sorghum")

combined_scaled_plot <- maize_combined_plot_scaled / sorghum_combined_plot_scaled
ggsave("graphs/test_scaled_hotspots_plot.png", plot = combined_scaled_plot, dpi = 300, width = 6.5, height = 5, units = "in", bg = "white")



# Combine Plots 
panelA_block <- panelA_plot_maize / panelA_plot_sorghum +
  plot_layout(heights = c(1, 1))

#panelBC_block <- panelB_gwas_stacked | panelC_eqtl_stacked +
#  plot_layout(heights = c(1, 1))

#panelABC_row <- free(panelA_block) + panelBC_block #+
 # plot_layout(widths = c(1.4, 1, heights = c(1, 1)))

panelABC_row <- wrap_plots(
  free(panelA_block),
  plot_spacer(),
  panelB_gwas_stacked,
  panelC_eqtl_stacked,
  nrow = 1,
  widths = c(2, 0.1, 1, 1)
) &
  theme(plot.margin = margin(2, 2, 2, 2)) 

panelABC_row

ggsave("graphs/test_fig5_panelABC_row.png", plot = panelABC_row, dpi = 300, units = "in", height = 3, width = 6.5, bg = "white") 



maize_peaks_and_histo <- maize_combined_plot + maize_res_plot + plot_layout(widths = c(2.25, 1))

sorghum_peaks_and_histo <- sorghum_combined_plot + sorgh_res_plot + plot_layout(widths = c(2.25, 1))



stacked_maize_sorghum_peaks_and_histo <- wrap_plots(panelABC_row, maize_peaks_and_histo, sorghum_peaks_and_histo, nrow = 3, heights = c(1, 0.9, 0.9)) # &
 # plot_annotation(tag_levels = "A") & 
 # theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("graphs/final_GWAS_eQTL_stacked_hotspot_plots.png", plot = stacked_maize_sorghum_peaks_and_histo, dpi = 600, units = "in", height = 9, width = 6.5, bg = "white") 















####################### Stats for MS ####################



# h2 DISTRIBUTION STATS

calc_h2_stats <- function(df, species_name) {
  df_clean <- df %>% filter(!is.na(h2))
  
  tibble(
    species = species_name,
    n_taxa = nrow(df_clean),
    
    median_h2 = median(df_clean$h2),
    mean_h2 = mean(df_clean$h2),
    
    h2_IQR_low = quantile(df_clean$h2, 0.25),
    h2_IQR_high = quantile(df_clean$h2, 0.75),
    
    h2_min = min(df_clean$h2),
    h2_max = max(df_clean$h2),
    
    h2_p05 = quantile(df_clean$h2, 0.05),
    h2_p95 = quantile(df_clean$h2, 0.95),
    
    # % taxa with non-zero heritability 
    prop_h2_gt0 = mean(df_clean$h2 > 0),
    prop_h2_gt01 = mean(df_clean$h2 > 0.1)
  )
}

maize_h2_stats <- calc_h2_stats(maize_fungal_h2_df, "Maize")
sorghum_h2_stats <- calc_h2_stats(sorghum_fungal_h2_df, "Sorghum")

bind_rows(maize_h2_stats, sorghum_h2_stats)




# h2 ~ ABUNDANCE CORRELATION

calc_cor_stats <- function(df, species_name) {
  df_clean <- df %>% filter(!is.na(h2), !is.na(mean_rel_abundance))
  
  spearman <- cor.test(df_clean$h2, df_clean$mean_rel_abundance, method = "spearman")
  
  tibble(
    species = species_name,
    spearman_rho = spearman$estimate,
    spearman_p = spearman$p.value
  )
}

bind_rows(
  calc_cor_stats(maize_fungal_h2_df, "Maize"),
  calc_cor_stats(sorghum_fungal_h2_df, "Sorghum")
)




# GWAS ARCHITECTURE STATS

calc_gwas_architecture <- function(gwas_counts_df, species_name) {
  
  df <- gwas_counts_df %>% filter(species == species_name)
  
  tibble(
    species = species_name,
    
    prop_0 = mean(df$n == 0),
    prop_1 = mean(df$n == 1),
    prop_2plus = mean(df$n >= 2),
    
    mean_peaks_per_taxon = mean(df$n),
    median_peaks_per_taxon = median(df$n),
    
    max_peaks = max(df$n)
  )
}

bind_rows(
  calc_gwas_architecture(gwas_counts, "Maize"),
  calc_gwas_architecture(gwas_counts, "Sorghum")
)




# eQTL COLOCALIZATION STATS

calc_eqtl_stats <- function(eqtl_summary_df, species_name) {
  
  df <- eqtl_summary_df %>% filter(species == species_name)
  
  tibble(
    species = species_name,
    
    prop_eqtl_near = df$prop_genes[df$category == "eQTL near GWAS"],
    prop_eqtl_far = df$prop_genes[df$category == "eQTL not near GWAS"],
    prop_no_eqtl = df$prop_genes[df$category == "No eQTL"]
  )
}

bind_rows(
  calc_eqtl_stats(eqtl_summary, "Maize"),
  calc_eqtl_stats(eqtl_summary, "Sorghum")
)




# ENRICHMENT STATS (already computed, just format cleanly)

format_enrichment <- function(res, species_name) {
  tibble(
    species = species_name,
    observed_overlaps = res$observed,
    mean_permuted = mean(res$perm_hits),
    fold_enrichment = res$fold,
    p_value = res$p_value
  )
}

bind_rows(
  format_enrichment(maize_res, "Maize"),
  format_enrichment(sorgh_res, "Sorghum")
)

#10k permutations
#============================
#  Proximity enrichment: maize 
#Observed overlaps: 63 
#Mean permuted: 9.8563 
#Fold enrichment: 6.39 
#Permutation p-value: 9.999e-05 
#============================

#============================
#  Proximity enrichment: sorghum 
#Observed overlaps: 153 
#Mean permuted: 46.3499 
#Fold enrichment: 3.3 
#Permutation p-value: 9.999e-05 
#============================


# HOTSPOT / CLUSTERING METRICS 


calc_hotspot_stats <- function(peaks_df, species_name) {
  
  # count peaks per 1Mb window (already computed earlier if using hotspot df)
  hotspot_df <- return_Hotspot_df_Continuous(peaks_df, species = tolower(species_name))
  
  hotspot_summary <- hotspot_df %>%
    group_by(peak_type) %>%
    summarise(
      mean_peaks_per_window = mean(peaks_count),
      max_peaks_window = max(peaks_count),
      windows_gt10 = mean(peaks_count >= 10),
      .groups = "drop"
    )
  
  hotspot_summary$species <- species_name
  
  hotspot_summary
}

bind_rows(
  calc_hotspot_stats(maize_combined_peaks, "Maize"),
  calc_hotspot_stats(sorgh_combined_peaks, "Sorghum")
)


# GINI stats
# compute Gini from hotspot df
build_binned_counts <- function(peaks_df, species_name, bin_size = 1e6) {
  
  peaks_df %>%
    mutate(
      bin = floor(POS / bin_size),
      chr = as.character(CHROM),
      species = species_name
    ) %>%
    group_by(species, peak_type, chr, bin) %>%
    summarise(
      peaks_count = n(),
      .groups = "drop"
    )
}

binned_counts <- bind_rows(
  build_binned_counts(maize_combined_peaks, "Maize"),
  build_binned_counts(sorgh_combined_peaks, "Sorghum")
)


calc_gini <- function(df) {
  
  df %>%
    group_by(species, peak_type) %>%
    summarise(
      gini = ineq(peaks_count, type = "Gini"),
      
      prop_peaks_top10 = {
        x <- peaks_count
        x_sorted <- sort(x, decreasing = TRUE)
        sum(x_sorted[1:ceiling(0.10 * length(x_sorted))]) / sum(x_sorted)
      },
      
      prop_peaks_top5 = {
        x <- peaks_count
        x_sorted <- sort(x, decreasing = TRUE)
        sum(x_sorted[1:ceiling(0.05 * length(x_sorted))]) / sum(x_sorted)
      },
      
      total_peaks = sum(peaks_count), 
      n_bins = n(),
      .groups = "drop"
    )
}

calc_gini(binned_counts)






########################      Functions         ##########################




run_proximity_enrichment <- function(gwas_df, eqtl_df,
                                     species = "maize",
                                     radius = 250000,
                                     n_perm = 10000) {
  
  library(GenomicRanges)
  
  
  # Genome definition
  
  genome_df <- if (species == "maize") {
    data.frame(
      chr = as.character(1:10),
      end = c(308452471, 243675191, 238017767, 250330460, 226353449,
              181357234, 185808916, 182411202, 163004744, 152435371)
    )
  } else if (species == "sorghum") {
    data.frame(
      chr = as.character(1:10),
      end = c(85112863, 79114963, 80873341, 71215609, 77058072,
              62713908, 68911884, 65779274, 63277606, 62870657)
    )
  } else {
    stop("Species must be 'maize' or 'sorghum'")
  }
  
  
  # Convert to GRanges
  
  gwas_gr <- GRanges(seqnames = gwas_df$CHROM,
                     ranges = IRanges(gwas_df$POS, gwas_df$POS))
  
  eqtl_gr <- GRanges(seqnames = eqtl_df$CHROM,
                     ranges = IRanges(eqtl_df$POS, eqtl_df$POS))
  
  
  # Observed overlaps
  
  gwas_expanded <- resize(gwas_gr, width = 2*radius + 1, fix = "center")
  observed_hits <- sum(countOverlaps(gwas_expanded, eqtl_gr) > 0)
  
  
  # Permutations
  
  perm_hits <- replicate(n_perm, {
    
    random_pos <- sapply(gwas_df$CHROM, function(chr) {
      chr_len <- genome_df$end[genome_df$chr == chr]
      sample(chr_len, 1)
    })
    
    gwas_rand <- GRanges(seqnames = gwas_df$CHROM,
                         ranges = IRanges(random_pos, random_pos))
    
    gwas_rand_exp <- resize(gwas_rand, width = 2*radius + 1, fix = "center")
    
    sum(countOverlaps(gwas_rand_exp, eqtl_gr) > 0)
  })
  
  
  # Stats
  
  pval <- (sum(perm_hits >= observed_hits) + 1) / (n_perm + 1)
  fold_enrichment <- observed_hits / mean(perm_hits)
  
  
  # Print results
  
  cat("\n============================\n")
  cat("Proximity enrichment:", species, "\n")
  cat("Observed overlaps:", observed_hits, "\n")
  cat("Mean permuted:", mean(perm_hits), "\n")
  cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
  cat("Permutation p-value:", pval, "\n")
  cat("============================\n")
  
  return(list(observed = observed_hits,
              perm_hits = perm_hits,
              p_value = pval,
              fold = fold_enrichment))
}


make_enrichment_plot <- function(res, species_label = NULL) {
  
  ggplot(data.frame(x = res$perm_hits), aes(x = x)) +
    geom_histogram(bins = 40, fill = "grey70", color = "black") +
    geom_vline(xintercept = res$observed,
               color = "red",
               linewidth = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    labs(x = "GWAS peaks near eQTL",
         y = "Permutations",
         title = species_label) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", size = 9),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(),
      axis.text = element_text(color = "black", size = 9),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank())
}


plot_MultiType_Hotspot_Continuous <- function(peaks_df, 
                                              window_size = 1e6, 
                                              step_size = 100, 
                                              main = NULL, 
                                              species = 'maize') {

  # Genome definition
  
  genome <- if (species == 'maize') {
    data.frame(
      chr = as.character(1:10),
      end = c(308452471, 243675191, 238017767, 250330460, 226353449,
              181357234, 185808916, 182411202, 163004744, 152435371),
      start = rep(1, 10)
    )
  } else if (species == 'sorghum') {
    data.frame(
      chr = as.character(1:10),
      end = c(85112863, 79114963, 80873341, 71215609, 77058072,
              62713908, 68911884, 65779274, 63277606, 62870657),
      start = rep(1, 10)
    )
  } else {
    stop("Species must be 'maize' or 'sorghum'")
  }
  
  genome <- genome %>%
    arrange(as.numeric(chr)) %>%
    mutate(cumstart = cumsum(c(0, head(end, -1))),
           midpoint = cumstart + end / 2)
  
  
  # Create sliding windows
  
  window_list <- mapply(make_windows,
                        chr = genome$chr,
                        chr_len = genome$end,
                        window_size = window_size,
                        step_size = step_size,
                        SIMPLIFY = FALSE)
  
  windows <- do.call(c, unname(window_list))
  
  
  # Split by peak_type
  
  peak_list <- split(peaks_df, peaks_df$peak_type)
  
  combined_df <- do.call(rbind, lapply(names(peak_list), function(ptype) {
    
    df <- peak_list[[ptype]]
    
    peaks_gr <- GRanges(seqnames = df$CHROM,
                        ranges = IRanges(start = df$POS, end = df$POS))
    
    counts <- countOverlaps(windows, peaks_gr)
    
    temp <- as.data.frame(windows) %>%
      mutate(
        chr = as.character(seqnames),
        peaks_count = counts,
        peak_type = ptype
      )
    
    temp
  }))
  
  
  # Add cumulative genome position
  
  combined_df <- combined_df %>%
    left_join(genome %>% select(chr, cumstart), by = "chr") %>%
    mutate(genome_pos = start + cumstart,
           genome_pos_Mb = genome_pos / 1e6)
  
  
  # Plot
  
  p1 <- ggplot(combined_df,
               aes(x = genome_pos_Mb,
                   y = peaks_count,
                   color = peak_type)) +
    geom_line(linewidth = 0.7, alpha = 0.7) +
    geom_vline(xintercept = genome$cumstart[-1] / 1e6,
               linetype = "dotted",
               color = "grey50") +
    
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 60, ymax = 80, fill = "white", color = NA) +  #Maize only
    
    scale_color_manual(values = c("GWAS" = "blue",
                                  "eQTL" = "darkorange")) +
    scale_x_continuous(
      name = paste0(toupper(substr(species,1,1)),
                    substr(species,2,nchar(species)),
                    " Chromosome"),
      breaks = genome$midpoint / 1e6,
      labels = genome$chr,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0))) +
    labs(
      y = "Number of Peaks",
      title = main,
      color = NULL
    ) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", size = 9),
      axis.text = element_text(color = "black", size = 9),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = c(0.5, 0.95), legend.direction = "horizontal",    #Maize Only
      #legend.position = "none"   #Sorghum only
    )
  
  return(p1)

  
}

plot_MultiType_Hotspot_scaled <- function(peaks_df, 
                                          window_size = 1e6, 
                                          step_size = 100, 
                                          main = NULL, 
                                          species = 'maize') {
  
  # Genome definition
  genome <- if (species == 'maize') {
    data.frame(
      chr = as.character(1:10),
      end = c(308452471, 243675191, 238017767, 250330460, 226353449,
              181357234, 185808916, 182411202, 163004744, 152435371),
      start = rep(1, 10)
    )
  } else if (species == 'sorghum') {
    data.frame(
      chr = as.character(1:10),
      end = c(85112863, 79114963, 80873341, 71215609, 77058072,
              62713908, 68911884, 65779274, 63277606, 62870657),
      start = rep(1, 10)
    )
  } else stop("Species must be 'maize' or 'sorghum'")
  
  genome <- genome %>%
    arrange(as.numeric(chr)) %>%
    mutate(cumstart = cumsum(c(0, head(end, -1))),
           midpoint = cumstart + end / 2)
  
  # Create sliding windows
  window_list <- mapply(make_windows,
                        chr = genome$chr,
                        chr_len = genome$end,
                        window_size = window_size,
                        step_size = step_size,
                        SIMPLIFY = FALSE)
  windows <- do.call(c, unname(window_list))
  
  # Split by peak_type
  peak_list <- split(peaks_df, peaks_df$peak_type)
  
  combined_df <- do.call(rbind, lapply(names(peak_list), function(ptype) {
    df <- peak_list[[ptype]]
    peaks_gr <- GRanges(seqnames = df$CHROM,
                        ranges = IRanges(start = df$POS, end = df$POS))
    counts <- countOverlaps(windows, peaks_gr)
    
    # Normalize by total peaks per type and optionally sqrt
    prop_counts <- counts / length(peaks_gr)
    prop_counts_sqrt <- sqrt(prop_counts)
    
    temp <- as.data.frame(windows) %>%
      mutate(
        chr = as.character(seqnames),
        peaks_count_scaled = prop_counts_sqrt,
        peak_type = ptype
      )
    temp
  }))
  
  # Add cumulative genome position
  combined_df <- combined_df %>%
    left_join(genome %>% select(chr, cumstart), by = "chr") %>%
    mutate(genome_pos = start + cumstart,
           genome_pos_Mb = genome_pos / 1e6)
  
  # Plot
  p1 <- ggplot(combined_df,
               aes(x = genome_pos_Mb,
                   y = peaks_count_scaled,
                   color = peak_type)) +
    geom_line(linewidth = 0.7, alpha = 0.7) +
    geom_vline(xintercept = genome$cumstart[-1] / 1e6,
               linetype = "dotted",
               color = "grey50") +
    scale_color_manual(values = c("GWAS" = "blue",
                                  "eQTL" = "darkorange")) +
    scale_x_continuous(
      name = paste0(toupper(substr(species,1,1)),
                    substr(species,2,nchar(species)),
                    " Chromosome"),
      breaks = genome$midpoint / 1e6,
      labels = genome$chr,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      name = "Sqrt(Proportion of Total Peaks)",
      expand = expansion(mult = c(0, 0)),
      limits = c(0, 0.5)
    ) +
    labs(
      title = main,
      color = NULL
    ) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", size = 9),
      axis.text = element_text(color = "black", size = 9),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = c(0.5, 0.95), legend.direction = "horizontal"
    )
  
  return(p1)
}

make_windows <- function(chr, chr_len, window_size, step_size) 
{
  starts <- seq(1, chr_len - window_size, by = step_size)
  ends <- starts + window_size - 1
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends),
          seqlengths = setNames(chr_len, chr))
}



return_Hotspot_df_Continuous <- function(peaks_df, 
                                              window_size = 1e6, 
                                              step_size = 100, 
                                              main = NULL, 
                                              species = 'maize') {
  
  # Genome definition
  
  genome <- if (species == 'maize') {
    data.frame(
      chr = as.character(1:10),
      end = c(308452471, 243675191, 238017767, 250330460, 226353449,
              181357234, 185808916, 182411202, 163004744, 152435371),
      start = rep(1, 10)
    )
  } else if (species == 'sorghum') {
    data.frame(
      chr = as.character(1:10),
      end = c(85112863, 79114963, 80873341, 71215609, 77058072,
              62713908, 68911884, 65779274, 63277606, 62870657),
      start = rep(1, 10)
    )
  } else {
    stop("Species must be 'maize' or 'sorghum'")
  }
  
  genome <- genome %>%
    arrange(as.numeric(chr)) %>%
    mutate(cumstart = cumsum(c(0, head(end, -1))),
           midpoint = cumstart + end / 2)
  
  
  # Create sliding windows
  
  window_list <- mapply(make_windows,
                        chr = genome$chr,
                        chr_len = genome$end,
                        window_size = window_size,
                        step_size = step_size,
                        SIMPLIFY = FALSE)
  
  windows <- do.call(c, unname(window_list))
  
  
  # Split by peak_type
  
  peak_list <- split(peaks_df, peaks_df$peak_type)
  
  combined_df <- do.call(rbind, lapply(names(peak_list), function(ptype) {
    
    df <- peak_list[[ptype]]
    
    peaks_gr <- GRanges(seqnames = df$CHROM,
                        ranges = IRanges(start = df$POS, end = df$POS))
    
    counts <- countOverlaps(windows, peaks_gr)
    
    temp <- as.data.frame(windows) %>%
      mutate(
        chr = as.character(seqnames),
        peaks_count = counts,
        peak_type = ptype
      )
    
    temp
  }))
  
  
  # Add cumulative genome position
  
  combined_df <- combined_df %>%
    left_join(genome %>% select(chr, cumstart), by = "chr") %>%
    mutate(genome_pos = start + cumstart,
           genome_pos_Mb = genome_pos / 1e6)
  
  
  # Plot
  
  p1 <- ggplot(combined_df,
               aes(x = genome_pos_Mb,
                   y = peaks_count,
                   color = peak_type)) +
    geom_line(linewidth = 0.7, alpha = 0.7) +
    geom_vline(xintercept = genome$cumstart[-1] / 1e6,
               linetype = "dotted",
               color = "grey50") +
    scale_color_manual(values = c("GWAS" = "blue",
                                  "eQTL" = "orange")) +
    scale_x_continuous(
      name = paste0(toupper(substr(species,1,1)),
                    substr(species,2,nchar(species)),
                    " Chromosome"),
      breaks = genome$midpoint / 1e6,
      labels = genome$chr,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(limits = c(0, 80)) +
    labs(
      y = "Number of Peaks",
      title = main,
      color = NULL
    ) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", size = 9),
      axis.text = element_text(color = "black", size = 9),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "bottom"
    )
  
  return(combined_df)
  
  
}
