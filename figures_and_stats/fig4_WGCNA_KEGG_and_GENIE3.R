library(tidyverse)
library(WGCNA)
library(ggplot2)
library(gprofiler2)
library(ggridges)
library(viridis)
library(stringr)
library(multcompView)
library(tidytext)
library(ggpubr)
library(patchwork)
library(scales)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(stats)

options(stringsAsFactors = FALSE)

setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")


#  WGCNA RESULTS - MAIZE


# Load WGCNA objects
load("processed_outputs/WGCNA/maize_WGCNA_results.RData")  # moduleLabels, MEs
maize_exp_data <- readRDS("expression_data/for_coexpression_analysis/maize_WGCNA_ready.rds")
maize_module_colors <- labels2colors(moduleLabels)

# Build gene-module map
maize_gene_module_map <- data.frame(
  GENE = colnames(maize_exp_data),
  ModuleLabel = moduleLabels,
  ModuleColor = maize_module_colors
)

# Create module lists for enrichment
maize_module_lists <- maize_gene_module_map %>%
  group_by(ModuleColor) %>%
  summarise(Genes = list(GENE)) %>%
  deframe()  # Named list: module color -> vector of genes



#  WGCNA ENRICHMENT - MAIZE


maize_enrichment_results <- map_dfr(names(maize_module_lists), function(mod) {
  
  genes <- maize_module_lists[[mod]]
  
  gost_res <- gost(
    query = genes,
    organism = "zmays",
    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
    correction_method = "fdr",
    user_threshold = 0.05,
    domain_scope = "custom",
    custom_bg = colnames(maize_exp_data)
  )
  
  if(is.null(gost_res$result)) return(tibble())
  
  gost_res$result %>%
    as_tibble() %>%
    mutate(Module = mod)
})

# Fold enrichment calculation
maize_enrichment_clean <- maize_enrichment_results %>%
  mutate(fold_enrichment = (intersection_size / query_size) /
           (term_size / effective_domain_size))

# Save results
write_csv(maize_enrichment_clean, "gprofiler_maize_WGCNA_module_enrichment.csv")



#  EXAMPLE DOTPLOT - ONE MODULE


mod_to_plot <- "black"
dotplot_data <- maize_enrichment_clean %>%
  filter(Module == mod_to_plot) %>%
  slice_max(order_by = -p_value, n = 20) %>%
  mutate(term_name = fct_reorder(term_name, -p_value))

ggplot(dotplot_data, aes(
  x = -log10(p_value),
  y = term_name,
  size = intersection_size,
  color = fold_enrichment
)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma") +
  labs(
    x = "-log10(p)",
    y = "GO/KEGG term",
    size = "Module genes",
    color = "fold_enrichment"
  ) +
  theme_classic(base_size = 12)



#  WGCNA ENRICHMENTS - MULTIPLE SPECIES


# Maize
maize_enrichment <- run_wgcna_enrichment(
  wgcna_rdata = "processed_outputs/WGCNA/maize_WGCNA_results.RData",
  expression_rds = "expression_data/for_coexpression_analysis/maize_WGCNA_ready.rds",
  organism = "zmays",
  species_name = "maize"
)
write_csv(maize_enrichment, "data/Enrichment_data/WGCNA_module_enrichments/maize_WGCNA_module_enrichments.csv")

# Sorghum mapping
sorghum_mapping <- read_tsv("data/gene_data/Sbicolor_730_v5.1.locus_transcript_name_map.txt") %>%
  select(`new-locusName`, `old-locusName`) %>%
  rename(V5_ID = `new-locusName`, V3_ID = `old-locusName`) %>%
  mutate(V3_ID = str_replace(V3_ID, "Sobic.", "SORBI_3"))

# Sorghum
sorghum_enrichment <- run_wgcna_enrichment(
  wgcna_rdata = "processed_outputs/WGCNA/sorghum_WGCNA_results.RData",
  expression_rds = "expression_data/for_coexpression_analysis/sorghum_WGCNA_ready.rds",
  organism = "sbicolor",
  species_name = "sorghum",
  id_conversion = "sorghum",
  mapping_df = sorghum_mapping
)
write_csv(sorghum_enrichment, "data/Enrichment_data/WGCNA_module_enrichments/sorghum_WGCNA_module_enrichments.csv")

# Soybean
soybean_enrichment <- run_wgcna_enrichment(
  wgcna_rdata = "processed_outputs/WGCNA/soybean_WGCNA_results.RData",
  expression_rds = "expression_data/for_coexpression_analysis/soybean_WGCNA_ready.rds",
  organism = "gmax",
  species_name = "soybean",
  id_conversion = "soybean"
)
write_csv(soybean_enrichment, "data/Enrichment_data/WGCNA_module_enrichments/soybean_WGCNA_module_enrichments.csv")



#  FUNGAL-MODULE CORRELATIONS
top_n_modules <- 8
r2_upper_limit <- 1

maize_df <- process_species(
  "Maize",
  "processed_outputs/WGCNA/maize_WGCNA_results.RData",
  "data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2020_maize_COMBINED_ALL_LEVELS.csv",
  "TWAS_results/sig_genes/maize_688_TWAS_all_sig_hits.csv",
  "data/Enrichment_data/WGCNA_module_enrichments/maize_WGCNA_module_enrichments.csv"
)

sorghum_df <- process_species(
  "Sorghum",
  "processed_outputs/WGCNA/sorghum_WGCNA_results.RData",
  "data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv",
  "TWAS_results/sig_genes/sorghum_736_TWAS_all_sig_hits.csv",
  "data/Enrichment_data/WGCNA_module_enrichments/sorghum_WGCNA_module_enrichments.csv"
)

soybean_df <- process_species(
  "Soybean",
  "processed_outputs/WGCNA/soybean_WGCNA_results.RData",
  "data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_soybean_COMBINED_ALL_LEVELS.csv",
  "TWAS_results/sig_genes/soybean_620_TWAS_all_sig_hits.csv",
  "data/Enrichment_data/WGCNA_module_enrichments/soybean_WGCNA_module_enrichments.csv"
)

# Combine all species for plotting
combined_df <- bind_rows(maize_df, sorghum_df, soybean_df) %>%
  group_by(Species) %>%
  mutate(Module_f = reorder_within(Module, R2, Species, fun = median, desc = TRUE)) %>%
  ungroup() %>%
  mutate(Species = factor(Species, levels = c("Sorghum", "Maize", "Soybean")))

# ANOVA letters
anova_letters_df <- combined_df %>%
  group_by(Species) %>%
  group_modify(~{
    aov_model <- aov(R2 ~ Module, data = .x)
    tukey <- TukeyHSD(aov_model)
    letters <- multcompLetters4(aov_model, tukey)
    tibble(Module = names(letters$Module$Letters),
           letters = letters$Module$Letters)
  }) %>%
  ungroup()

# N labels
n_labels_df <- combined_df %>% distinct(Species, Module, n_genes)
kegg_labels_df <- combined_df %>% distinct(Species, Module, top_term)

# Add factor info
anova_letters_df <- anova_letters_df %>%
  left_join(combined_df %>% select(Species, Module, Module_f) %>% distinct(), by = c("Species", "Module"))
n_labels_df <- n_labels_df %>%
  left_join(combined_df %>% select(Species, Module, Module_f) %>% distinct(), by = c("Species", "Module"))
kegg_labels_df <- kegg_labels_df %>%
  left_join(combined_df %>% select(Species, Module, Module_f) %>% distinct(), by = c("Species", "Module"))



#  VIOLIN PLOT


violin_plot <- ggplot(combined_df,
                      aes(y = Module_f, x = R2, fill = prop_TWAS)) +
  
  geom_violin(color = NA, alpha = 0.9, width = 0.9) +
  geom_boxplot(width = 0.15, color = "black") +
  
  facet_wrap(~Species, ncol = 1, scales = "free_y") +
  scale_y_reordered() +
  scale_fill_viridis_c(option = "C", limits = c(0, 1), name = "Proportion TWAS genes") +
  
  # Labels
  geom_text(data = n_labels_df, aes(x = -0.01, y = Module_f, label = paste0("N=", n_genes)),
            inherit.aes = FALSE, hjust = 1, size = 3) +
  geom_text(data = kegg_labels_df, aes(x = 0.52, y = Module_f, label = top_term),
            inherit.aes = FALSE, hjust = 0, size = 3) +
  geom_text(data = anova_letters_df, aes(x = 0.47, y = Module_f, label = letters),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  
  coord_cartesian(xlim = c(-0.05, 1), clip = "off") +
  
  theme_classic() +
  theme(strip.text = element_text(face = "bold"),
        plot.margin = margin(5.5, 5, 5.5, 1),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.33, 0.5),
        legend.text = element_text(size = 9, color = "black"),
        legend.title = element_text(size = 9, color = "black", hjust = 0.5),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.4, "cm"),
        legend.spacing.y = unit(2, "pt"),
        legend.margin = margin(2, 2, 2, 2)) +
  
  labs(y = NULL, x = expression("Module–fungal association (" * R^2 * ")"))

ggsave("graphs/test_WGCNA_corr_violins_with_KEGG.png",
       plot = violin_plot, units = "in", dpi = 300, height = 5.5, width = 6.5, bg = "white")



#  GENIE3 CENTRALITY


top_percent <- 0.01
sorghum_TWAS_genes <- read_lines("TWAS_results/sig_genes/sorghum_TWAS_sig_gene_list.txt")
maize_TWAS_genes <- read_lines("TWAS_results/sig_genes/maize_TWAS_sig_gene_list.txt")
soybean_TWAS_genes <- read_lines("TWAS_results/sig_genes/soybean_TWAS_sig_gene_list.txt")

# Function to process one species
process_genie3 <- function(weight_matrix, twas_genes, species_name) {
  # Comment: Convert weight matrix to long edge list and calculate centrality
  edge_df <- as.data.frame(as.table(weight_matrix)) %>%
    rename(Regulator = Var1, Target = Var2, Weight = Freq) %>%
    mutate(Regulator = as.character(Regulator), Target = as.character(Target)) %>%
    filter(Regulator != Target)
  
  weight_threshold <- quantile(edge_df$Weight, probs = 1 - top_percent)
  edge_df <- edge_df %>% filter(Weight >= weight_threshold)
  
  out_strength <- edge_df %>% group_by(Regulator) %>% summarise(out_strength = sum(Weight), .groups = "drop") %>% rename(Gene = Regulator)
  in_strength  <- edge_df %>% group_by(Target)    %>% summarise(in_strength  = sum(Weight), .groups = "drop") %>% rename(Gene = Target)
  
  centrality_df <- full_join(out_strength, in_strength, by = "Gene") %>%
    mutate(out_strength = replace_na(out_strength, 0),
           in_strength  = replace_na(in_strength, 0),
           total_strength = out_strength + in_strength,
           TWAS_status = ifelse(Gene %in% twas_genes, "TWAS", "Non-TWAS"),
           Species = species_name)
  return(centrality_df)
}

# Load GENIE3 matrices
sorghum_weight_matrix <- readRDS("processed_outputs/GENIE3/sorghum_GENIE3_weight_matrix_all_genes.rds")
maize_weight_matrix   <- readRDS("processed_outputs/GENIE3/maize_GENIE3_weight_matrix_all_genes.rds")
soybean_weight_matrix <- readRDS("processed_outputs/GENIE3/soybean_GENIE3_weight_matrix_all_genes.rds")

# Process all species
sorghum_genie3_df <- process_genie3(sorghum_weight_matrix, sorghum_TWAS_genes, "Sorghum")
maize_genie3_df   <- process_genie3(maize_weight_matrix,   maize_TWAS_genes,   "Maize")
soybean_genie3_df <- process_genie3(soybean_weight_matrix, soybean_TWAS_genes, "Soybean")

combined_genie3_df <- bind_rows(sorghum_genie3_df, maize_genie3_df, soybean_genie3_df) %>%
  mutate(out_strength_plot = pmin(out_strength, quantile(out_strength, 1)),
         in_strength_plot  = pmin(in_strength, quantile(in_strength, 1))) %>%
  mutate(Species = factor(Species, levels = c("Sorghum", "Maize", "Soybean")))

# Plot GENIE3 out- and in-strength
p_out <- ggplot(combined_genie3_df, aes(x = TWAS_status, y = out_strength, fill = TWAS_status)) +
  geom_violin(alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.2, color = "black", outlier.shape = NA) +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("TWAS","Non-TWAS")), label = "p.signif", hide.ns = FALSE) +
  facet_wrap(~Species, ncol = 3) +
  theme_classic() +
  labs(x = NULL, y = "Out-Strength (Top 1% edges)") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none",
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1))

p_in <- ggplot(combined_genie3_df, aes(x = TWAS_status, y = in_strength, fill = TWAS_status)) +
  geom_violin(alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.2, color = "black", outlier.shape = NA) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  scale_y_log10(labels = label_number()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("TWAS","Non-TWAS")), label = "p.signif", hide.ns = FALSE) +
  facet_wrap(~Species, ncol = 3) +
  theme_classic() +
  labs(x = NULL, y = "In-Strength (Top 1% edges)") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none",
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1))


# Generating some stats
genie3_stats <- combined_genie3_df %>%
  group_by(Species) %>%
  summarise(
    mean_out_TWAS = mean(out_strength[TWAS_status == "TWAS"], na.rm = TRUE),
    mean_out_NonTWAS = mean(out_strength[TWAS_status == "Non-TWAS"], na.rm = TRUE),
    mean_in_TWAS = mean(in_strength[TWAS_status == "TWAS"], na.rm = TRUE),
    mean_in_NonTWAS = mean(in_strength[TWAS_status == "Non-TWAS"], na.rm = TRUE),
    wilcox_out_p = wilcox.test(
      out_strength[TWAS_status == "TWAS"],
      out_strength[TWAS_status == "Non-TWAS"]
    )$p.value,
    wilcox_in_p = wilcox.test(
      in_strength[TWAS_status == "TWAS"],
      in_strength[TWAS_status == "Non-TWAS"]
    )$p.value,
    .groups = "drop"
  )

# Printing said stats
genie3_stats %>%
  mutate(
    wilcox_out_p = signif(wilcox_out_p, 3),
    wilcox_in_p = signif(wilcox_in_p, 3)
  ) %>%
  print()


p_genie3_out_and_in <- p_out | p_in
p_genie3_out_and_in


# Example GRN from soybean


#soybean_gene_module_map <- build_module_map(
#  wgcna_rdata = "processed_outputs/WGCNA/soybean_WGCNA_results.RData",
#  expression_rds = "expression_data/for_coexpression_analysis/soybean_WGCNA_ready.rds",
#  species_name = "Soybean",
#  twas_genes = soybean_TWAS_genes
#  )


# For soybean:
#soybean_TF_list <- read_table("expression_data/for_coexpression_analysis/Gma_TF_list.txt") %>%
#  pull(Gene_ID)
#soybean_weight_matrix <- readRDS(
#  "processed_outputs/GENIE3/soybean_GENIE3_weight_matrix_all_genes.rds"
#)

#soybean_module_plot <- plot_module_grn(
#  species_name = "Soybean",
#  weight_matrix = soybean_weight_matrix,
#  gene_module_map = soybean_gene_module_map,
#  tfs_list = soybean_TF_list,
#  twas_genes = soybean_TWAS_genes,
#  module_color = "paleturquoise",   # replace with desired module
#  top_k = 3
#)

#soybean_module_plot





#  COMBINED FIGURE 4
bottom_fig4_row <- (p_genie3_out_and_in) 

fig4_plot <- wrap_plots(
  violin_plot,
  bottom_fig4_row,
  ncol = 1, heights = c(11,7)
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.margin = margin(1,1,1,1))

fig4_plot

ggsave("graphs/fig4_plot_no_legend.svg", plot = fig4_plot, units = "in", dpi = 300, height = 9, width = 6.5, bg = "white")
#ggsave("graphs/fig4_plot_panel_legend.svg", plot = soybean_module_plot, units = "in", dpi = 300, height = 9, width = 6.5, bg = "white")











### Stats for MS  ##

### Stats for module–fungal associations ###

# Module-Fungal Correlation Stats per module
module_correlation_stats <- combined_df %>%
  group_by(Species, Module) %>%
  summarise(
    median_R2 = median(R2, na.rm = TRUE),          # Median across all fungi in the module
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    max_R2 = max(R2, na.rm = TRUE),               # Max R2 for any fungus in the module
    n_significant = sum(p_adj < 0.05, na.rm = TRUE),
    n_fungi_tested = n(),
    prop_significant = n_significant / n_fungi_tested,
    n_genes = first(n_genes),
    prop_TWAS = first(prop_TWAS),
    top_term = first(top_term),
    top_fungal_R2 = max(R2, na.rm = TRUE),       # Explicit column for top individual fungus
    .groups = "drop"
  )

# Summary by species
module_correlation_summary <- module_correlation_stats %>%
  group_by(Species) %>%
  summarise(
    median_module_R2 = median(median_R2),
    mean_module_R2 = mean(mean_R2),
    max_module_R2 = max(max_R2),                  # Max across all modules
    median_top_fungal_R2 = median(top_fungal_R2),# Median of module top fungus R2
    median_prop_significant = median(prop_significant),
    mean_prop_significant = mean(prop_significant),
    n_modules = n(),
    n_modules_with_KEGG = sum(top_term != "No pathway significantly enriched"),
    prop_modules_with_KEGG = n_modules_with_KEGG / n_modules,
    .groups = "drop"
  )

# Correlation between TWAS enrichment and module R2
twas_vs_module_R2 <- module_correlation_stats %>%
  group_by(Species) %>%
  summarise(
    spearman_rho = cor(prop_TWAS, median_R2, method = "spearman"),
    spearman_p = cor.test(prop_TWAS, median_R2, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )

# GENIE3 Centrality Stats
genie3_centrality_stats <- combined_genie3_df %>%
  group_by(Species, TWAS_status) %>%
  summarise(
    mean_out_strength = mean(out_strength, na.rm = TRUE),
    median_out_strength = median(out_strength, na.rm = TRUE),
    mean_in_strength = mean(in_strength, na.rm = TRUE),
    median_in_strength = median(in_strength, na.rm = TRUE),
    max_out_strength = max(out_strength, na.rm = TRUE),
    max_in_strength = max(in_strength, na.rm = TRUE),
    .groups = "drop"
  )

# Wilcoxon tests for in- and out-strength, plus direction
genie3_wilcox <- combined_genie3_df %>%
  group_by(Species) %>%
  summarise(
    wilcox_out_p = wilcox.test(out_strength[TWAS_status=="TWAS"],
                               out_strength[TWAS_status=="Non-TWAS"])$p.value,
    wilcox_out_direction = ifelse(
      median(out_strength[TWAS_status=="TWAS"]) > median(out_strength[TWAS_status=="Non-TWAS"]),
      "TWAS > Non-TWAS", "TWAS < Non-TWAS"),
    wilcox_in_p  = wilcox.test(in_strength[TWAS_status=="TWAS"],
                               in_strength[TWAS_status=="Non-TWAS"])$p.value,
    wilcox_in_direction = ifelse(
      median(in_strength[TWAS_status=="TWAS"]) > median(in_strength[TWAS_status=="Non-TWAS"]),
      "TWAS > Non-TWAS", "TWAS < Non-TWAS"),
    .groups = "drop"
  )

# Top module-fungal associations per species (explicitly highlighting strongest individual associations)
top_module_fungal <- combined_df %>%
  arrange(desc(R2)) %>%
  group_by(Species) %>%
  slice_head(n = 10) %>%
  select(Species, Module, Fungal_Taxon, R2, p_adj, n_genes, prop_TWAS, top_term)

# Per-species R2 distribution for reporting
species_R2_dist <- combined_df %>%
  group_by(Species) %>%
  summarise(
    median_R2 = median(R2, na.rm = TRUE),
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    max_R2 = max(R2, na.rm = TRUE),
    min_R2 = min(R2, na.rm = TRUE),
    .groups = "drop"
  )

# Combine into one summary list
figure4_stats <- list(
  module_correlation_stats = module_correlation_stats,
  module_correlation_summary = module_correlation_summary,
  twas_vs_module_R2 = twas_vs_module_R2,
  genie3_centrality_stats = genie3_centrality_stats,
  genie3_wilcox = genie3_wilcox,
  top_module_fungal = top_module_fungal,
  species_R2_dist = species_R2_dist
)

# Print for quick inspection
figure4_stats




####################


#### FUNCTIONS
# Run these before the rest of the script

# Convert sorghum Phytozome IDs to SORBI format
convert_sorghum_ids <- function(genes, mapping_df){
  
  tibble(V5_ID = genes) %>%
    left_join(mapping_df, by = "V5_ID") %>%
    pull(V3_ID)
}

# Convert soybean Glyma.xxxx to GLYMA_xxxx
convert_soybean_ids <- function(genes){
  genes %>%
    str_replace("^Glyma\\.", "GLYMA_")
}


run_wgcna_enrichment <- function(
    wgcna_rdata,
    expression_rds,
    organism,
    species_name,
    id_conversion = NULL,
    mapping_df = NULL
){
  
  # Load WGCNA results
  load(wgcna_rdata)  # loads moduleLabels
  exp_data <- readRDS(expression_rds)
  module_colors <- labels2colors(moduleLabels)
  
  # Build gene-module map
  gene_module_map <- data.frame(
    GENE = colnames(exp_data),
    ModuleLabel = moduleLabels,
    ModuleColor = module_colors
  )
  
  # Optional ID conversion
  if(!is.null(id_conversion)){
    
    if(id_conversion == "sorghum"){
      
      gene_module_map <- gene_module_map %>%
        left_join(mapping_df, by = c("GENE" = "V5_ID")) %>%
        mutate(GENE = V3_ID) %>%
        select(-V3_ID)
    }
    
    if(id_conversion == "soybean"){
      
      gene_module_map <- gene_module_map %>%
        mutate(GENE = str_replace(GENE, "^Glyma\\.", "GLYMA_"))
    }
    
    # Drop genes that failed to map
    gene_module_map <- gene_module_map %>%
      filter(!is.na(GENE))
  }
  
  # Remove NA genes after conversion
  gene_module_map <- gene_module_map %>%
    filter(!is.na(GENE))
  
  # Create module lists
  module_lists <- gene_module_map %>%
    group_by(ModuleColor) %>%
    summarise(Genes = list(GENE)) %>%
    deframe()
  
  background_genes <- unique(gene_module_map$GENE)
  
  # Run g:Profiler
  enrichment_results <- map_dfr(names(module_lists), function(mod){
    
    genes <- module_lists[[mod]]
    
    gost_res <- gost(
      query = genes,
      organism = organism,
      sources = c("GO:BP","GO:MF","GO:CC","KEGG"),
      correction_method = "fdr",
      user_threshold = 0.05,
      domain_scope = "custom",
      custom_bg = background_genes
    )
    
    if(is.null(gost_res$result)) return(tibble())
    
    gost_res$result %>%
      as_tibble() %>%
      mutate(
        Module = mod,
        Species = species_name
      )
  })
  
  # Calculate fold enrichment
  enrichment_results <- enrichment_results %>%
    mutate(
      fold_enrichment =
        (intersection_size / query_size) /
        (term_size / effective_domain_size)
    )
  
  
  return(enrichment_results)
}






compute_correlations <- function(ME_mat, fungal_mat){
  
  corMat <- cor(ME_mat, fungal_mat,
                use = "pairwise.complete.obs",
                method = "pearson")
  
  nMat <- matrix(NA, nrow = nrow(corMat), ncol = ncol(corMat))
  
  for(i in seq_len(nrow(corMat))){
    for(j in seq_len(ncol(corMat))){
      nMat[i,j] <- sum(complete.cases(ME_mat[,i], fungal_mat[,j]))
    }
  }
  
  df <- nMat - 2
  
  pMat <- 2 * pt(
    -abs(corMat) * sqrt(df / (1 - corMat^2)),
    df = df
  )
  
  pAdj <- apply(pMat, 2, p.adjust, method = "fdr")
  
  cor_df <- as.data.frame(corMat) %>%
    rownames_to_column("Module") %>%
    pivot_longer(-Module,
                 names_to = "Fungal_Taxon",
                 values_to = "r") %>%
    mutate(
      R2 = r^2,
      p_value = as.vector(pMat),
      p_adj = as.vector(pAdj)
    )
  
  return(cor_df)
}


get_top_term <- function(enrich_df, module_color){
  
  candidate <- enrich_df %>%
    filter(Module == module_color,
           source == "KEGG",
           p_value <= 0.05,
           term_name != "Viral life cycle - HIV-1") %>%  # skip artifact) %>%
    arrange(p_value, desc(fold_enrichment)) %>%
    #arrange(desc(fold_enrichment)) %>%
    slice_head(n = 1)
  
  if(nrow(candidate) == 0){
    return("No pathway significantly enriched")
  }
  
  return(candidate$term_name)
}


process_species <- function(
    species_name,
    wgcna_rdata,
    fungal_file,
    twas_file,
    enrichment_file
){
  
  message("Processing: ", species_name)
  
  ################ Load TWAS
  twas_genes <- read_csv(twas_file, show_col_types = FALSE) %>%
    pull(unique(SNP))
  
  ################ Load WGCNA
  load(wgcna_rdata)  # loads moduleLabels + MEs
  
  gene_module_map <- tibble(
    GENE = names(moduleLabels),
    ModuleLabel = moduleLabels,
    ModuleColor = labels2colors(moduleLabels),
    ME = paste0("ME", moduleLabels)
  ) %>%
    mutate(TWAS_sig = GENE %in% twas_genes)
  
  ################ Load fungal abundance
  fungal_data <- read_csv(fungal_file, show_col_types = FALSE)
  
  MEs_df <- MEs %>%
    rownames_to_column("GENOTYPE") %>%
    mutate(GENOTYPE = str_remove(GENOTYPE, "^4\\d{3}_")) %>%
    arrange(GENOTYPE)
  
  fungal_data <- fungal_data %>%
    arrange(GENOTYPE)
  
  ME_mat <- MEs_df %>%
    select(-GENOTYPE) %>%
    as.matrix()
  
  fungal_mat <- fungal_data %>%
    select(-GENOTYPE) %>%
    as.matrix()
  
  ################ Correlations
  cor_results <- compute_correlations(ME_mat, fungal_mat)
  
  sig_results <- cor_results %>%
    filter(p_adj < 0.05)
  
  module_rank <- sig_results %>%
    group_by(Module) %>%
    summarise(median_R2 = median(R2, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(median_R2)) %>%
    slice_head(n = top_n_modules)
  
  top_modules <- module_rank$Module
  
  ################ Module TWAS summary
  module_summary <- gene_module_map %>%
    group_by(ME, ModuleColor) %>%
    summarise(
      n_genes = n(),
      n_TWAS = sum(TWAS_sig),
      prop_TWAS = n_TWAS / n_genes,
      .groups = "drop"
    ) %>%
    rename(Module = ME)
  
  ################ Load enrichment
  enrichment_df <- read_csv(enrichment_file, show_col_types = FALSE)
  
  ################ Attach enrichment labels
  module_summary <- module_summary %>%
    rowwise() %>%
    mutate(
      top_term = get_top_term(enrichment_df, ModuleColor)
    ) %>%
    ungroup()
  
  ################ Final dataset
  plot_df <- cor_results %>%
    filter(Module %in% top_modules) %>%
    left_join(module_summary, by = "Module") %>%
    left_join(module_rank, by = "Module") %>%
    mutate(
      Species = species_name,
      Module = factor(Module,
                      levels = module_rank$Module)
    )
  
  return(plot_df)
}



#  EXAMPLE GRN PLOT FUNCTION

plot_module_grn <- function(
    species_name,
    weight_matrix,          # GENIE3 weight matrix
    gene_module_map,        # gene-module map (with ModuleColor, TWAS_sig)
    tfs_list,               # list of TFs in species
    twas_genes,             # vector of TWAS-significant genes
    module_color,           # which module to plot
    top_k = 3               # top k incoming edges per target
){
  
  message("Building GRN for module: ", module_color)
  
  # Select genes in module
  module_genes <- gene_module_map %>%
    filter(ModuleColor == module_color) %>%
    pull(GENE)
  
  # Convert weight matrix to long format and filter to module genes
  edges <- as.data.frame(as.table(weight_matrix)) %>%
    rename(regulator = Var1, target = Var2, weight = Freq) %>%
    filter(
      regulator %in% module_genes,
      target %in% module_genes,
      weight > 0
    ) %>%
    rename(Weight = weight)
  
  # Keep top k edges per target
  edges_topk <- edges %>%
    group_by(target) %>%
    slice_max(Weight, n = top_k) %>%
    ungroup()
  
  # Build igraph
  net <- graph_from_data_frame(edges_topk, directed = TRUE)
  
  # Convert to tbl_graph and annotate nodes
  graph_tbl <- as_tbl_graph(net) %>%
    activate(nodes) %>%
    mutate(
      is_TF = name %in% tfs_list,
      TWAS_sig = name %in% twas_genes,
      degree = centrality_degree(mode = "all")
    ) %>%
    mutate(
      TWAS_status = ifelse(TWAS_sig, "Associated", "Background")
    ) 
  
  # Plot GRN
  grn_plot <- ggraph(graph_tbl, layout = "kk") +
    
    # Edges
    geom_edge_link(aes(edge_alpha = Weight), color = "grey70") +
    
    # Non-TF nodes
    geom_node_point(
      data = function(x) dplyr::filter(x, !is_TF),
      aes(size = degree, color = TWAS_status),
      shape = 16,
      alpha = 0.9
    ) +
    
    # TF nodes
    geom_node_point(
      data = function(x) dplyr::filter(x, is_TF),
      aes(size = degree, color = TWAS_status),
      shape = 17,
      alpha = 0.95
    ) +
    
    # Labels only for TWAS genes
    #    geom_text_repel(
    #      data = function(x) dplyr::filter(x, TWAS_sig),
    #      aes(x = x, y = y, label = name),
    #      size = 3,
    #      box.padding = 1,
    #      point.padding = 0.5,
    #      segment.color = "black",
    #      segment.size = 0.25,
    #      max.overlaps = Inf,
    #      show.legend = FALSE
    #    ) +
    
    # Color scale for TWAS
    scale_color_manual(
      values = c("Background" = "#009E73",
                 "Associated" = "#E69F00"),
      name = "TWAS status"
    ) +
    
    # Node size scale
    scale_size_continuous(
      range = c(2, 8),
      name = "Degree"
    ) +
    
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9, color = "black"),
      legend.text = element_text(size = 9, color = "black")
    )
  #return(edges_topk)
  return(grn_plot)
}


build_module_map <- function(
    wgcna_rdata,
    expression_rds,
    species_name,
    twas_genes = NULL,
    id_conversion = NULL,
    mapping_df = NULL
){
  # Load WGCNA results
  load(wgcna_rdata)  # expects object moduleLabels
  exp_data <- readRDS(expression_rds)
  
  # Convert numeric labels to colors
  module_colors <- labels2colors(moduleLabels)
  
  # Base gene-module map
  gene_module_map <- tibble(
    GENE = colnames(exp_data),
    ModuleLabel = moduleLabels,
    ModuleColor = module_colors
  )
  
  # Optional ID conversion
  if(!is.null(id_conversion)){
    if(id_conversion == "sorghum"){
      gene_module_map <- gene_module_map %>%
        left_join(mapping_df, by = c("GENE" = "V5_ID")) %>%
        mutate(GENE = V3_ID) %>%
        select(-V3_ID)
    }
    if(id_conversion == "soybean"){
      gene_module_map <- gene_module_map %>%
        mutate(GENE = str_replace(GENE, "^Glyma\\.", "GLYMA_"))
    }
    gene_module_map <- gene_module_map %>% filter(!is.na(GENE))
  }
  
  # Add TWAS annotation if provided
  if(!is.null(twas_genes)){
    gene_module_map <- gene_module_map %>%
      mutate(TWAS_sig = GENE %in% twas_genes)
  } else {
    gene_module_map <- gene_module_map %>%
      mutate(TWAS_sig = FALSE)
  }
  
  return(gene_module_map)
}

