library(tidyverse)
library(ggplot2)
library(gprofiler2)
library(viridis)
library(stringr)
library(multcompView)
library(tidytext)
library(ggpubr)
library(forcats)
library(ComplexUpset)
library(patchwork)
library(cowplot)
library(clusterProfiler)
library(GOSemSim)



setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")



########### TWAS Orthogroups #######
# LOAD ORTHOGROUPS 
orthogroups <- read_tsv("data/orthofinder/Orthogroups.tsv")

# LOAD TWAS GENE LISTS 
maize_twas <- read_csv("TWAS_results/sig_genes/maize_688_TWAS_all_sig_hits.csv") %>%
  select(SNP, Phenotype) %>% distinct(SNP) %>% pull(SNP)

sorghum_twas <- read_csv("TWAS_results/sig_genes/sorghum_736_TWAS_all_sig_hits.csv") %>%
  select(SNP, Phenotype) %>% distinct(SNP) %>% pull(SNP)

soybean_twas <- read_csv("TWAS_results/sig_genes/soybean_620_TWAS_all_sig_hits.csv") %>%
  select(SNP, Phenotype) %>% distinct(SNP) %>% pull(SNP)

#write_lines(maize_twas, "TWAS_results/sig_genes/maize_TWAS_sig_gene_list.txt")
#write_lines(sorghum_twas, "TWAS_results/sig_genes/sorghum_TWAS_sig_gene_list.txt")
#write_lines(soybean_twas, "TWAS_results/sig_genes/soybean_TWAS_sig_gene_list.txt")

maize_twas_background_genes <- read_lines("TWAS_results/sig_genes/maize_TWAS_background_list.txt")
sorghum_twas_background_genes <- read_lines("TWAS_results/sig_genes/sorghum_TWAS_background_list.txt")
soybean_twas_background_genes <- read_lines("TWAS_results/sig_genes/soybean_TWAS_background_list.txt")



# ORTHOGROUP → LONG =
og_long <- orthogroups %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Genes") %>%
  separate_rows(Genes, sep = ", ") %>%
  mutate(Genes = str_remove(Genes, "_P\\d+$")) %>%
  mutate(Genes = str_remove(Genes, ".1.p|.2.p|.3.p|.4.p|.5.p|.6.p|.7.p|.8.p|.9.p|.10.p"))

# MARK TWAS HITS 
og_twas <- og_long %>%
  mutate(
    TWAS = case_when(
      Species == "Zmays_v5_protein_clean" & Genes %in% maize_twas ~ TRUE,
      Species == "Sbicolor_v5_protein_clean" & Genes %in% sorghum_twas ~ TRUE,
      Species == "Gmax_Wm82_v2_protein_clean" & Genes %in% soybean_twas ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(TWAS)

# Print og_TWAS
og_twas_print <- og_twas %>%
  mutate(Species = str_remove(Species, "_protein_clean")) %>%
  select(!TWAS)

write_csv(og_twas_print, "data/orthofinder/orthogroups_of_TWAS_genes.csv")

# PRESENCE MATRIX 
og_matrix <- og_twas %>%
  select(Orthogroup, Species) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Species, values_from = present, values_fill = 0) %>%
  dplyr::rename(Maize = Zmays_v5_protein_clean,
                Sorghum = Sbicolor_v5_protein_clean,
                Soybean = Gmax_Wm82_v2_protein_clean)

# OVERLAP STATS 
overlap_counts <- og_matrix %>%
  mutate(n_species = rowSums(across(-Orthogroup))) %>%
  count(n_species, name = "n_orthogroups")

print(overlap_counts)

 
# Clean UpSet plot
# Define all intersections explicitly so the 3-way intersection appears even if size = 0
all_intersections <- list(
  c("Maize"),
  c("Sorghum"),
  c("Soybean"),
  c("Maize", "Sorghum"),
  c("Maize", "Soybean"),
  c("Sorghum", "Soybean"),
  c("Maize", "Sorghum", "Soybean")
)

# Upset plot of TWAS ortholog overlap
upset_plot <- ComplexUpset::upset(
  og_matrix,
  intersect = c("Maize", "Sorghum", "Soybean"),
  intersections = all_intersections,
  
  base_annotations = list(
    "TWAS gene orthogroups" = intersection_size(
      fill = "#E69F00",      # clean green bars
      text = list(size = 3),
      theme = theme(
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9, color = "black")
      )
    )
  ),
  
  set_sizes = FALSE,        # removes horizontal set size bars
  width_ratio = 0.2
) +
  
  labs(
    x = NULL,
    y = "TWAS gene orthogroups"
  ) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_blank(),
    panel.border = element_blank()
  )



upset_plot


############### TWAS Volcano Plots ##########
# LOAD FULL TWAS TABLES 
maize_full <- read_csv("TWAS_results/sig_genes/maize_688_TWAS_all_sig_hits.csv") %>%
  mutate(SPECIES = "Maize")

sorghum_full <- read_csv("TWAS_results/sig_genes/sorghum_736_TWAS_all_sig_hits.csv") %>%
  mutate(SPECIES = "Sorghum")

soybean_full <- read_csv("TWAS_results/sig_genes/soybean_620_TWAS_all_sig_hits.csv") %>%
  mutate(SPECIES = "Soybean")

combined_full <- bind_rows(maize_full, sorghum_full, soybean_full)

TWAS_to_print <- combined_full %>%
  rename(Gene = SNP,
         Rsquare_of_Model_without_Gene = Rsquare.of.Model.without.SNP,
         Rsquare_of_Model_with_Gene = Rsquare.of.Model.with.SNP,
         Taxon = Phenotype) %>%
  mutate(taxon_ID = str_extract(Taxon, "[^_]+$")) %>%
  select(SPECIES, Taxon, taxon_ID, Gene, Chromosome, Position, P.value, nobs, Rsquare_of_Model_without_Gene, Rsquare_of_Model_with_Gene, effect, FDR) 

write_csv(TWAS_to_print, "processed_outputs/TWAS/all_species_sig_TWAS_hits_supplemental.csv")


# TWAS for Sobic.004G214900
Sobic_004G214900_TWAS <- TWAS_to_print %>%
  filter(Gene == "Sobic.004G214900")
write_csv(Sobic_004G214900_TWAS, "processed_outputs/TWAS/Sobic_004G214900_sig_TWAS_hits_supplemental.csv")

#
species_colors <- c("Maize"="#F0A202","Soybean"="#4DAF4A","Sorghum"="#7E57C2")

# TAXA COUNTS 
maize_taxa <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2020_maize_COMBINED_ALL_LEVELS.csv")
sorghum_taxa <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv")
soybean_taxa <- read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_soybean_COMBINED_ALL_LEVELS.csv")

maize_taxa_n <- ncol(maize_taxa)-1
sorghum_taxa_n <- ncol(sorghum_taxa)-1
soybean_taxa_n <- ncol(soybean_taxa)-1

# Make sure SNP is a character vector
combined_full <- combined_full %>%
  mutate(SNP = as.character(SNP))

# Then count taxa per gene
gene_counts <- combined_full %>%
  group_by(SNP) %>%   # include SPECIES if needed
  summarise(n_taxa_sig = n(), .groups="drop")

combined_counts <- combined_full %>%
  left_join(gene_counts, by="SNP") %>%
  mutate(prop_taxa_sig = case_when(
    SPECIES=="Maize" ~ n_taxa_sig/maize_taxa_n,
    SPECIES=="Sorghum" ~ n_taxa_sig/sorghum_taxa_n,
    SPECIES=="Soybean" ~ n_taxa_sig/soybean_taxa_n
  )) %>%
  distinct(SNP, SPECIES, prop_taxa_sig)

# PANEL B
p_taxa <- ggplot(combined_counts,
                 aes(SPECIES, prop_taxa_sig, fill=SPECIES)) +
  geom_violin(trim=FALSE, alpha=0.5) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  #geom_jitter(width=0.15, alpha=0.3, size=1) +
  scale_y_log10() +
  scale_fill_manual(values=species_colors) +
  theme_classic() +
  labs(y="Proportion of taxa affected", x="") +
  theme(legend.position = "none")
p_taxa

# Stats
combined_counts %>%
  group_by(SPECIES) %>%
  summarise(median=median(prop_taxa_sig),
            mean=mean(prop_taxa_sig),
            sd=sd(prop_taxa_sig),
            n=n()) %>% print()

print(kruskal.test(prop_taxa_sig ~ SPECIES, data=combined_counts))

# EFFECT SIZE delta R2 
combined_full <- combined_full %>%
  mutate(
    R2_gain_with_SNP = Rsquare.of.Model.with.SNP - Rsquare.of.Model.without.SNP,
    neglog10FDR = -log10(FDR)
  )

# PANEL C — delta R2 DISTRIBUTION 
p_effect <- ggplot(combined_full,
                   aes(SPECIES, R2_gain_with_SNP, fill=SPECIES)) +
  geom_violin(trim=FALSE, alpha=0.6) +
  geom_boxplot(width=0.15, outlier.shape=NA) +
  scale_fill_manual(values=species_colors) +
  theme_classic() +
  labs(y=expression(Delta*R^2~"(variance explained gain)"), x="") +
  theme(legend.position = "none")
p_effect
# Stats for delta R2
combined_full %>%
  group_by(SPECIES) %>%
  summarise(median=median(R2_gain_with_SNP, na.rm=TRUE),
            mean=mean(R2_gain_with_SNP, na.rm=TRUE),
            sd=sd(R2_gain_with_SNP, na.rm=TRUE),
            n=n()) %>% print()

# Merge in proportion of taxa affected
volcano_df <- combined_full %>%
  left_join(combined_counts, by = c("SNP", "SPECIES")) %>%
  mutate(
    effect_sign = case_when(
      effect > 0 ~ "Positive",
      effect < 0 ~ "Negative",
      TRUE ~ "neutral"
    ),
    # Cap very small prop_taxa_sig to avoid invisible points
    prop_taxa_sig = ifelse(prop_taxa_sig < 0.001, 0.001, prop_taxa_sig)
  )

# Volcano/trident plot
p_volcano <- ggplot(volcano_df, aes(x = effect, y = -log10(FDR))) +
  geom_point(aes(color = effect_sign, size = prop_taxa_sig), alpha = 0.6) +
  facet_wrap(~SPECIES, scales = "free_y") +
  scale_color_manual(values = c("Positive" = "#E41A1C",
                                "Negative" = "#377EB8",
                                "neutral"  = "#999999")) +
  scale_size_continuous(range = c(1, 4), name = "Prop. taxa affected") +
  theme_classic() +
  labs(x = "TWAS effect size", y = expression(-log[10](FDR)),
       color = "Effect sign") +
  theme(strip.text = element_text(face = "bold", size = 9, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9, color = "black"),
        legend.title = element_text(size = 9, color = "black"),
        legend.position = c(0.08, 0.6),
        legend.background = element_blank())

print(p_volcano)
#ggsave("graphs/transparent_bg_versions/TWAS_volcano_effect_vs_FDR.png", plot = p_trident, dpi = 600, units = "in", height = 6, width = 12, bg = "white")
#  SUMMARY STATS
volcano_summary <- volcano_df %>%
  filter(FDR < 0.1) %>%  # only significant-ish points
  group_by(SPECIES, effect_sign) %>%
  summarise(
    n_genes = n(),
    median_effect = median(effect),
    mean_prop_taxa = mean(prop_taxa_sig, na.rm = TRUE),
    max_neglogFDR = max(-log10(FDR)),
    .groups = "drop"
  )

print(volcano_summary)


# EXTRA STATS 
combined_full %>% count(SPECIES) %>% print()

combined_counts %>%
  left_join(combined_full, by=c("SNP","SPECIES")) %>%
  group_by(SPECIES) %>%
  summarise(correlation = cor(effect, prop_taxa_sig, method="spearman")) %>%
  print()  

# Does variance explained relate to breadth of taxa affected?
combined_counts %>%
  left_join(combined_full, by=c("SNP","SPECIES")) %>%
  group_by(SPECIES) %>%
  summarise(
    cor_R2_breadth = cor(R2_gain_with_SNP, prop_taxa_sig, method="spearman"),
    cor_effect_breadth = cor(effect, prop_taxa_sig, method="spearman"),
    cor_R2_sig = cor(R2_gain_with_SNP, neglog10FDR, method="spearman")
  ) %>% print()





############  More stats


summarize_twas_architecture <- function(df, fdr_thresh = 0.05) {
  
  library(dplyr)
  
  # Gini helper 
  gini <- function(x) {
    if(length(x) == 0 || sum(x) == 0) return(NA_real_)
    x <- sort(x)
    n <- length(x)
    sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
  }
  
  # Filter + compute delta R2 
  df_sig <- df %>%
    filter(FDR < fdr_thresh) %>%
    mutate(delta_R2 = Rsquare.of.Model.with.SNP - Rsquare.of.Model.without.SNP)
  
  if(nrow(df_sig) == 0) {
    return(tibble(
      n_genes = 0,
      mean_taxa_per_gene = NA,
      median_taxa_per_gene = NA,
      max_taxa_per_gene = NA,
      gini_taxa_per_gene = NA,
      prop_top5 = NA,
      mean_delta_R2 = NA,
      median_delta_R2 = NA
    ))
  }
  
  #  Gene-level summary 
  gene_summary <- df_sig %>%
    group_by(SNP) %>%
    summarise(
      n_taxa = n(),
      mean_delta_R2 = mean(delta_R2),
      .groups = "drop"
    )
  
  #  Top 5% contribution 
  top5_prop <- gene_summary %>%
    arrange(desc(n_taxa)) %>%
    mutate(rank = row_number(),
           prop_genes = rank / n()) %>%
    summarise(
      top5_taxa = sum(n_taxa[prop_genes <= 0.05]),
      total_taxa = sum(n_taxa),
      prop = top5_taxa / total_taxa
    ) %>%
    pull(prop)
  
  #  Final summary 
  tibble(
    n_genes = nrow(gene_summary),
    mean_taxa_per_gene = mean(gene_summary$n_taxa),
    median_taxa_per_gene = median(gene_summary$n_taxa),
    max_taxa_per_gene = max(gene_summary$n_taxa),
    gini_taxa_per_gene = gini(gene_summary$n_taxa),
    prop_top5 = top5_prop,
    mean_delta_R2 = mean(df_sig$delta_R2),
    median_delta_R2 = median(df_sig$delta_R2)
  )
}


maize_stats <- summarize_twas_architecture(maize_full)
sorghum_stats <- summarize_twas_architecture(sorghum_full)
soybean_stats <- summarize_twas_architecture(soybean_full)

# Combine into one table
bind_rows(
  Maize = maize_stats,
  Sorghum = sorghum_stats,
  Soybean = soybean_stats,
  .id = "Species"
)











################# TWAS Gene Distribution plots ########################



















##############  Running enrichments on TWAS gene sets   #################
# Maize
maize_enrichment <- run_twas_enrichment(
  twas_gene_file = "TWAS_results/sig_genes/maize_TWAS_sig_gene_list.txt",
  background_genes = "TWAS_results/sig_genes/maize_TWAS_background_list.txt",
  organism = "zmays",
  species_name = "Maize"
)

# Sorghum
sorghum_enrichment <- run_twas_enrichment(
  twas_gene_file = "TWAS_results/sig_genes/sorghum_TWAS_sig_gene_list.txt",
  background_genes = "TWAS_results/sig_genes/sorghum_TWAS_background_list.txt",
  organism = "sbicolor",
  species_name = "Sorghum",
  id_conversion = "sorghum",
  mapping_file = "data/gene_data/Sbicolor_730_v5.1.locus_transcript_name_map.txt"
)

# Soybean
soybean_enrichment <- run_twas_enrichment(
  twas_gene_file = "TWAS_results/sig_genes/soybean_TWAS_sig_gene_list.txt",
  background_genes = "TWAS_results/sig_genes/soybean_TWAS_background_list.txt",
  organism = "gmax",
  species_name = "Soybean",
  id_conversion = "soybean"
)


# Combine enrichment results from all species


all_enrichment1 <- bind_rows(
  maize_enrichment,
  sorghum_enrichment,
  soybean_enrichment
)

write_csv(all_enrichment1, "data/Enrichment_data/TWAS_gene_enrichments_all_species.csv")

# Collapse GO categories into one group (GO) and keep KEGG
all_enrichment <- read_csv("data/Enrichment_data/TWAS_gene_enrichments_all_species.csv")

all_enrichment <- all_enrichment %>%
  mutate(
    category = case_when(
      source %in% c("GO:BP","GO:MF","GO:CC") ~ "GO",
      source == "KEGG" ~ "KEGG",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category))# %>%
  #filter(source %in% c("GO:BP", "KEGG"))







top_df <- all_enrichment %>%
  group_by(Species, source) %>%
  slice_min(p_value, n = 10, with_ties = FALSE) %>% #Could select top X terms per species
  ungroup()



# Determine how many species share each enriched term


sharing_df <- top_df %>%
  distinct(category, term_name, Species) %>%
  group_by(category, term_name) %>%
  summarise(
    n_species = n(),
    .groups = "drop"
  )

plot_df <- top_df %>%
  left_join(sharing_df,
            by = c("category","term_name"))



# Add visual indicators for shared terms
# ● shared by all species
# ○ shared by two species


plot_df <- plot_df %>%
  mutate(
    share_label = case_when(
      n_species == 3 ~ "● ",
      n_species == 2 ~ "○ ",
      TRUE ~ ""
    ),
    
    term_display = paste0(
      share_label,
      stringr::str_wrap(term_name, 60)
    )
  )



# Rank terms by mean significance across species
# so ordering is consistent within each panel


ranking_df <- plot_df %>%
  group_by(category, term_display) %>%
  summarise(
    mean_fdr = mean(p_value),
    mean_neglog = -log10(mean_fdr),
    .groups = "drop"
  )

plot_df <- plot_df %>%
  left_join(ranking_df,
            by = c("category","term_display"))

plot_df <- plot_df %>%
  group_by(category) %>%
  mutate(
    term_display = forcats::fct_reorder(
      term_display,
      mean_neglog,
      .desc = FALSE
    )
  ) %>%
  ungroup()

fe_limits <- range(plot_df$fold_enrichment, na.rm = TRUE)
# Split enrichment data into GO and KEGG panels

go_plot_df <- plot_df %>%
  filter(category == "GO")

kegg_plot_df <- plot_df %>%
  filter(category == "KEGG")

species_shapes <- c("Maize" = 15,    # filled square
                    "Sorghum" = 17,  # filled triangle
                    "Soybean" = 16)  # filled circle

# Define the three ontologies
ontologies <- c("BP", "MF", "CC")

for (ont in ontologies) {
  
  # Filter for current ontology
  go_df_ont <- go_plot_df %>%
    mutate(source = str_replace(source, "GO:", "")) %>%  
    filter(source == ont)    
  
  # Plot
  p_enrich_go <- ggplot(
    go_df_ont,
    aes(x = -log10(p_value),
        y = term_display,
        color = fold_enrichment,
        shape = Species)
  ) +
    geom_point(size = 4.5) +
    scale_color_viridis_c(
      name = "Fold Enrichment",
      option = "C",
      limits = fe_limits,
      alpha = 0.7
    ) +  
    scale_shape_manual(values = species_shapes) +
    labs(
      x = expression(-log[10]("FDR")),
      y = NULL,
      title = paste0("GO Enrichment - ", ont)
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 9, color = "black")
    )
  
  # Save file
  ggsave(
    filename = paste0("graphs/combined_species_TWAS_GO_enrichment_", ont, ".png"),
    plot = p_enrich_go,
    units = "in",
    width = 9,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

# KEGG enrichment plot

p_enrich_kegg <- ggplot(
  kegg_plot_df,
  aes(x = -log10(p_value),
      y = term_display,
      color = fold_enrichment,
      shape = Species)
) +
  geom_point(size = 4.5) +
  scale_color_viridis_c(
    name = "Fold Enrichment",
    option = "C",
    limits = fe_limits,
    alpha = 0.7
  ) +  
  scale_shape_manual(values = species_shapes) +
  labs(
    x = expression(-log[10]("FDR")),
    y = NULL,
    title = "KEGG Enrichment"
  ) +
  scale_x_log10() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position ="bottom",
    legend.title = element_text(size = 9, color = "black")
  )
p_enrich_kegg
ggsave("graphs/combined_species_TWAS_KEGG_enrichment.png", plot = p_enrich_kegg, units = "in", width = 9, height = 6.5, dpi = 300, bg = "white")


############################################################



# Align enrichment plots
aligned_enrich <- align_plots(
  p_enrich_go,
  p_enrich_kegg,
  align = "h",
  axis = "l"
)
p_enrich_go  <- aligned_enrich[[1]]
p_enrich_kegg <- aligned_enrich[[2]]



# Align volcano + upset plots
aligned_top <- align_plots(
  p_volcano,
  upset_plot,
  align = "h",
  axis = "l"
)
p_volcano <- aligned_top[[1]]
upset_plot <- aligned_top[[2]]

# Top row with labels for each plot
top_row <- plot_grid(
  p_volcano,
  upset_plot,
  ncol = 1,
  rel_widths = c(1, 1),
  rel_heights = c(2, 1),
  labels = c("A","B"),
  label_size = 16,
  label_fontface = "bold"
)

# Bottom row with labels for each plot
bottom_row <- plot_grid(
  p_enrich_go,
  p_enrich_kegg,
  ncol = 1,
  rel_widths = c(1, 1),
  rel_heights = c(2, 1),
  labels = c("C","D"),
  label_size = 16,
  label_fontface = "bold"
)

# Combine rows into final figure
fig3_combined <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(4, 4)
)


fig3_combined

ggsave("graphs/temp_fig3_combined.svg", plot = top_row, dpi = 600, units = "in", height = 6, width = 6.5, bg = "white")


#









#############  Functions (Run first)  ############

run_twas_enrichment <- function(
    twas_gene_file,
    background_genes,
    organism,
    species_name,
    id_conversion = NULL,
    mapping_file = NULL
){
  
  
  # Load TWAS genes
  
  twas_genes <- read_lines(twas_gene_file)
  
  
  # Load expression data (for background)
  
  exp_data <- read_lines(background_genes)
  
  gene_df <- tibble(GENE = exp_data)
  
  
  # Optional ID conversion
  
  if(!is.null(id_conversion)){
    
    if(id_conversion == "sorghum"){
      
      mapping_df <- read_tsv(mapping_file) %>%
        select(`new-locusName`, `old-locusName`) %>%
        rename(V5_ID = `new-locusName`,
               V3_ID = `old-locusName`) %>%
        mutate(V3_ID = str_replace(V3_ID, "Sobic.", "SORBI_3"))
      
      # Convert background
      gene_df <- gene_df %>%
        left_join(mapping_df, by = c("GENE" = "V5_ID")) %>%
        mutate(GENE = V3_ID) %>%
        select(GENE)
      
      # Convert TWAS genes
      twas_genes <- tibble(GENE = twas_genes) %>%
        left_join(mapping_df, by = c("GENE" = "V5_ID")) %>%
        mutate(GENE = V3_ID) %>%
        pull(GENE)
    }
    
    if(id_conversion == "soybean"){
      
      gene_df <- gene_df %>%
        mutate(GENE = str_replace(GENE, "^Glyma\\.", "GLYMA_"))
      
      twas_genes <- str_replace(twas_genes, "^Glyma\\.", "GLYMA_")
    }
  }
  
  # Remove NAs after conversion
  gene_df <- gene_df %>% filter(!is.na(GENE))
  twas_genes <- twas_genes[!is.na(twas_genes)]
  
  background_genes <- unique(gene_df$GENE)
  
  
  # Run g:Profiler
  
  gost_res <- gost(
    query = twas_genes,
    organism = organism,
    sources = c("GO:BP","GO:MF","GO:CC","KEGG"),
    correction_method = "fdr",
    user_threshold = 0.05,
    domain_scope = "custom",
    custom_bg = background_genes
  )
  
  if(is.null(gost_res$result)) return(tibble())
  
  enrichment_results <- gost_res$result %>%
    as_tibble() %>%
    mutate(
      Species = species_name,
      fold_enrichment =
        (intersection_size / query_size) /
        (term_size / effective_domain_size),
      neg_log10_fdr = -log10(p_value)
    )
  
  return(enrichment_results)
}





