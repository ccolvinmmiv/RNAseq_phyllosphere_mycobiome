library(tidyverse)
library(vegan)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(svglite)
library(patchwork)



setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")


#################### Sorghum #####################

# Input file (Counts data, all samples with min read depth, no normalization, very rare taxa removed)
sorghum_family_counts_df <- read_csv("data/bracken_intermediates/filtered_total_and_rel_abund_table_ne_2021_sorgh_F.csv") %>%
  select(NAME, GENOTYPE, NUMBER_READS) %>%
  pivot_wider(names_from = GENOTYPE, values_from = NUMBER_READS) %>%
  mutate(NAME = str_split_i(NAME, "-", 1))

# Read data and build phyloseq object (COUNTS)
sorghum_family_otu_mat <- sorghum_family_counts_df %>%
  column_to_rownames("NAME") %>%
  as.matrix()

sorghum_family_ps <- phyloseq(
  otu_table(sorghum_family_otu_mat, taxa_are_rows = TRUE)
)


# Alpha diversity (COUNTS ONLY)
sorghum_family_alpha_wide <- estimate_richness(
  sorghum_family_ps,
  measures = c(
    "Observed",
    "Chao1",
    "Shannon",
    "Simpson",
    "InvSimpson",
    "Fisher"
  )
) %>%
  rownames_to_column("SampleID") %>%
  select(-starts_with("se."))   # drop SE columns

# Long format version
sorghum_family_alpha_long <- sorghum_family_alpha_wide %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Metric",
    values_to = "Value"
  )

# Sorghum alpha stats summary
sorghum_family_alpha_summary <- sorghum_family_alpha_wide %>%
  select(-SampleID) %>%
  summarise(across(everything(),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE),
                        min  = ~min(.x, na.rm = TRUE),
                        max  = ~max(.x, na.rm = TRUE)
                   )
  )) %>%
  pivot_longer(cols = everything(),
               names_to = c("Metric", ".value"),
               names_sep = "_") %>% 
  mutate(host = "Sorghum") %>%
  relocate(host)


# Relative abundance transformation
sorghum_family_ps_rel <- transform_sample_counts(
  sorghum_family_ps,
  function(x) x / sum(x)
)


# Bray–Curtis clustering of samples
sorghum_family_otu_rel_mat <- as(
  otu_table(sorghum_family_ps_rel),
  "matrix"
)

if (taxa_are_rows(sorghum_family_ps_rel)) {
  sorghum_family_otu_rel_mat <- t(sorghum_family_otu_rel_mat)
}

# Remove any empty samples
sorghum_family_otu_rel_mat <- sorghum_family_otu_rel_mat[
  rowSums(sorghum_family_otu_rel_mat) > 0,
]

sorghum_family_sample_order <- hclust(
  vegdist(sorghum_family_otu_rel_mat, method = "bray")
)$order %>%
  { rownames(sorghum_family_otu_rel_mat)[.] }


# Prepare data for community bar plot
sorghum_family_melt <- psmelt(sorghum_family_ps_rel)


sorghum_family_top_taxa <- sorghum_family_melt %>%
  group_by(OTU) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  slice_max(mean_abund, n = 9) %>%
  pull(OTU)

# Collapse rare taxa into "Other" and set sample order
sorghum_family_melt <- sorghum_family_melt %>%
  mutate(
    OTU = if_else(OTU %in% sorghum_family_top_taxa, OTU, "Other"),
    Sample = factor(Sample, levels = sorghum_family_sample_order)
  )

# Compute mean abundance ordering (descending)
sorghum_family_taxa_order <- sorghum_family_melt %>%
  group_by(OTU) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  filter(OTU != "Other") %>%
  arrange(desc(mean_abundance)) %>%
  pull(OTU)
# Explicit stacking order:
sorghum_family_melt$OTU <- factor(
  sorghum_family_melt$OTU,
  levels = c(
    sort(setdiff(unique(sorghum_family_melt$OTU), "Other")),
    "Other"
  )
)

# Community composition bar plot
sorghum_taxa_colors <- c(
  "#ED665DFF",   
  "#FF9E4AFF",  
  "#CDCC5DFF",   
  "#6DCCDAFF",  
  "#729ECEFF",  
  "#ED97CAFF",
  "#AD8BC9FF",  
  "#A8786EFF",  
  "#67BF5CFF",  
  "grey70"      
)

sorghum_otu_levels <- levels(sorghum_family_melt$OTU)
sorghum_taxa_colors_named <- setNames(sorghum_taxa_colors, sorghum_otu_levels)

sorghum_family_barplot <- ggplot(sorghum_family_melt, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = sorghum_taxa_colors_named) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title = NULL,
    x = NULL,
    y = "Relative Abundance",
    fill = NULL
  ) +
  annotate("text", x = Inf, y = 1, label = "Sorghum", hjust = 1.1, vjust = 1.3, size = 3.3, fontface = "bold") +
  guides(fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.5), nrow = 3, title.position = "top")) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0.25, 0.1, 0.25, 0.1), "lines"),
    legend.position = c(0.52, 0.12),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9, color = "black", margin = margin(b = 1)),
    legend.key.size = unit(0.6, "lines"),
    legend.background = element_rect(fill = NA, color = NA),
    plot.title = element_blank(),
    legend.box.spacing = unit(1, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.spacing.y = unit(1, "pt")
  )

#sorghum_family_barplot

# Objects created (for reference)

# sorghum_family_ps              -> counts phyloseq object
# sorghum_family_ps_rel          -> relative abundance phyloseq object
# sorghum_family_alpha_wide      -> samples x diversity metrics
# sorghum_family_alpha_long      -> long-format diversity table
# sorghum_family_sample_order    -> Bray–Curtis sample order
# sorghum_family_barplot         -> community composition figure












#################### Maize #####################

# Input file (Counts data, all samples with min read depth, no normalization, very rare taxa removed)
Maize_family_counts_df <- read_csv("data/bracken_intermediates/filtered_total_and_rel_abund_table_ne_2020_maize_F.csv") %>%
  select(NAME, GENOTYPE, NUMBER_READS) %>%
  pivot_wider(names_from = GENOTYPE, values_from = NUMBER_READS) %>%
  mutate(NAME = str_split_i(NAME, "-", 1))

# Read data and build phyloseq object (COUNTS)
Maize_family_otu_mat <- Maize_family_counts_df %>%
  column_to_rownames("NAME") %>%
  as.matrix()

Maize_family_ps <- phyloseq(
  otu_table(Maize_family_otu_mat, taxa_are_rows = TRUE)
)


# Alpha diversity (COUNTS ONLY)
Maize_family_alpha_wide <- estimate_richness(
  Maize_family_ps,
  measures = c(
    "Observed",
    "Chao1",
    "Shannon",
    "Simpson",
    "InvSimpson",
    "Fisher"
  )
) %>%
  rownames_to_column("SampleID") %>%
  select(-starts_with("se."))   # drop SE columns

# Long format version
Maize_family_alpha_long <- Maize_family_alpha_wide %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Metric",
    values_to = "Value"
  )

# Maize alpha stats summary
Maize_family_alpha_summary <- Maize_family_alpha_wide %>%
  select(-SampleID) %>%
  summarise(across(everything(),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE),
                        min  = ~min(.x, na.rm = TRUE),
                        max  = ~max(.x, na.rm = TRUE)
                   )
  )) %>%
  pivot_longer(cols = everything(),
               names_to = c("Metric", ".value"),
               names_sep = "_") %>% 
  mutate(host = "Maize") %>%
  relocate(host)

# Relative abundance transformation
Maize_family_ps_rel <- transform_sample_counts(
  Maize_family_ps,
  function(x) x / sum(x)
)


# Bray–Curtis clustering of samples
Maize_family_otu_rel_mat <- as(
  otu_table(Maize_family_ps_rel),
  "matrix"
)

if (taxa_are_rows(Maize_family_ps_rel)) {
  Maize_family_otu_rel_mat <- t(Maize_family_otu_rel_mat)
}

# Remove any empty samples
Maize_family_otu_rel_mat <- Maize_family_otu_rel_mat[
  rowSums(Maize_family_otu_rel_mat) > 0,
]

Maize_family_sample_order <- hclust(
  vegdist(Maize_family_otu_rel_mat, method = "bray")
)$order %>%
  { rownames(Maize_family_otu_rel_mat)[.] }


# Prepare data for community bar plot
Maize_family_melt <- psmelt(Maize_family_ps_rel)


Maize_family_top_taxa <- Maize_family_melt %>%
  group_by(OTU) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  slice_max(mean_abund, n = 9) %>%
  pull(OTU)

# Collapse rare taxa into "Other" and set sample order
Maize_family_melt <- Maize_family_melt %>%
  mutate(
    OTU = if_else(OTU %in% Maize_family_top_taxa, OTU, "Other"),
    Sample = factor(Sample, levels = Maize_family_sample_order)
  )

# Compute mean abundance ordering (descending)
Maize_family_taxa_order <- Maize_family_melt %>%
  group_by(OTU) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  filter(OTU != "Other") %>%
  arrange(desc(mean_abundance)) %>%
  pull(OTU)
# Explicit stacking order:
Maize_family_melt$OTU <- factor(
  Maize_family_melt$OTU,
  levels = c(
    sort(setdiff(unique(Maize_family_melt$OTU), "Other")),
    "Other"
  )
)

# Community composition bar plot
Maize_taxa_colors <- c(
  "#ED665DFF",   
  "#FF9E4AFF",  
  "#CDCC5DFF",   
  "#6DCCDAFF",  
  "#729ECEFF",  
  "#AD8BC9FF",  
  "#67BF5CFF",
  "#ED97CAFF",  
  "#A8786EFF",  
  "grey70"      
)

Maize_otu_levels <- levels(Maize_family_melt$OTU)
Maize_taxa_colors_named <- setNames(Maize_taxa_colors, Maize_otu_levels)

Maize_family_barplot <- ggplot(Maize_family_melt, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = Maize_taxa_colors_named) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title = NULL,
    x = NULL,
    y = "Relative Abundance",
    fill = NULL
  ) +
  annotate("text", x = Inf, y = 1, label = "Maize", hjust = 1.1, vjust = 1.3, size = 3.3, fontface = "bold") +
  guides(fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.5), nrow = 3, title.position = "top")) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 9),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0.25, 0.1, 0.25, 0.1), "lines"),
    legend.position = c(0.52, 0.12),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9, color = "black", margin = margin(b = 1)),
    legend.key.size = unit(0.6, "lines"),
    legend.background = element_rect(fill = NA, color = NA),
    plot.title = element_blank(),
    legend.box.spacing = unit(1, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.spacing.y = unit(1, "pt")
  )

#Maize_family_barplot










#################### Soybean #####################

# Input file (Counts data, all samples with min read depth, no normalization, very rare taxa removed)
Soybean_family_counts_df <- read_csv("data/bracken_intermediates/filtered_total_and_rel_abund_table_soybean_F.csv") %>%
  select(NAME, GENOTYPE, NUMBER_READS) %>%
  pivot_wider(names_from = GENOTYPE, values_from = NUMBER_READS) %>%
  mutate(NAME = str_split_i(NAME, "-", 1))

# Read data and build phyloseq object (COUNTS)
Soybean_family_otu_mat <- Soybean_family_counts_df %>%
  column_to_rownames("NAME") %>%
  as.matrix()

Soybean_family_ps <- phyloseq(
  otu_table(Soybean_family_otu_mat, taxa_are_rows = TRUE)
)


# Alpha diversity (COUNTS ONLY)
Soybean_family_alpha_wide <- estimate_richness(
  Soybean_family_ps,
  measures = c(
    "Observed",
    "Chao1",
    "Shannon",
    "Simpson",
    "InvSimpson",
    "Fisher"
  )
) %>%
  rownames_to_column("SampleID") %>%
  select(-starts_with("se."))   # drop SE columns

# Long format version
Soybean_family_alpha_long <- Soybean_family_alpha_wide %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Metric",
    values_to = "Value"
  )


# Soybean alpha stats summary
Soybean_family_alpha_summary <- Soybean_family_alpha_wide %>%
  select(-SampleID) %>%
  summarise(across(everything(),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE),
                        min  = ~min(.x, na.rm = TRUE),
                        max  = ~max(.x, na.rm = TRUE)
                   )
  )) %>%
  pivot_longer(cols = everything(),
               names_to = c("Metric", ".value"),
               names_sep = "_") %>% 
  mutate(host = "Soybean") %>%
  relocate(host)

# Relative abundance transformation
Soybean_family_ps_rel <- transform_sample_counts(
  Soybean_family_ps,
  function(x) x / sum(x)
)


# Bray–Curtis clustering of samples
Soybean_family_otu_rel_mat <- as(
  otu_table(Soybean_family_ps_rel),
  "matrix"
)

if (taxa_are_rows(Soybean_family_ps_rel)) {
  Soybean_family_otu_rel_mat <- t(Soybean_family_otu_rel_mat)
}

# Remove any empty samples
Soybean_family_otu_rel_mat <- Soybean_family_otu_rel_mat[
  rowSums(Soybean_family_otu_rel_mat) > 0,
]

Soybean_family_sample_order <- hclust(
  vegdist(Soybean_family_otu_rel_mat, method = "bray")
)$order %>%
  { rownames(Soybean_family_otu_rel_mat)[.] }


# Prepare data for community bar plot
Soybean_family_melt <- psmelt(Soybean_family_ps_rel)


Soybean_family_top_taxa <- Soybean_family_melt %>%
  group_by(OTU) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  slice_max(mean_abund, n = 9) %>%
  pull(OTU)

# Collapse rare taxa into "Other" and set sample order
Soybean_family_melt <- Soybean_family_melt %>%
  mutate(
    OTU = if_else(OTU %in% Soybean_family_top_taxa, OTU, "Other"),
    Sample = factor(Sample, levels = Soybean_family_sample_order)
  )

# Compute mean abundance ordering (descending)
Soybean_family_taxa_order <- Soybean_family_melt %>%
  group_by(OTU) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  filter(OTU != "Other") %>%
  arrange(desc(mean_abundance)) %>%
  pull(OTU)
# Explicit stacking order:
Soybean_family_melt$OTU <- factor(
  Soybean_family_melt$OTU,
  levels = c(
    sort(setdiff(unique(Soybean_family_melt$OTU), "Other")),
    "Other"
  )
)

# Community composition bar plot
Soybean_taxa_colors <- c(
  "#ED665DFF",   
  "#FF9E4AFF",  
  "#CDCC5DFF",   
  "#6DCCDAFF",  
  "#729ECEFF",  
  "#AD8BC9FF",  
  "#67BF5CFF",
  "#ED97CAFF",  
  "#A8786EFF",  
  "grey70"      
)

Soybean_otu_levels <- levels(Soybean_family_melt$OTU)
Soybean_taxa_colors_named <- setNames(Soybean_taxa_colors, Soybean_otu_levels)

Soybean_family_barplot <- ggplot(Soybean_family_melt, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = Soybean_taxa_colors_named) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title = NULL,
    x = NULL,
    y = "Relative Abundance",
    fill = NULL
  ) +
  annotate("text", x = Inf, y = 1, label = "Soybean", hjust = 1.1, vjust = 1.3, size = 3.3, fontface = "bold") +
  guides(fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.5), nrow = 3, title.position = "top")) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0.25, 0.1, 0.25, 0.1), "lines"),
    legend.position = c(0.52, 0.12),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9, color = "black", margin = margin(b = 1)),
    legend.key.size = unit(0.6, "lines"),
    legend.background = element_rect(fill = NA, color = NA),
    plot.title = element_blank(),
    legend.box.spacing = unit(1, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.spacing.y = unit(1, "pt")
  )

#Soybean_family_barplot

# Objects created (for reference)

# Soybean_family_ps              -> counts phyloseq object
# Soybean_family_ps_rel          -> relative abundance phyloseq object
# Soybean_family_alpha_wide      -> samples x diversity metrics
# Soybean_family_alpha_long      -> long-format diversity table
# Soybean_family_sample_order    -> Bray–Curtis sample order
# Soybean_family_barplot         -> community composition figure

sorghum_otu_levels
Maize_otu_levels
Soybean_otu_levels

# Combining final community barplots
stacked_3_species_barplot <- (sorghum_family_barplot + labs(y = NULL)) / Maize_family_barplot / (Soybean_family_barplot + labs(y = NULL))
stacked_3_species_barplot 


ggsave("graphs/top_families_all_3_species_combined_5_high.png", plot = stacked_3_species_barplot, dpi = 300, units = "in", width = 6.5, height = 6, bg = "white")


# Combining final alpha diversity stats tables

final_combined_alpha_div_stats_summary <- bind_rows(sorghum_family_alpha_summary, Maize_family_alpha_summary, Soybean_family_alpha_summary)
write_csv(final_combined_alpha_div_stats_summary, "graphs/supplemental_tables/final_combined_alpha_div_summary_table.csv")

  



