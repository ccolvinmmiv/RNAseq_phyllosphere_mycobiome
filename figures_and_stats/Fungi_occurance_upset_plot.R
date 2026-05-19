library(tidyverse)
library(ComplexUpset)
library(patchwork)

setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")

# Helper to extract taxa names from file
get_taxa <- function(file) {
  read_csv(file) %>%
    select(!GENOTYPE) %>%
    colnames()
}

# Build presence/absence dataframe for one taxonomic level
make_upset_df <- function(maize_vec, sorghum_vec, soybean_vec, level_name) {
  all_taxa <- unique(c(maize_vec, sorghum_vec, soybean_vec))
  
  tibble(ID = all_taxa) %>%
    mutate(
      Maize = ID %in% maize_vec,
      Sorghum = ID %in% sorghum_vec,
      Soybean = ID %in% soybean_vec,
      TL = level_name
    )
}

# Load taxa per level
maize_classes   <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2020_maize_C.csv")
maize_orders    <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2020_maize_O.csv")
maize_families  <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2020_maize_F.csv")
maize_generas   <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2020_maize_G.csv")
maize_species   <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2020_maize_S.csv")

sorghum_classes  <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_C.csv")
sorghum_orders   <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_O.csv")
sorghum_families <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_F.csv")
sorghum_generas  <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_G.csv")
sorghum_species  <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_S.csv")

soybean_classes  <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_soybean_C.csv")
soybean_orders   <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_soybean_O.csv")
soybean_families <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_soybean_F.csv")
soybean_generas  <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_soybean_G.csv")
soybean_species  <- get_taxa("data/bracken_intermediates/GWAS_genos_winsor_pheno_file_soybean_S.csv")

# Build per-level dataframes
class_df   <- make_upset_df(maize_classes,  sorghum_classes,  soybean_classes,  "Class")
order_df   <- make_upset_df(maize_orders,   sorghum_orders,   soybean_orders,   "Order")
family_df  <- make_upset_df(maize_families, sorghum_families, soybean_families, "Family")
genus_df   <- make_upset_df(maize_generas,  sorghum_generas,  soybean_generas,  "Genus")
species_df <- make_upset_df(maize_species,  sorghum_species,  soybean_species,  "Species")

# Define intersections
all_intersections <- list(
  c("Maize"),
  c("Sorghum"),
  c("Soybean"),
  c("Maize", "Sorghum"),
  c("Maize", "Soybean"),
  c("Sorghum", "Soybean"),
  c("Maize", "Sorghum", "Soybean")
)

# Function to create upset plot
make_taxa_upset <- function(df, y_label) {
  ComplexUpset::upset(
    df,
    intersect = c("Maize", "Sorghum", "Soybean"),
    intersections = all_intersections,
    
    base_annotations = list(
      ' ' = intersection_size(
        fill = "grey10",      
        text = list(size = 3),
        theme = theme(
          axis.text = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 9, color = "black")
        ))),
    
    set_sizes = FALSE,      
    width_ratio = 0.2
  ) +
    labs(
      x = NULL,
      y = NULL
      ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      strip.text = element_text(size = 9, face = "bold"),
      axis.text.x = element_blank(),
      panel.border = element_blank()
    )
}

# Add TL label as a separate ggdraw object
add_tl_label <- function(plot, tl_name) {
  cowplot::ggdraw(plot) +
    cowplot::draw_label(tl_name, x = 0.525, y = 0.95, hjust = 1,
                        fontface = "bold", size = 9, color = "black")
}

# Wrap each plot with label
p_class   <- add_tl_label(make_taxa_upset(class_df), "Class")
p_order   <- add_tl_label(make_taxa_upset(order_df), "Order")
p_family  <- add_tl_label(make_taxa_upset(family_df), "Family")
p_genus   <- add_tl_label(make_taxa_upset(genus_df), "Genus")
p_species <- add_tl_label(make_taxa_upset(species_df), "Species")

# Combine vertically
combined_plot <- cowplot::plot_grid(
  p_class, p_order, p_family, p_genus, p_species,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1,1,1,1,1)
)

combined_plot

ggsave(
  "graphs/all_TL_taxa_upset_plots.png",
  plot = combined_plot,
  dpi = 300,
  bg = "white",
  units = "in",
  height = 9,
  width = 6.5
)
