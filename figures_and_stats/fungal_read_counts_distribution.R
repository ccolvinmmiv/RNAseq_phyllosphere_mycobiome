library(tidyverse)
library(patchwork)

setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")



#
#exclude_genos <- c("B73HTRHM", "DK84QAB1", "HP72-11", "PHT69")
exclude_genos <- c()


maize_read_dist <- get_read_dist(
  "data/bracken_intermediates/filtered_total_and_rel_abund_table_ne_2020_maize_C.csv",
  "Maize",
  exclude_genos
)

soybean_read_dist <- get_read_dist(
  "data/bracken_intermediates/filtered_total_and_rel_abund_table_soybean_C.csv",
  "Soybean",
  exclude_genos
)

sorghum_read_dist <- get_read_dist(
  "data/bracken_intermediates/filtered_total_and_rel_abund_table_NE_2021_sorgh_C.csv",
  "Sorghum",
  exclude_genos
)

combined_df <- bind_rows(maize_read_dist, soybean_read_dist, sorghum_read_dist)


# Summary stats by dataset (environment)
read_summary <- combined_df %>%
  group_by(ENVIRONMENT) %>%
  summarise(
    n_samples = n(),
    mean_reads = mean(TOTAL_READS, na.rm = TRUE),
    median_reads = median(TOTAL_READS, na.rm = TRUE),
    sd_reads = sd(TOTAL_READS, na.rm = TRUE),
    min_reads = min(TOTAL_READS, na.rm = TRUE),
    max_reads = max(TOTAL_READS, na.rm = TRUE),
    q1_reads = quantile(TOTAL_READS, 0.25, na.rm = TRUE),
    q3_reads = quantile(TOTAL_READS, 0.75, na.rm = TRUE)
  )

read_summary


p <- plot_read_hist(combined_df)
p



ggsave(
  "graphs/all_datasets_combined_read_dist_histo.png",
  plot = p,
  dpi = 300,
  bg = "white",
  units = "in",
  height = 4.5,
  width = 6.5
)



#####

tax_ID_list <- read_csv("data/raw_bracken_outputs/full_NCBI_taxID_conversion_table.csv")
fungal_only_tax_ID_table <- tax_ID_list %>%
  filter(kingdom == "Fungi")
fungal_only_tax_ID_list <- fungal_only_tax_ID_table %>%
  pull(taxid)


res_maize <- summarize_fungal_reads(
  "data/raw_bracken_outputs/combined_ne_maize_2020_bracken_C.tsv",
  fungal_only_tax_ID_list
)

res_sorghum <- summarize_fungal_reads(
  "data/raw_bracken_outputs/combined_ne_sorgh_2021_bracken_C.tsv",
  fungal_only_tax_ID_list
)

res_soybean <- summarize_fungal_reads(
  "data/raw_bracken_outputs/combined_soybean_bracken_C.tsv",
  fungal_only_tax_ID_list
)


all_datasets_read_stats <- bind_rows(
  maize = res_maize$dataset_level,
  sorghum = res_sorghum$dataset_level,
  soybean = res_soybean$dataset_level,
  .id = "dataset"
)

all_datasets_read_stats





# Functions (call first)

#
# Function to read + summarize fungal read counts
get_read_dist <- function(file,
                          env_name,
                          exclude_genotypes = NULL,
                          min_reads = NULL) {
  
  df <- read_csv(file) %>%
    select(GENOTYPE, TOTAL_READS)
  
  # Remove unwanted genotypes
  if (!is.null(exclude_genotypes)) {
    df <- df %>% filter(!GENOTYPE %in% exclude_genotypes)
  }
  
  # Optional read filtering
  if (!is.null(min_reads)) {
    df <- df %>% filter(TOTAL_READS > min_reads)
  }
  
  # Summarize per genotype
  df %>%
    group_by(GENOTYPE) %>%
    summarise(TOTAL_READS = mean(TOTAL_READS, na.rm = TRUE), .groups = "drop") %>%
    mutate(ENVIRONMENT = env_name)
}

# Function to plot histogram
plot_read_hist <- function(df,
                           cutoff_display = 200000,
                           cutoff_line = 2500,
                           binwidth = 1000,
                           x_break_step = 25000) {
  
  # Count omitted samples
  omit_summary <- df %>%
    group_by(ENVIRONMENT) %>%
    summarise(num_omitted = sum(TOTAL_READS > cutoff_display), .groups = "drop")
  # Median per environment
  median_summary <- df %>%
    group_by(ENVIRONMENT) %>%
    summarise(median_reads = median(TOTAL_READS, na.rm = TRUE), .groups = "drop")
  
  
  df_labeled <- df %>%
    left_join(omit_summary, by = "ENVIRONMENT") %>%
    left_join(median_summary, by = "ENVIRONMENT") %>%
    mutate(
      facet_label = paste0(
        ENVIRONMENT,
        " (", num_omitted, " samples >", "200,000", " reads not shown)"
      )
    )
  
  ggplot(df_labeled, aes(x = TOTAL_READS)) +
    geom_histogram(binwidth = binwidth, color = "black", fill = "grey70") +
    geom_vline(xintercept = cutoff_line, color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_vline(aes(xintercept = median_reads), color = "blue", linetype = "dashed", linewidth = 0.8) +
    facet_wrap(vars(facet_label), ncol = 1, scales = "free_y") +
    scale_x_continuous(
      limits = c(0, cutoff_display),
      breaks = seq(0, cutoff_display, by = x_break_step),  # Finer ticks
      labels = function(x) format(x, scientific = FALSE, trim = TRUE),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Total Fungal Reads at Class Level",
      y = "Number of Samples"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      # Axis lines
      axis.line = element_line(color = "black"),
      
      # Tick marks
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      
      # Text styling
      strip.text = element_text(size = 9, color = "black"),
      strip.background = element_rect(fill = "white", color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      
      # Clean background
      panel.grid = element_blank(),
      
      plot.margin = margin(10, 20, 10, 10)
    )
}



summarize_fungal_reads <- function(bracken_file, fungal_tax_ids) {
  
  # Read full dataset
  df <- read_tsv(bracken_file)
  
  # Identify sample columns (those ending in "_num")
  sample_cols <- df %>%
    select(contains("_num")) %>%
    colnames()
  
  # TOTAL reads per sample (all taxa)
  total_reads <- df %>%
    select(all_of(sample_cols)) %>%
    summarise(across(everything(), ~sum(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "total_reads")
  
  # FUNGAL reads per sample
  fungal_reads <- df %>%
    filter(taxonomy_id %in% fungal_tax_ids) %>%
    select(all_of(sample_cols)) %>%
    summarise(across(everything(), ~sum(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "fungal_reads")
  
  # Merge + compute proportion
  summary_df <- total_reads %>%
    left_join(fungal_reads, by = "sample") %>%
    mutate(
      fungal_reads = replace_na(fungal_reads, 0),
      prop_fungal = fungal_reads / total_reads
    )
  
  # Dataset-level summary stats
  dataset_summary <- summary_df %>%
    summarise(
      total_fungal_reads = sum(fungal_reads),
      median_fungal_reads = median(fungal_reads),
      IQR_fungal_reads = IQR(fungal_reads),
      median_prop_fungal = median(prop_fungal),
      IQR_prop_fungal = IQR(prop_fungal)
    )
  
  return(list(
    sample_level = summary_df,
    dataset_level = dataset_summary
  ))
}
