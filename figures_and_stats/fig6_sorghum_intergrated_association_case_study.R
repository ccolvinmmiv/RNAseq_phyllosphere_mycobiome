library(tidyverse)
library(ggrepel)
library(ggtext)
library(ggpubr)
library(patchwork)
library(snpStats)
library(rstatix)
library(coin)

setwd("C:/Users/ccolv/OneDrive - The Pennsylvania State University/Chopra lab/Mycobiome")



sorghum_Clav_GWAS <- read_csv("data/case_studies/Clavicipitaceae_34397.MLM.csv") %>%
  mutate(Trait = "Clavicipitaceae Relative Abundance") %>%
  rename(PVAL = Clavicipitaceae_34397.MLM) %>%
  mutate(SNP_to_color = if_else(PVAL < 0.05 / 50150.02, "Clavicipitaceae Relative Abundance", NA_character_))


sorghum_Clav_gene_eQTL <- read_csv("data/case_studies/Sobic.004G214900.MLM.csv") %>%
  mutate(Trait = "Expression of Sobic.004G214900") %>%
  rename(PVAL = Sobic.004G214900.MLM) %>%
  mutate(SNP_to_color = if_else(PVAL < 0.05 / 50150.02, "Expression of Sobic.004G214900", NA_character_)) 


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

filtered_VCF_Chr04_58988275 <- mangal_2025_genos_vcf %>%
  filter(ID == "Chr04_58988275") 
population_expression_data <- read_csv("expression_data/raw_expression_data/ne_2021_sorgh_merged_gene_tpms.csv")


sorgh_used_samples_names <- read_csv("genotype_files/ne_sorgh_2021_representative_samples_by_TPMs.csv")
sorgh_used_genos_736 <- read_table("genotype_files/sorgh_final_filtered_736_GWAS_genotypes_list.txt", col_names = "GENOTYPE")


sorgh_genos_to_keep_for_coexp <- sorgh_used_samples_names %>%
  filter(GENOTYPE %in% sorgh_used_genos_736$GENOTYPE) %>%
  pull(SAMPLE)

sorghum_raw_count_data_correct_samples <- population_expression_data %>% #Need to subset to the same 736 genotypes/samples used for GWAS/TWAS
  select(TranscriptID, any_of(sorgh_genos_to_keep_for_coexp)) 





flipped_sorghum_raw_tpm_data <- sorghum_raw_count_data_correct_samples %>%
  pivot_longer(cols = -TranscriptID, names_to = "Sample", values_to = "expression_count") %>%
  pivot_wider(id_cols = Sample, names_from = TranscriptID, values_from = expression_count) %>%
  mutate(GENOTYPE = str_remove(Sample, "^4\\d{3}_")) %>%
  select(GENOTYPE, Sobic.004G214900.1) %>%
  relocate(GENOTYPE) %>%
  left_join(filtered_VCF_Chr04_58988275, join_by(GENOTYPE)) %>%
  filter(GENOTYPE_CODE %in% c("0|0", "1|1"))

# Above contains all genotypes, their allele at peak SNP (filtered to homozygous only), and their TPM expression data for Sobic.004G214900.1




#Plotting exp boxplots
expr_boxplot <- flipped_sorghum_raw_tpm_data %>%
  ggplot(aes(x = snp_type, y = Sobic.004G214900.1, fill = snp_type)) +
  
  geom_violin(alpha = 0.3, width = 0.8, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  
  scale_fill_manual(values = c("#E69F00", "#009E73")) +
  
  stat_compare_means(
    comparisons = list(c("C", "G")),
    method = "wilcox.test",
    label = "p.signif",
    size = 4,
    step.increase = 0.1
  ) +
  
  stat_summary(
    fun.data = function(x) {
      data.frame(y = 4050, label = paste0("n=", length(x)))
    },
    geom = "text",
    vjust = -0.5,
    size = 3
  ) +
  
  labs(
    x = "Peak SNP genotype",
    y = "Sobic.004G214900\nExpression (TPM)"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 9)
  )

expr_boxplot



# Plot rel abund boxplots
clav_abund_df <-  read_csv("data/processed_normalized_phenotype_files_for_GWAS/GWAS_genos_winsor_pheno_file_ne_2021_sorgh_COMBINED_ALL_LEVELS.csv") %>%
  select(GENOTYPE, Clavicipitaceae_34397)

clav_geno_abund <- clav_abund_df %>%
  left_join(filtered_VCF_Chr04_58988275, by = "GENOTYPE") %>%
  filter(GENOTYPE_CODE %in% c("0|0", "1|1"))

abund_boxplot <- clav_geno_abund %>%
  ggplot(aes(x = snp_type, y = Clavicipitaceae_34397, fill = snp_type)) +
  scale_y_continuous(limits = c(0, 0.03)) +
  geom_violin(alpha = 0.3, width = 0.8, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  
  scale_fill_manual(values = c("#E69F00", "#009E73")) +
  
  stat_compare_means(
    comparisons = list(c("C", "G")),
    method = "wilcox.test",
    label = "p.signif",
    size = 4,
    step.increase = 0.1,
    bracket.size = 0
  ) +
  
  stat_summary(
    fun.data = function(x) {
      data.frame(y = 0.026, label = paste0("n=", length(x)))
    },
    geom = "text",
    vjust = -0.5,
    size = 3
  ) +
  
  labs(
    x = "Peak SNP genotype",
    y = "Clavicipitaceae\nRelative abundance"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 9)
  )

abund_boxplot


################### PAV Heatmaps   ###########
Final_Sbicolor_gene_conversion_table <- read_csv("data/gene_data/Final_Sbicolor_gene_conversion_table.csv") %>%
  mutate(QUERY_GENOTYPE = str_replace(QUERY_GENOTYPE, "730", "BTx623"))

Final_Sbicolor_gene_conversion_table_hotspot_4_snps <- Final_Sbicolor_gene_conversion_table %>%
  left_join(filtered_VCF_Chr04_58988275, join_by(QUERY_GENOTYPE == GENOTYPE)) %>%
  filter(GENOTYPE_CODE %in% c("0|0", "1|1")) %>%
  filter(ref_chrom == 4) %>% 
  filter(ref_gene_center > 58500891) %>%
  filter(ref_gene_center < 59207864) 

ggplot(Final_Sbicolor_gene_conversion_table_hotspot_4_snps, aes(ref_gene_center, query_gene_center, color = snp_type)) +
  geom_point() +
  facet_wrap(~QUERY_GENOTYPE) +
  ggtitle("Chromosome_4_hotspot") +
  geom_vline(xintercept = 58988275)

presence_absence_matrix_1 <- Final_Sbicolor_gene_conversion_table_hotspot_4_snps %>%
  select(QUERY_GENOTYPE, Reference_Gene, snp_type) %>%  # Keep the column
  count(QUERY_GENOTYPE, Reference_Gene, snp_type) %>%
  mutate(present = 1) %>%
  select(-n) %>%
  pivot_wider(
    names_from = Reference_Gene,
    values_from = present,
    values_fill = list(present = 0)
  ) %>%
  pivot_longer(
    cols = -c(QUERY_GENOTYPE, snp_type),  # exclude both from pivoting
    names_to = "Reference_Gene",
    values_to = "Presence"
  ) %>%
  arrange(Reference_Gene)

presence_absence_matrix_1 <- presence_absence_matrix_1 %>%
  mutate(
    fill_color = case_when(
      Presence == 0 ~ "white",
      snp_type == "C" ~ "#FF5A61",
      snp_type == "G" ~ "#569BFF",
      TRUE ~ "grey"  # fallback for unexpected values
    )
  ) %>%
  arrange(snp_type) %>% 
  rename(snp_type = snp_type) %>%
  mutate(QUERY_GENOTYPE = factor(QUERY_GENOTYPE, levels = unique(QUERY_GENOTYPE)))

pam1 <- ggplot(presence_absence_matrix_1, aes(x = Reference_Gene, y = QUERY_GENOTYPE, fill = fill_color)) +
  geom_tile(color = "white") +
  scale_fill_identity(
    guide = "legend",
    labels = c("C", "G"),
    breaks = c("#FF5A61", "#569BFF"),
    name = "Chr04"
  ) +
  theme_set(theme_classic(base_size = 19)) +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 9, hjust = 1, colour = "black"),
    axis.title = element_blank(),
    axis.line = element_blank(),       # removes axis lines
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.box.margin = margin(t = 0, r = 0, b = -25, l = 0),
    legend.text = element_text(size = 9, colour = "black"),
    legend.title = element_text(size = 9, colour = "black"),
    legend.key.size = unit(0.15, "in"),
    legend.background = element_blank()
  )
pam1

###############  local region plot   ##########
peak_chr <- 4
peak_start <- 58663165   # BTx623 V5 coordinates
peak_end   <- 59133836
peak_snp   <- 58988275  # from GWAS peak
# TWAS significance threshold
FDR_cutoff <- 0.05

# FILTER GENES IN REGION
# ANNOTATE TWAS SIGNIFICANCE
sorghum_Clav_TWAS <- read_csv("data/case_studies/TWAS.CMLM_Clavicipitaceae_34397.csv") %>%
  mutate(Taxa = "Clavicipitaceae")

genes_region <- Final_Sbicolor_gene_conversion_table %>%
  filter(ref_chrom == peak_chr,
         start >= peak_start,
         end <= peak_end) %>%
  # Annotate TWAS significance
  left_join(
    sorghum_Clav_TWAS %>% 
      filter(FDR < 0.05) %>%
      select(SNP) %>% 
      mutate(TWAS_sig = TRUE),
    by = c("Reference_Gene" = "SNP")
  ) %>%
  mutate(TWAS_sig = if_else(is.na(TWAS_sig), FALSE, TRUE))



# Filter TWAS hits in region
twas_hits <- sorghum_Clav_TWAS %>%
  filter(Chromosome == 4,
         Position >= peak_start,
         Position <= peak_end,
         FDR < FDR_cutoff) %>%
  pull(SNP)

genes_region <- genes_region %>%
  mutate(TWAS_sig = Reference_Gene %in% twas_hits) %>%
  filter(QUERY_GENOTYPE == "BTx623")

# Assign y-track for plotting
genes_region <- genes_region %>%
  mutate(track = row_number())

# Keep SNPs in peak region
mangal_2025_genos_vcf_filtered <- mangal_2025_genos_vcf %>%
  filter(str_detect(ID, "Chr04"))

mangal_2025_genos_vcf_filtered1 <- mangal_2025_genos_vcf_filtered %>%
  # Split ID into CHROM and POS
  separate(ID, into = c("CHROM", "POS"), sep = "_", remove = FALSE) %>%
  mutate(
    CHROM = as.numeric(str_remove(CHROM, "Chr")),  # Chr01 -> 1
    POS = as.integer(POS)
  )

vcf_region <- mangal_2025_genos_vcf_filtered1 %>%
  filter(POS >= peak_start & POS <= peak_end) %>%
  select(ID, GENOTYPE, GENOTYPE_CODE)

geno_wide <- vcf_region %>%
  pivot_wider(names_from = GENOTYPE, values_from = GENOTYPE_CODE) %>%
  column_to_rownames("ID")

# Vectorized conversion using lookup
lookup <- c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)
# geno_numeric: rows = SNPs, columns = genotypes
geno_numeric <- matrix(
  lookup[as.matrix(geno_wide)],
  nrow = nrow(geno_wide),
  ncol = ncol(geno_wide),
  dimnames = list(rownames(geno_wide), colnames(geno_wide))
)

peak_snp_id <- paste0("Chr04_", peak_snp)
peak_index <- which(rownames(geno_numeric) == peak_snp_id)

# Calculate r² (squared correlation) between peak SNP and all other SNPs
r2_peak <- cor(t(geno_numeric), t(geno_numeric[peak_index, , drop = FALSE]), use = "pairwise.complete.obs")^2

# r2_peak is a vector with names
ld_df <- tibble(
  SNP = rownames(geno_numeric),
  r2_peak = as.numeric(r2_peak)
)

# Quick check
head(ld_df)

# Define region coordinates
region_start <- peak_start
region_end <- peak_end

# Filter genes in region
genes_in_region <- Final_Sbicolor_gene_conversion_table %>%
  filter(ref_chrom == 4,
         start >= region_start,
         end <= region_end) %>%
  left_join(
    sorghum_Clav_TWAS %>%
      filter(FDR < 0.05) %>%
      select(SNP, FDR),
    by = c("Reference_Gene" = "SNP")
  ) %>%
  mutate(TWAS_sig = if_else(!is.na(FDR), TRUE, FALSE))

# Merge LD with positions
ld_df <- ld_df %>%
  mutate(
    CHROM = as.integer(sub("Chr04_", "", SNP)),
    POS = as.integer(sub("Chr04_", "", SNP))
  )

ld_region <- ld_df %>%
  filter(SNP %in% rownames(geno_numeric)) %>%
  mutate(POS = as.integer(str_extract(SNP, "\\d+$"))) %>%  # extract numeric pos
  filter(POS >= peak_start & POS <= peak_end)

genes_region <- genes_region %>%
  mutate(label_face = ifelse(Query_Gene == "Sobic.004G214900", "bold", "plain"))

ld_gene_plot <- ggplot() +
  # LD points
  geom_point(data = ld_region, aes(x = POS, y = r2_peak), color = "#D55E00", size = 1.5) +
  geom_segment(
    aes(x = peak_snp, xend = peak_snp, y = 0, yend = 1),
    linetype = "dashed",
    color = "#0072B2"
  ) +
 # geom_vline(xintercept = 58954950, linetype = "dashed", color = "darkorange") + #eQTL peak SNP for Sobic.004G214900
  # Gene rectangles
  geom_rect(
    data = genes_region,
    aes(xmin = start, xmax = end, ymin = -0.05, ymax = 0),
    fill = ifelse(genes_region$TWAS_sig, "#E69F00", "#56B4E9"),
    color = ifelse(genes_region$Reference_Gene == "Sobic.004G214900",
                   "black", "grey30"),
    linewidth = ifelse(genes_region$Reference_Gene == "Sobic.004G214900", 0.8, 0.3)) +
  # Gene labels
  geom_text_repel(
    data = subset(genes_region, Reference_Gene == "Sobic.004G214900"),
    aes(x = (start + end)/2, y = 0, label = Reference_Gene),
    #direction = "y",
    nudge_x = 60000,
    nudge_y = 0.12,
    size = 3,
    fontface = "bold",
    segment.color = "black",
    min.segment.length = 0
  ) +
  scale_y_continuous(
    limits = c(-0.08, 1),
    breaks = seq(0, 1, by = 0.2),
    name = expression(r^2~with~peak~SNP),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    name = paste0("Chr0", peak_chr, " position (bp)"),
    limits = c(peak_start, peak_end)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(color = "black", size = 9),
    axis.text = element_text(color = "black", size = 9)
  )

ld_gene_plot



#################  Manhattan plots
# Getting values into correct format
sorghum_Clav_GWAS$CHROM <- as.character(sorghum_Clav_GWAS$CHROM)
sorghum_Clav_GWAS$CHROM <- gsub("^chr", "", sorghum_Clav_GWAS$CHROM, ignore.case = TRUE)
sorghum_Clav_GWAS$CHROM <- as.integer(sorghum_Clav_GWAS$CHROM)

#Transformation of p values
transformed_sorghum_Clav_GWAS <- sorghum_Clav_GWAS %>%
  mutate(PVAL = if_else(
    PVAL > 0,
    -log10(PVAL),
    -1*(log10(-1*PVAL))
  ))

##


# Getting values into correct format
sorghum_Clav_gene_eQTL$CHROM <- as.character(sorghum_Clav_gene_eQTL$CHROM)
sorghum_Clav_gene_eQTL$CHROM <- gsub("^chr", "", sorghum_Clav_gene_eQTL$CHROM, ignore.case = TRUE)
sorghum_Clav_gene_eQTL$CHROM <- as.integer(sorghum_Clav_gene_eQTL$CHROM)

#Transformation of p values
transformed_sorghum_Clav_gene_eQTL <- sorghum_Clav_gene_eQTL %>%
  mutate(PVAL = if_else(
    PVAL > 0,
    -log10(PVAL),
    -1*(log10(-1*PVAL))
  ))

##

single_manhattan_GWAS <- plot_manhattan_single(
  data = transformed_sorghum_Clav_GWAS,
  species = "sorghum",
  chr_col = CHROM,
  bp_col = POS,
  pval_col = PVAL,
  xaxis_lab = "",
  test_type = "Clavicipitaceae GWAS",
 # main = "Clavicipitaceae"
)

single_manhattan_eQTL <- plot_manhattan_single(
  data = transformed_sorghum_Clav_gene_eQTL,
  species = "sorghum",
  chr_col = CHROM,
  bp_col = POS,
  pval_col = PVAL,
  xaxis_lab = "Chromosome",
  test_type = "Sobic.004G214900 eQTL",
 # main = "Expression of Sobic.004G214900"
)






sorghum_Clav_TWAS <- read_csv("data/case_studies/TWAS.CMLM_Clavicipitaceae_34397.csv") %>%
  mutate(Taxa = "Clavicipitaceae")

chr_sizes <- sorghum_Clav_TWAS %>%
  group_by(Chromosome) %>%
  summarise(chr_len = max(Position, na.rm = TRUE)) %>%
  mutate(chr_start = cumsum(lag(chr_len, default = 0)))
sorghum_Clav_TWAS <- sorghum_Clav_TWAS %>%
  left_join(chr_sizes, by = "Chromosome") %>%
  mutate(BPcum = Position + chr_start)



Clav_start_01 <- chr_sizes %>%
  filter(Chromosome == 4) %>%
  pull(chr_start)
Clav_peak_01 <- Clav_start_01 + 58988275

Clavicipitaceae_TWAS_plot <- ggplot() +
  ggrastr::geom_point_rast(
    data = sorghum_Clav_TWAS %>% filter(FDR >= 0.05),
    aes(x = BPcum, y = -log10(FDR), color = factor(Chromosome)),
    size = 1,
    alpha = 0.6
  ) +
  ggrastr::geom_point_rast(size = 1.5, alpha = 0.6) +
  geom_segment(aes(x = Clav_peak_01, xend = Clav_peak_01, y = 0, yend = (1/2 * -log10(8.99638355331639e-81))), linetype = "21", color = "#67A9CF", alpha = 0.5, linewidth = 2) +
  scale_color_manual(
    values = rep(c("#001F5B", "#67A9CF"), 10),
    guide = "none"  # hide legend for chromosome coloring
  ) +
  ggrastr::geom_point_rast(
    data = sorghum_Clav_TWAS %>% filter(FDR < 0.05),
    aes(x = BPcum, y = -log10(FDR), fill = Taxa, color = Taxa),
    shape = 21,  
    size = 1.5,
    alpha = 0.6
  ) +
  ggrastr::geom_point_rast(size = 1.5, alpha = 0.6) +
  scale_fill_manual(
    values = c(
      `Clavicipitaceae` = "#001F5B"
    )
  ) +
  
  # Labels: use Taxa fill color
  geom_text_repel(
    data = sorghum_Clav_TWAS %>% filter(SNP == "Sobic.004G214900"),
    aes(x = BPcum, y = -log10(FDR), label = SNP),
    color = "black",
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    segment.size = 0.2,
    min.segment.length = 0
  ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_x_continuous(
    breaks = chr_sizes %>%
      mutate(center = chr_start + chr_len / 2) %>%
      pull(center),
    labels = chr_sizes$Chromosome
  ) +
  scale_y_continuous(
    name = expression(-log[10](FDR)),
    sec.axis = sec_axis(~ . / (1/2), name = expression(log[10](GWAS~p)))) +
  theme_bw(base_size = 12) +
  labs(
    x = "",
    fill = "Taxa",  # fill is used for sig SNPs
    color = "Taxa",
    title = ""
  ) +
  theme(
    axis.text.x = element_text(size = 9, color = 'black'),
    axis.text.y = element_text(size = 9, color = 'black'),
    axis.title.y = element_text(size = 9, color = 'black'),
    axis.title.x = element_text(size = 9, color = 'black'),
    legend.position = "hidden",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) 

Clavicipitaceae_TWAS_plot





### Combining plots
boxplots_right <- expr_boxplot / abund_boxplot

tight_theme <- theme(
  plot.margin = margin(0, 0, 0, 0),
  panel.spacing = unit(0, "pt")
)

manhattan_stack <- (single_manhattan_GWAS + tight_theme) /
  (Clavicipitaceae_TWAS_plot + tight_theme) /
  (single_manhattan_eQTL + tight_theme)

top_layer_plots <- (manhattan_stack | (boxplots_right + tight_theme)) +
  plot_layout(widths = c(2.5, 1))

final_plot <- (top_layer_plots / (ld_gene_plot + tight_theme)) +
  plot_layout(heights = c(3, 1)) +
  #plot_annotation(tag_levels = "A") &
  theme(
    #plot.tag = element_text(size = 16, face = "bold"),
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave("graphs/case_studies_examples/intergrated_association_sorghum_clav.svg", plot = final_plot, dpi = 300, height = 9, width = 6.5, units = "in")



#### Stats for MS ####

# EXPRESSION STATS

expr_stats <- flipped_sorghum_raw_tpm_data %>%
  group_by(snp_type) %>%
  summarise(
    n = n(),
    median_TPM = median(Sobic.004G214900.1, na.rm = TRUE),
    mean_TPM = mean(Sobic.004G214900.1, na.rm = TRUE),
    .groups = "drop"
  )

expr_wilcox <- wilcox.test(
  Sobic.004G214900.1 ~ snp_type,
  data = flipped_sorghum_raw_tpm_data
)
expr_wilcox$p.value

expr_fc <- expr_stats %>%
  summarise(
    log2FC = log2(median_TPM[snp_type == "G"] / median_TPM[snp_type == "C"]),
    fold_change = median_TPM[snp_type == "G"] / median_TPM[snp_type == "C"],
    diff_TPM = median_TPM[snp_type == "G"] - median_TPM[snp_type == "C"]
  )

# ABUNDANCE STATS

abund_stats <- clav_geno_abund %>%
  group_by(snp_type) %>%
  summarise(
    n = n(),
    median_abund = median(Clavicipitaceae_34397, na.rm = TRUE),
    mean_abund = mean(Clavicipitaceae_34397, na.rm = TRUE),
    .groups = "drop"
  )

abund_wilcox <- wilcox.test(
  Clavicipitaceae_34397 ~ snp_type,
  data = clav_geno_abund
)
abund_wilcox$p.value

abund_fc <- abund_stats %>%
  summarise(
    log2FC = log2(median_abund[snp_type == "G"] / median_abund[snp_type == "C"]),
    fold_change = median_abund[snp_type == "G"] / median_abund[snp_type == "C"],
    diff_abund = median_abund[snp_type == "G"] - median_abund[snp_type == "C"],
    percent_change = 100 * (median_abund[snp_type == "G"] - median_abund[snp_type == "C"]) /
      median_abund[snp_type == "C"]
  )

merged_df <- flipped_sorghum_raw_tpm_data %>%
  select(GENOTYPE, Sobic.004G214900.1) %>%
  inner_join(
    clav_geno_abund %>%
      select(GENOTYPE, Clavicipitaceae_34397),
    by = "GENOTYPE"
  )
cor.test(
  merged_df$Sobic.004G214900.1,
  merged_df$Clavicipitaceae_34397,
  method = "spearman"
)

# PRINT CLEAN OUTPUT
cat("\n--- Expression (Sobic.004G214900) ---\n")
print(expr_stats)
print(expr_fc)
print(expr_wilcox)

cat("\n--- Clavicipitaceae Abundance ---\n")
print(abund_stats)
print(abund_fc)
print(abund_wilcox)












## Function to make manhattan plot ##

plot_manhattan_single <- function(
    data,
    species = "maize",
    chr_col = CHROM,
    bp_col = POS,
    pval_col = PVAL,
    test_type = NULL,
    xaxis_lab = NULL,
    highlight_col = NULL,         # pass a column of TRUE/FALSE to highlight top SNPs
    main = NULL
) {
  
  # Clean chromosome numbers
  data <- data %>%
    mutate(
      CHROM_num = as.integer(gsub("chr", "", as.character({{ chr_col }})))
    ) %>%
    arrange(CHROM_num, {{ bp_col }})
  
  # Chromosome lengths
  if (species == "maize") {
    chromLength <- tibble(
      CHROM_num = 1:10,
      max_bp = c(
        308452471, 243675191, 238017767, 250330460, 226353449,
        181357234, 185808916, 182411202, 163004744, 152435371
      )
    )
  } else if (species == "sorghum") {
    chromLength <- tibble(
      CHROM_num = 1:10,
      max_bp = c(
        85112863, 79114963, 80873341, 71215609, 77058072,
        62713908, 68911884, 65779274, 63277606, 62870657
      )
    )
  } else stop("species must be 'maize' or 'sorghum'")
  
  # cumulative offsets
  offsets <- chromLength %>%
    arrange(CHROM_num) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(CHROM_num, bp_add)
  
  # add cumulative positions
  df <- data %>%
    left_join(offsets, by = "CHROM_num") %>%
    mutate(LOC = {{ bp_col }} + bp_add)
  
  # Colors (alternating blue shades)
  chr_colors <- rep(c("#001F5B", "#67A9CF"), length.out = length(unique(df$CHROM_num)))
  names(chr_colors) <- unique(df$CHROM_num)
  
  df$color_group <- factor(df$CHROM_num)
  
  # Significance threshold
  bonf <- -log10(0.05 / 50150.02)
  
  # Axis tick centers
  centers <- offsets %>%
    arrange(CHROM_num) %>%
    mutate(
      chr_start = bp_add,
      chr_end = lead(bp_add, default = max(df$LOC)),
      center = (chr_start + chr_end) / 2
    )
  x_pos <- max(df$LOC, na.rm = TRUE) * 0.98    # 98% of max x
  y_pos <- max(df[[rlang::as_name(ensym(pval_col))]], na.rm = TRUE) * 0.98  # 98% of max y
  # Base plot
  p <- ggplot(df, aes(x = LOC, y = {{ pval_col }})) +
    
    # SNP points — alternating blue by chromosome
    ggrastr::geom_point_rast(
      aes(color = color_group),
      size = 1.2, alpha = 0.75
    ) +
    
    # Optional highlighting of top SNPs
    { if (!is.null(highlight_col))
      geom_point(
        data = df %>% filter({{ highlight_col }} == TRUE),
        color = "red", size = 2, alpha = 0.9
      )
    } +
    
    geom_hline(yintercept = bonf, linetype = "dashed") +
    annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = test_type,
      hjust = 1, vjust = 1,   # anchor to top-right
      size = 3,
      fontface = "bold"
    ) +
    scale_color_manual(values = chr_colors, guide = "none") +
    
    scale_x_continuous(
      breaks = centers$center,
      labels = centers$CHROM_num
    ) +
    
    labs(
      title = main,
      x = xaxis_lab,
      y = expression(-log[10](p))
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, color = 'black'),
      axis.text.y = element_text(size = 9, color = 'black'),
      axis.title.y = element_text(size = 9, color = 'black'),
      axis.title.x = element_text(size = 9, color = 'black'),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  #ggsave("graphs/case_studies_examples/sorghum_Sobic.004G214900_eQTL.png", units = "in", dpi = 300, width = 10, height = 5)
  return(p)
}

