library(tidyverse)
library(devtools)
library(SpiecEasi)
library(igraph)
library(Matrix)
library(phyloseq)
library(vegan)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(viridis)
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

sorghum_family_SE_output <- spiec.easi(sorghum_family_ps, 
                                    method = 'mb', 
                                    lambda.min.ratio = 1e-2, 
                                    nlambda = 30, 
                                    pulsar.params = list(rep.num = 100, ncores = 1))

# Stability scores
getStability(sorghum_family_SE_output)
# Chosen number of edges
getOptLambda(sorghum_family_SE_output)


#Plotting
sorghum_taxa_names <- taxa_names(sorghum_family_ps)
sorghum_adj <- getRefit(sorghum_family_SE_output)
sorghum_net <- adj2igraph(sorghum_adj)

V(sorghum_net)$family <- sorghum_taxa_names
V(sorghum_net)$degree <- degree(sorghum_net)

# Convert to tidygraph and annotate top 10 nodes
sorghum_graph_tbl <- as_tbl_graph(sorghum_net) %>%
  mutate(
    name = as.character(family),
    degree = centrality_degree(),
    label = if_else(rank(-degree) <= 10, name, NA_character_)  # only top 10 get labels
  )

sorghum_graph_tbl_filtered <- sorghum_graph_tbl %>%
  filter(degree > 0)

# Plot with labels on top 10 nodes by degree

sorghum_top10_nodes <- sorghum_graph_tbl_filtered %>%
  arrange(desc(degree)) %>%
  slice(1:10) %>%
  pull(name)

# Assign color only to top 10 nodes; rest are grey
sorghum_layout_coords <- create_layout(sorghum_graph_tbl_filtered, layout = "kk")
sorghum_layout_coords <- sorghum_layout_coords %>%
  mutate(
    is_top10 = name %in% sorghum_top10_nodes,
    node_group = ifelse(is_top10, name, "Other"),  # For coloring nodes
    label = ifelse(is_top10, name, NA_character_)  # Only label top 10
  )

# Build color palette: top 10 + "Other" in grey
sorghum_viridis_colors <- viridis(length(sorghum_top10_nodes), option = "D", end = 0.5)
names(sorghum_viridis_colors) <- sorghum_top10_nodes
sorghum_palette <- c(sorghum_viridis_colors, Other = "grey70")

scale_color_viridis_c(option = "viridis", begin = 0, end = 0.5, direction = -1)


# Get top 10 node names by degree
sorghum_top10_names <- sorghum_graph_tbl %>%
  as_tibble() %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10) %>%
  pull(name)

sorghum_top10_colors <- viridis(10, option = "viridis", begin = 0, end = 0.6)

sorghum_graph_tbl_filtered <- sorghum_graph_tbl %>%
  filter(degree > 0)


sorghum_graph_tbl_filtered <- sorghum_graph_tbl_filtered %>%
  activate(nodes) %>%
  mutate(
    node_group = if_else(name %in% sorghum_top10_names, name, "Other"),
    node_group = factor(node_group, levels = c(sorghum_top10_names, "Other"))
  ) %>%
  mutate(
    label = if_else(name %in% sorghum_top10_names, name, NA_character_)
  )




sorghum_network <- ggraph(sorghum_layout_coords) +
  geom_edge_link(alpha = 0.5, color = "grey70") +
  geom_node_point(aes(color = node_group, size = degree), show.legend = c(color = FALSE, size = TRUE)) +
  geom_text_repel(
    aes(x = x, y = y, label = label, color = node_group),
    size = 3.15,
    box.padding = 1,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = sorghum_palette) +
  scale_size_continuous(range = c(1, 8), name = "Edges") +
  theme_void() +
  theme(
    legend.position = "none",
    legend.justification = c(0, 0),
    plot.margin = margin(0, 0, 0, 0), 
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, color = "black"),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  annotate("text", x = Inf, y = Inf, label = "Sorghum", hjust = 1.1, vjust = 1.3, size = 3.3, fontface = "bold")


sorghum_top10_taxa_df <- sorghum_graph_tbl %>%
  as_tibble() %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10)

sorghum_top10_taxa_df$name <- factor(sorghum_top10_taxa_df$name, levels = sorghum_top10_taxa_df$name[order(-sorghum_top10_taxa_df$degree)])

sorghum_bar <- ggplot(sorghum_top10_taxa_df, aes(x = reorder(name, degree), y = degree, fill = name)) +
  geom_col() +
  scale_fill_manual(values = sorghum_palette) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
  coord_flip() + 
  labs(x = "", y = "", fill = "Taxa") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(size = 9, color = 'black'),
        axis.text.y = element_text(size = 9, color = 'black'),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        legend.text = element_text(size = 9)) 

combined_sorghum_network_plot <- sorghum_bar + sorghum_network + plot_layout(widths = c(1, 4), guides = "keep") &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  ) 
#combined_sorghum_network_plot

#ggsave("graphs/combined_sorghum_network_plot.png", plot = combined_sorghum_network_plot, dpi = 300, units = "in", width = 6.5, height = 3, bg = "white")




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

Maize_family_SE_output <- spiec.easi(Maize_family_ps, 
                                     method = 'mb', 
                                     lambda.min.ratio = 1e-2, 
                                     nlambda = 30, 
                                     pulsar.params = list(rep.num = 100, ncores = 1))

# Stability scores
getStability(Maize_family_SE_output)
# Chosen number of edges
getOptLambda(Maize_family_SE_output)


#Plotting
Maize_taxa_names <- taxa_names(Maize_family_ps)
Maize_adj <- getRefit(Maize_family_SE_output)
Maize_net <- adj2igraph(Maize_adj)

V(Maize_net)$family <- Maize_taxa_names
V(Maize_net)$degree <- degree(Maize_net)

# Convert to tidygraph and annotate top 10 nodes
Maize_graph_tbl <- as_tbl_graph(Maize_net) %>%
  mutate(
    name = as.character(family),
    degree = centrality_degree(),
    label = if_else(rank(-degree) <= 10, name, NA_character_)  # only top 10 get labels
  )

Maize_graph_tbl_filtered <- Maize_graph_tbl %>%
  filter(degree > 0)

# Plot with labels on top 10 nodes by degree

Maize_top10_nodes <- Maize_graph_tbl_filtered %>%
  arrange(desc(degree)) %>%
  slice(1:10) %>%
  pull(name)

# Assign color only to top 10 nodes; rest are grey
Maize_layout_coords <- create_layout(Maize_graph_tbl_filtered, layout = "kk")
Maize_layout_coords <- Maize_layout_coords %>%
  mutate(
    is_top10 = name %in% Maize_top10_nodes,
    node_group = ifelse(is_top10, name, "Other"),  # For coloring nodes
    label = ifelse(is_top10, name, NA_character_)  # Only label top 10
  )

# Build color palette: top 10 + "Other" in grey
Maize_viridis_colors <- viridis(length(Maize_top10_nodes), option = "D", end = 0.5)
names(Maize_viridis_colors) <- Maize_top10_nodes
Maize_palette <- c(Maize_viridis_colors, Other = "grey70")

scale_color_viridis_c(option = "viridis", begin = 0, end = 0.5, direction = -1)


# Get top 10 node names by degree
Maize_top10_names <- Maize_graph_tbl %>%
  as_tibble() %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10) %>%
  pull(name)

Maize_top10_colors <- viridis(10, option = "viridis", begin = 0, end = 0.6)

Maize_graph_tbl_filtered <- Maize_graph_tbl %>%
  filter(degree > 0)


Maize_graph_tbl_filtered <- Maize_graph_tbl_filtered %>%
  activate(nodes) %>%
  mutate(
    node_group = if_else(name %in% Maize_top10_names, name, "Other"),
    node_group = factor(node_group, levels = c(Maize_top10_names, "Other"))
  ) %>%
  mutate(
    label = if_else(name %in% Maize_top10_names, name, NA_character_)
  )




Maize_network <- ggraph(Maize_layout_coords) +
  geom_edge_link(alpha = 0.5, color = "grey70") +
  geom_node_point(aes(color = node_group, size = degree), show.legend = c(color = FALSE, size = TRUE)) +
  geom_text_repel(
    aes(x = x, y = y, label = label, color = node_group),
    size = 3.15,
    box.padding = 1,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = Maize_palette) +
  scale_size_continuous(range = c(1, 8), name = "Edges") +
  theme_void() +
  theme(
    legend.position = "none",
    legend.justification = c(0, 0),
    plot.margin = margin(0, 0, 0, 0), 
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, color = "black"),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  annotate("text", x = Inf, y = Inf, label = "Maize", hjust = 1.1, vjust = 1.3, size = 3.3, fontface = "bold")


Maize_top10_taxa_df <- Maize_graph_tbl %>%
  as_tibble() %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10)

Maize_top10_taxa_df$name <- factor(Maize_top10_taxa_df$name, levels = Maize_top10_taxa_df$name[order(-Maize_top10_taxa_df$degree)])

Maize_bar <- ggplot(Maize_top10_taxa_df, aes(x = reorder(name, degree), y = degree, fill = name)) +
  geom_col() +
  scale_fill_manual(values = Maize_palette) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
  coord_flip() + 
  labs(x = "", y = "", fill = "Taxa") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(size = 9, color = 'black'),
        axis.text.y = element_text(size = 9, color = 'black'),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        legend.text = element_text(size = 9)) 

combined_Maize_network_plot <- Maize_bar + Maize_network + plot_layout(widths = c(1, 4), guides = "keep") &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  ) 
#combined_Maize_network_plot

#ggsave("graphs/combined_Maize_network_plot.png", plot = combined_Maize_network_plot, dpi = 300, units = "in", width = 6.5, height = 3, bg = "white")






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

Soybean_family_SE_output <- spiec.easi(Soybean_family_ps, 
                                       method = 'mb', 
                                       lambda.min.ratio = 1e-2, 
                                       nlambda = 30, 
                                       pulsar.params = list(rep.num = 100, ncores = 1))

# Stability scores
getStability(Soybean_family_SE_output)
# Chosen number of edges
getOptLambda(Soybean_family_SE_output)


#Plotting
Soybean_taxa_names <- taxa_names(Soybean_family_ps)
Soybean_adj <- getRefit(Soybean_family_SE_output)
Soybean_net <- adj2igraph(Soybean_adj)

V(Soybean_net)$family <- Soybean_taxa_names
V(Soybean_net)$degree <- degree(Soybean_net)

# Convert to tidygraph and annotate top 10 nodes
Soybean_graph_tbl <- as_tbl_graph(Soybean_net) %>%
  mutate(
    name = as.character(family),
    degree = centrality_degree(),
    label = if_else(rank(-degree) <= 10, name, NA_character_)  # only top 10 get labels
  )

Soybean_graph_tbl_filtered <- Soybean_graph_tbl %>%
  filter(degree > 0)

# Plot with labels on top 10 nodes by degree

Soybean_top10_nodes <- Soybean_graph_tbl_filtered %>%
  arrange(desc(degree)) %>%
  slice(1:10) %>%
  pull(name)

# Assign color only to top 10 nodes; rest are grey
Soybean_layout_coords <- create_layout(Soybean_graph_tbl_filtered, layout = "kk")
Soybean_layout_coords <- Soybean_layout_coords %>%
  mutate(
    is_top10 = name %in% Soybean_top10_nodes,
    node_group = ifelse(is_top10, name, "Other"),  # For coloring nodes
    label = ifelse(is_top10, name, NA_character_)  # Only label top 10
  )

# Build color palette: top 10 + "Other" in grey
Soybean_viridis_colors <- viridis(length(Soybean_top10_nodes), option = "D", end = 0.5)
names(Soybean_viridis_colors) <- Soybean_top10_nodes
Soybean_palette <- c(Soybean_viridis_colors, Other = "grey70")

scale_color_viridis_c(option = "viridis", begin = 0, end = 0.5, direction = -1)


# Get top 10 node names by degree
Soybean_top10_names <- Soybean_graph_tbl %>%
  as_tibble() %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10) %>%
  pull(name)

Soybean_top10_colors <- viridis(10, option = "viridis", begin = 0, end = 0.6)

Soybean_graph_tbl_filtered <- Soybean_graph_tbl %>%
  filter(degree > 0)


Soybean_graph_tbl_filtered <- Soybean_graph_tbl_filtered %>%
  activate(nodes) %>%
  mutate(
    node_group = if_else(name %in% Soybean_top10_names, name, "Other"),
    node_group = factor(node_group, levels = c(Soybean_top10_names, "Other"))
  ) %>%
  mutate(
    label = if_else(name %in% Soybean_top10_names, name, NA_character_)
  )




Soybean_network <- ggraph(Soybean_layout_coords) +
  geom_edge_link(alpha = 0.5, color = "grey70") +
  geom_node_point(aes(color = node_group, size = degree), show.legend = c(color = FALSE, size = TRUE)) +
  geom_text_repel(
    aes(x = x, y = y, label = label, color = node_group),
    size = 3.15,
    box.padding = 1,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = Soybean_palette) +
  scale_size_continuous(range = c(1, 8), name = "Edges") +
  theme_void() +
  theme(
    legend.position = "none",
    legend.justification = c(0, 0),
    plot.margin = margin(0, 0, 0, 0), 
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, color = "black"),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  annotate("text", x = Inf, y = Inf, label = "Soybean", hjust = 1.1, vjust = 1.3, size = 3.3, fontface = "bold")


Soybean_top10_taxa_df <- Soybean_graph_tbl %>%
  as_tibble() %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10)

Soybean_top10_taxa_df$name <- factor(Soybean_top10_taxa_df$name, levels = Soybean_top10_taxa_df$name[order(-Soybean_top10_taxa_df$degree)])

Soybean_bar <- ggplot(Soybean_top10_taxa_df, aes(x = reorder(name, degree), y = degree, fill = name)) +
  geom_col() +
  scale_fill_manual(values = Soybean_palette) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
  coord_flip() + 
  labs(x = "", y = "Edges", fill = "Taxa") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(size = 9, color = 'black'),
        axis.text.y = element_text(size = 9, color = 'black'),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        legend.text = element_text(size = 9)) 

combined_Soybean_network_plot <- Soybean_bar + Soybean_network + plot_layout(widths = c(1, 4), guides = "keep") &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  ) 
#combined_Soybean_network_plot

#ggsave("graphs/combined_Soybean_network_plot.png", plot = combined_Soybean_network_plot, dpi = 300, units = "in", width = 6.5, height = 3, bg = "white")



stacked_3_species_family_network <- combined_sorghum_network_plot / combined_Maize_network_plot / combined_Soybean_network_plot +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16, face = "bold"))



ggsave("graphs/stacked_3_hosts_fungal_family_network_plot.svg", plot = stacked_3_species_family_network, dpi = 300, units = "in", width = 6.5, height = 9, bg = "white")












