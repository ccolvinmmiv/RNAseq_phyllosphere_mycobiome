library(GENIE3)
library(dplyr)
library(readr)

#sink("/storage/work/cfc5873/Ufo1_eQTL/logs/coexpression/GENIE3_no_ufo_predictor_insideR.log", split = TRUE)

# Load expression data (samples x genes, all numeric; logCPM TMM normalized)
expr <- readRDS("/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/input_data/maize_GENIE3_exp_mat.rds")
# colnames(expr) are V5 gene IDs
# Convert to a *pure numeric matrix*
expr <- as.matrix(expr)
storage.mode(expr) <- "double"

expr <- t(expr)

# Check:
cat(dim(expr), "\n")        # Should be # Genes as rows × # of genotypes/samples as columns
head(rownames(expr))        # Should be gene IDs
head(colnames(expr))        # Should be sample names

cat("Expression matrix loaded: ", ncol(expr), "samples,", nrow(expr), "genes\n")

# Load TF predictor list (already in V5 format)
tf_table <- read.csv("/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/input_data/Grassius_maize_V5_TF_list.csv")
# Assume TF IDs are in a column TF_ID
regulators <- tf_table$`gene.ID`

cat("Initial regulator count (TFs): ", length(regulators), "\n")

# Keep only regulators present in expression data
regulators <- intersect(regulators, rownames(expr))
cat("Regulators present in expression matrix: ", length(regulators), "\n")

# Save list for record
#write.table(regulators, "outputs/regulator_list_used.tsv",
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Run GENIE3
cat("Running GENIE3...\n")

weight_matrix <- GENIE3(expr,
                        regulators = regulators,
                        nCores = 24,
                        verbose = TRUE)

cat("GENIE3 finished. Saving results...\n")

saveRDS(weight_matrix,
        file = "/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/outputs/maize_GENIE3_weight_matrix_all_genes.rds")

# Convert to ranked edge list:
link_list <- getLinkList(weight_matrix)
write.table(link_list,
            "/storage/work/cfc5873/Ufo1_eQTL/wu_2022_data/scripts/coexpression/outputs/maize_GENIE3_links_all_genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("All done.\n")
