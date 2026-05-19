library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()   # multithreading

#sink("/storage/work/cfc5873/Ufo1_eQTL/logs/coexpression/WGCNA_log_insideR.log", split = TRUE)
# Load data: genes x samples matrix
datExpr <- readRDS("/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/input_data/maize_WGCNA_ready.rds")
  
# WGCNA expects samples x genes

# 1. Check good samples/genes
gsg <- goodSamplesGenes(datExpr, verbose=3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 2. Soft-threshold selection
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5)

power <- sft$powerEstimate  # choose automatically or manually

# 3. Network construction using blockwise modules
net <- blockwiseModules(
  datExpr,
  power = power,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/outputs/tom_block",
  maxBlockSize = 6000,
  nThreads = 24,  # match HPC cpus-per-task
  verbose = 5
)

# Save module results
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs <- net$MEs

save(moduleLabels, moduleColors, MEs, net, file="/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/outputs/maize_WGCNA_results.RData")

# Write gene-module table
annot <- data.frame(Gene=colnames(datExpr),
                    Module=moduleColors)
write.table(annot, "/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/coexpression/outputs/maize_WGCNA_gene_module_assignment.tsv", sep="\t", quote=FALSE, row.names=FALSE)
