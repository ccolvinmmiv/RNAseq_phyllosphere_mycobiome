#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript parse_ldak_reml.R <reml_dir> <output_csv>")
}

reml_dir <- args[1]
out_csv  <- args[2]

reml_files <- list.files(
  reml_dir,
  pattern = "\\.reml$",
  full.names = TRUE
)

if (length(reml_files) == 0) {
  stop("No .reml files found")
}

parse_reml <- function(file) {

  lines <- readLines(file, warn = FALSE)

  get_val <- function(pattern, idx) {
    line <- grep(pattern, lines, value = TRUE)
    if (length(line) == 0) return(NA_real_)
    as.numeric(strsplit(trimws(line[1]), "\\s+")[[1]][idx])
  }

  her_line <- grep("^Her_All", lines, value = TRUE)

  tibble(
    Trait = basename(file) |>
      sub("_reml\\.reml$", "", x = _) |>
      sub("^maize_fungal_", "", x = _),

    h2 = if (length(her_line) > 0)
      as.numeric(strsplit(trimws(her_line), "\\s+")[[1]][2])
    else NA_real_,

    h2_SE = if (length(her_line) > 0)
      as.numeric(strsplit(trimws(her_line), "\\s+")[[1]][3])
    else NA_real_,

    LRT = get_val("^LRT_Stat", 2),
    LRT_P = get_val("^LRT_P", 2),
    Converged = any(grepl("^Converged YES", lines)),
    N = get_val("^With_Phenotypes", 2),
    File = basename(file)
  )
}

res <- bind_rows(lapply(reml_files, parse_reml))

write.csv(res, out_csv, row.names = FALSE)

cat("Parsed", nrow(res), "REML files\n")
cat("Output:", out_csv, "\n")
