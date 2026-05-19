library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
infile <- str_remove(args[length(args) - 1], fixed('-'))
outfile <- str_remove(args[length(args)], fixed('-'))

summarisePeaks <- function(path) {
  files <- Sys.glob(path)
  peaks <- tibble()
  for(f in files) {
    base_file <- basename(f)
    
    df <- read_csv(f) %>%
      mutate(
        n_environments = if_else(str_detect(base_file, 'twoEnvs'), 2L, 3L),
        # Remove suffixes and add an underscore before peak_id to avoid concatenation issues
        trait = base_file,
        # convert columns to numeric (safe way)
        CHROM = as.numeric(CHROM),
        POS = as.numeric(POS),
        top_Pvalue = as.numeric(top_Pvalue),
        pStart = as.numeric(pStart),
        pStop = as.numeric(pStop),
        pLength = as.numeric(pLength),
        num_SNPs = as.numeric(num_SNPs),
        # Add trait info and keep peak_id numeric as is
        peak_id = as.character(peak_id)
      ) %>%
      mutate(
        transcript = word(trait, 1, sep = fixed('.')),       # before first dot
        gene_model = word(trait, 1, sep = fixed('_')),       # before first underscore
        param = word(trait, 2, sep = fixed('.')) %>%         # between first and second dot, or NA if none
                word(1, sep = fixed('_'))                     # first part before underscore if multiple parts
      )
    
    peaks <- bind_rows(peaks, df)
  }
  return(peaks)
}

write_csv(summarisePeaks(infile), outfile, quote = "needed")