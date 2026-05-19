library(rMVP)
library(tidyverse)

args <- commandArgs(TRUE)
vcf_path <- args[1]
output_prefix <- args[2]

MVP.Data(
  fileVCF = vcf_path,
  fileKin = TRUE,
  filePC  = TRUE,
  priority = "memory",
  maxLine = 10000,
  out = output_prefix
)