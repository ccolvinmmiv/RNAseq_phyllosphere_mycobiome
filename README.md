# RNAseq_phyllosphere_mycobiome

This repository contains scripts and analytical workflows used to characterize host genetic control of leaf-associated fungal communities across maize (*Zea mays*), sorghum (*Sorghum bicolor*), and soybean (*Glycine max*). The project integrates RNA-seq–derived fungal abundance profiles with host genomic, transcriptomic, GWAS, TWAS, eQTL, coexpression, and gene regulatory network analyses.

---

# Overview

We analyze publicly available field-collected RNA-seq datasets to:

- Profile fungal communities across host genotypes and environments
- Estimate narrow-sense heritability of fungal abundance traits
- Perform genome-wide association studies (GWAS)
- Conduct transcriptome-wide association studies (TWAS)
- Map expression quantitative trait loci (eQTLs)
- Identify coexpression modules and gene regulatory networks
- Infer cross-species orthologous fungal-associated genes

---

# Repository Structure

```
scripts/
├── data_processing/
│   ├── maize_process_bracken_outputs_loop.R
│   ├── sorghum_process_bracken_outputs_loop.R
│   ├── soybean_process_bracken_outputs_loop.R
│   ├── prep_all_data_for_TWAS.R
│   ├── prep_data_for_eQTL.R
│   ├── prep_data_for_hertability_calcs.R
│   ├── prep_fungal_supp_datasets.R
│   └── prepare_all_expression_data_for_coexpression_analysis.R
│
├── figures_and_stats/
│   ├── fig1_community_comp_graphs.R
│   ├── fig2_community_cooccurance_networks.R
│   ├── fig3_TWAS_ortho_GO_KEGG.R
│   ├── fig4_WGCNA_KEGG_and_GENIE3.R
│   ├── fig5_genome_wide_peak_distribution.R
│   ├── fig6_sorghum_intergrated_association_case_study.R
│   ├── fungal_read_counts_distribution.R
│   └── Fungi_occurance_upset_plot.R
│
├── HPC_scripts/
│   ├── GWAS/
│   ├── TWAS/
│   ├── coexpression/
│   ├── read_processing/
│   ├── tax_assignment/
│   ├── orthofinder/
│   └── heritability/
```

---

# Data Sources

We analyzed RNA-seq datasets from three previously published field-grown plant diversity panels:

| Species | Panel | Location | Samples | Primary publication |
|----------|------|----------|---------|----------------------|
| Maize (*Zea mays*) | WiDiv panel | Lincoln, NE, USA | ~750 | Torres-Rodríguez et al., 2024, *The Plant Journal*  |
| Sorghum (*Sorghum bicolor*) | SAP + SDP panels | Lincoln, NE, USA | ~822 | Mangal et al., 2025, *The Plant Journal* |
| Soybean (*Glycine max*) | Diversity panel | Sanya, China | ~622 | Li et al., 2024, *Plant Communications* |

Raw sequencing data are publicly available from:

- ENA: PRJEB83049 (sorghum), PRJEB67964 (maize)  
- GSA: CRA009979 (soybean)

Genotype data were available for maize and sorghum only; soybean was excluded from GWAS, eQTL, and heritability analyses.

---

# Supplementary Materials

## Supplementary Figures

**S1.** Distribution of fungal read counts per sample across datasets.  
**S2.** Overlap of fungal taxa detected across maize, sorghum, and soybean across taxonomic levels.  
**S3.** Top ≤10 KEGG terms enriched among TWAS-significant genes (maize, sorghum, soybean).  
**S4.** Top ≤10 GO Cellular Component terms enriched among TWAS-significant genes.  
**S5.** Top ≤10 GO Biological Process terms enriched among TWAS-significant genes.  
**S6.** Top ≤10 GO Molecular Function terms enriched among TWAS-significant genes.  
**S7.** Genome-wide distribution of GWAS and eQTL peaks across maize and sorghum genomes.

---

## Supplementary Tables

**S1.** Number of fungal taxa retained per dataset across taxonomic levels following filtering thresholds.  
**S2.** Number of genetic markers, effective markers (Meff), and significance thresholds used in GWAS.  
**S3.** Narrow-sense heritability (h²) estimates for fungal abundance traits in maize and sorghum.  
**S4.** Significant TWAS associations for fungal abundance traits across maize, sorghum, and soybean.  
**S5.** Significant GWAS peaks associated with fungal abundance traits in maize and sorghum.  
**S6.** eQTL peaks for TWAS-significant genes in maize and sorghum.  
**S7.** GWAS and eQTL peaks in sorghum chromosome 4 hotspot (Chr04: 57.99–59.96 Mb).  
**S8.** Fungal taxa associated with expression of Sobic.004G214900.  
**S9.** GWAS and eQTL peaks in sorghum chromosome 9 hotspot (Chr09: 61.62–63.28 Mb).  
**S10.** GWAS and eQTL peaks in sorghum chromosome 10 eQTL hotspot (Chr10: 58.88–60.74 Mb).

---

## Supplementary Datasets

**S1.** Maize fungal relative abundance matrix. Includes all samples passing the 2500 fungal-read threshold prior to winsorization and genotype-level collapsing for GWAS/TWAS analyses.  

**S2.** Sorghum fungal relative abundance matrix. Includes all samples passing the 2500 fungal-read threshold prior to winsorization and genotype-level collapsing for GWAS/TWAS analyses.  

**S3.** Soybean fungal relative abundance matrix. Includes all samples passing the 2500 fungal-read threshold prior to winsorization.  

**S4.** Orthogroup assignments for TWAS-significant genes across maize, sorghum, and soybean.

---

# Analysis Overview

Key analytical components include:

- RNA-seq processing (Trimmomatic)
- Taxonomic classification (Kraken2 + Bracken)
- Heritability estimation (LDAK REML)
- GWAS (rMVP MLM framework)
- TWAS (GAPIT cMLM framework)
- eQTL mapping (GWAS framework applied to expression traits)
- Coexpression analysis (WGCNA)
- Gene regulatory network inference (GENIE3)
- Orthology inference (OrthoFinder)

---

# License

Recommended: MIT License (for reusable pipelines and analysis scripts)

---

# Notes

- Intermediate files are not included unless required for reproducibility
- Representative scripts are provided for each major analysis stage
- Dataset-specific redundancies were removed where workflows are identical across species
