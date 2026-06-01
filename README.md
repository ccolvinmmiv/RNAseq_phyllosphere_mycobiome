# RNAseq_phyllosphere_mycobiome

This repository contains scripts and analytical workflows used to characterize host genetic control of leaf-associated fungal communities across maize (*Zea mays*), sorghum (*Sorghum bicolor*), and soybean (*Glycine max*). The project integrates RNA-seq–derived fungal abundance profiles with host genomic and transcriptomic data to perform GWAS, TWAS, eQTL mapping, coexpression analysis, and gene regulatory network inference.

These scripts represent the original analysis code used in the study. File paths are system-specific and will need to be updated by users to match their local environment. All input datasets required to reproduce the analyses are available in the supplementary materials and public repositories referenced in the manuscript.

---
## Citation

If you use this code, please cite:

Colvin C, Chopra S. Population-scale transcriptomics reveals host genetic control of phyllosphere fungal communities. bioRxiv. 2026. https://doi.org/10.64898/2026.05.21.727028

---
# Data Sources

RNA-seq datasets were obtained from previously published field-grown diversity panels:

| Species | Panel | Location | Samples | Reference |
|----------|------|----------|---------|----------|
| Maize (*Zea mays*) | WiDiv | Lincoln, NE, USA | ~750 | Torres-Rodríguez et al., 2024, *The Plant Journal* |
| Sorghum (*Sorghum bicolor*) | SAP + SDP | Lincoln, NE, USA | ~822 | Mangal et al., 2025, *The Plant Journal* |
| Soybean (*Glycine max*) | Diversity panel | Sanya, China | ~622 | Li et al., 2024, *Plant Communications* |

Raw sequencing data:
- ENA: PRJEB67964 (maize), PRJEB83049 (sorghum)  
- GSA: CRA009979 (soybean)

---

## Genotype and Reference Data

Genotype data were available for maize and sorghum only:

- Maize genotypes: https://doi.org/10.5061/dryad.bnzs7h4f1  
- Sorghum genotypes: https://doi.org/10.6084/m9.figshare.27936195  

Reference taxonomy database:
- Kraken2 prebuilt database: https://benlangmead.github.io/aws-indexes/k2  

Soybean was not included in GWAS, eQTL, and heritability analyses due to lack of available genotype data.

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
    ├── GWAS/
    ├── TWAS/
    ├── coexpression/
    ├── read_processing/
    ├── tax_assignment/
    ├── orthofinder/
    └── heritability/
```


---

# Analysis Overview

Core analyses include:

- RNA-seq read processing and taxonomic classification (Trimmomatic, Kraken2, Bracken)
- Fungal abundance quantification and filtering
- Narrow-sense heritability estimation (LDAK)
- Genome-wide association studies (rMVP)
- Transcriptome-wide association studies (GAPIT)
- eQTL mapping
- Coexpression network analysis (WGCNA)
- Gene regulatory network inference (GENIE3)
- Orthology inference (OrthoFinder)
- Functional enrichment (gProfiler2)

---


# Software

Key tools used:

- Trimmomatic  
- Kraken2 + Bracken  
- Kallisto  
- rMVP  
- GAPIT  
- LDAK  
- WGCNA  
- GENIE3  
- OrthoFinder  
- gProfiler2  
- VCFtools  

R (≥4.5) and Python (≥3.8) were used for custom analyses.

---

# License

MIT License recommended for reusable analysis pipelines.

---

# Notes

- All input data required for reproduction are either publicly available or included in supplementary datasets.
- Intermediate files are not included unless necessary for reproducibility of downstream analyses.
- Species-specific workflows share a common analytical framework implemented across scripts.
