# **PRS/TRS in Risk Prediction of Nine Common Diseases in Vietnamese Population**

## Project Overview

This project implements and evaluates Polygenic Risk Score (PRS) methods for predicting the risk of nine common diseases in the Vietnamese population. The project uses various PRS calculation methods including PRSice2, LDpred2, PRS-CS, PRS-CSx and shaPRS (combine with LDpred2 and PRS-CS) on Vietnamese genomic data to assess disease risk prediction performance.

## Diseases Analyzed

The project analyzes the following nine diseases:

- Breast Cancer
- Colorectal Cancer
- Gastric Cancer
- Parkinson's Disease (PD)
- Chronic Kidney Disease (CKD)
- Osteoporosis
- Coronary Artery Disease (CAD)
- Hyperlipidemia
- Osteoarthritis

## Project Structure

```
VGP_10_diseases/
├── metadata/                    # Metadata files for the Vietnamese cohort
├── methods/                     # Implementation of various PRS methods
│   ├── LDpred2/
│   ├── PRScs/
│   ├── PRScsx/
│   ├── PRSice2/
│   └── shaprs/ 
├── pgs_catalog_results/         # Results from PGS Catalog analysis
├── reference_data/              # Reference data
│   ├── 1kgp/                    # 1000 Genomes Project data
│   └── ukbb/                    # UK Biobank data
├── results/                     # Results for each disease
├── sum_stats_data/              # Summary statistics for each disease
└── target_data/                 # Target genetic data for the Vietnamese population
│   ├── case_data/               # Genetic data of nine diseases
│   └── control_data/            # Genetic data of VN1K project (healthy individuals)
```

## Overall Workflow

1. **Data Preprocessing**:

   - **Summary Statistics Preprocessing**: Process GWAS summary statistics for different PRS methods (**preprocess_sumstats.ipynb**)
   - **Target Data Preprocessing**: Process VN genomic data (both case and control) including filtering, merging, and quality control (**preprocess_target_data.ipynb**)
2. **PRS and TRS Calculation**:

   - Multiple PRS methods are applied: PRSice2, LDpred2, PRS-CS, PRS-CSx and shaPRS (in **methods/** folder)
   - Each method is run with various parameters to assess in the Vietnamese population (results in **results/** folder)
   - PGS Catalog analysis (results in **pgs_catalog_results/** folder)
   - TRS...
3. **Evaluation**:

   - Performance assessment using metrics such as liability R² and deltaAUC (**predict.acc.PRS.R** and **predict.deltaAUC.R**)
   - Comparing PRS performance across different methods and diseases
   - Analysis of covariates including age, sex, 10 first principal components
4. **Visualization**:

   - Visualizing results using R and Python for better understanding of PRS performance
   - Generating plots for R² values, AUC, deltaAUC, and other relevant analysis (**visualization.results.ipynb**)

## Setup and Requirements

### Dependencies

- Python 3.11.9
- R 4.2.1
- CrossMap 0.7.0
- Bioinformatics tools: PLINK v1.9, bcftools v1.13
- Python packages: pandas, numpy, matplotlib, seaborn, scipy
- R packages: ggplot2, pROC, dplyr, bigsnpr

### Data Requirements

- Target Vietnamese genetic data (VN1K) in PLINK format
- GWAS summary statistics for each disease
- Reference genome data: 504 EAS samples in 1KGP
- HapMap3+.rds
- LD reference panels: 1kgp and/or ukbb

## Version History

- **v1.0 - May 29, 2025**: Initial release with key information of VGP.

## References

### Github Repositories:

* https://github.com/getian107/PRScs
* https://github.com/getian107/PRScsx
* https://github.com/mkelcb/shaprs

### Tutorials and Documentation:

* https://choishingwan.github.io/PRS-Tutorial/
* https://choishingwan.github.io/PRSice/
* https://privefl.github.io/bigsnpr/articles/LDpred2.html

### Research Papers:

* https://doi.org/10.1093/bioinformatics/btaa1029
* https://doi.org/10.1101/2021.12.10.21267272
* https://doi.org/10.1038/s41467-019-09718-5
* https://doi.org/10.1038/s41588-022-01054-7
