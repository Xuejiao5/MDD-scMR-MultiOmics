# Multi-Omics Mendelian Randomization Identifies CKAP2 as a Causal Protective Factor in Major Depressive Disorder

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## Overview

This repository contains the complete analytical pipeline for a comprehensive multi-omics Mendelian randomization (MR) study investigating the causal relationship between gene expression, protein levels, and Major Depressive Disorder (MDD). The study integrates bulk tissue eQTL data, plasma pQTL data, and single-cell eQTL validation to identify causal protective factors.

### Key Findings
- **CKAP2** identified as a causal protective factor for MDD through multiple omics layers
- Cross-validated using bulk tissue (eQTLGen) and protein (deCODE) datasets
- Single-cell eQTL validation across brain cell types
- Brain imaging mediation analysis supporting causal pathways

## Study Design

### Analysis Components
1. **Bulk Tissue eQTL-MR**: Using eQTLGen consortium data
2. **Plasma pQTL-MR**: Using deCODE genetics protein QTL data
3. **Single-cell eQTL-MR**: Brain cell type-specific validation
4. **Brain Imaging Mediation**: IDP → MDD pathway analysis
5. **Sensitivity Analyses**: Reverse MR, pleiotropy tests, heterogeneity assessment

### Statistical Methods
- **Two-Sample MR**: IVW, MR-Egger, Weighted Median, Mode-based methods
- **GSMR**: Generalized Summary-data-based MR with HEIDI outlier detection
- **MR-PRESSO**: Pleiotropy RESidual Sum and Outlier detection
- **Steiger Filtering**: Directional testing for causal inference
- **Mediation Analysis**: Sobel test with bootstrap confidence intervals

## Installation

### Prerequisites
- R ≥ 4.0.0
- RStudio (recommended)
- Git

### Quick Installation
```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yourusername/MDD-scMR-MultiOmics")

# Or clone and install locally
git clone https://github.com/yourusername/MDD-scMR-MultiOmics.git
cd MDD-scMR-MultiOmics
```

### Using renv for Reproducibility
```r
# Install renv if not already installed
if (!require("renv")) install.packages("renv")

# Restore the exact package environment
renv::restore()
```

### Manual Package Installation
```r
# Core MR packages
install.packages(c("TwoSampleMR", "gsmr2", "data.table"))

# Visualization and general packages  
install.packages(c("ggplot2", "dplyr", "tidyverse", "ggrepel", "patchwork"))

# Additional packages for specific analyses
install.packages(c("foreach", "progress", "boot", "logging"))

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("MendelianRandomization")
```

## Quick Start

### Basic Usage Example
```r
library(MddScMrMultiOmics)

# Load example data
data(example_exposure_data)
data(example_outcome_data)

# Run basic MR analysis
mr_results <- perform_mr_analysis(
  exposure_data = example_exposure_data,
  outcome_data = example_outcome_data,
  clump_r2 = 0.001,
  pval_threshold = 5e-8
)

# Generate summary plots
plot_mr_results(mr_results)
```

### Complete Analysis Workflow
```r
# 1. eQTL-MR Analysis
source("scripts/01_eqtl_mr_analysis.R")

# 2. pQTL-MR Analysis  
source("scripts/02_pqtl_mr_analysis.R")

# 3. Single-cell validation
source("scripts/03_single_cell_eqtl_mr.R")

# 4. Mediation analysis
source("scripts/04_mediation_analysis.R")

# 5. Sensitivity analyses
source("scripts/05_reverse_mr_sensitivity.R")

# 6. Generate publication figures
source("scripts/06_generate_figures.R")
```

## Data Requirements

### Input Data Formats
All GWAS summary statistics should be in standard format with the following columns:

**Exposure Data:**
- `SNP`: SNP identifier (rsID)
- `beta`: Effect size
- `se`: Standard error
- `pval`: P-value
- `effect_allele`: Effect allele
- `other_allele`: Other allele
- `eaf`: Effect allele frequency
- `samplesize`: Sample size

**Outcome Data:**
- Same format as exposure data
- Must include overlapping SNPs with exposure

### Data Sources Used
- **eQTL Data**: eQTLGen Consortium
- **pQTL Data**: deCODE genetics
- **MDD GWAS**: PGC-MDD 2025 (no 23andMe, no UKBB)
- **Brain Imaging**: UK Biobank imaging-derived phenotypes
- **Single-cell eQTL**: Brain cell type-specific datasets

### Example Data
Example datasets are provided in the `data/` directory:
```
data/
├── example_gwas_summary_stats.tsv     # Synthetic GWAS data
├── brain_cell_annotations.csv         # Cell type annotations
└── data_download_instructions.md      # Links to public datasets
```

## Directory Structure

```
MDD-scMR-MultiOmics/
├── R/                          # Modular R functions
├── scripts/                    # Main analysis workflows
├── data/                       # Example and reference data
├── results/                    # Analysis outputs (gitignored)
├── docs/                       # Additional documentation
└── tests/                      # Unit tests
```

## Analysis Pipeline

### 1. eQTL-MR Analysis (`scripts/01_eqtl_mr_analysis.R`)
- Loads eQTLGen expression data
- Performs LD clumping
- Conducts two-sample MR with multiple methods
- Generates forest plots and scatter plots

### 2. pQTL-MR Analysis (`scripts/02_pqtl_mr_analysis.R`) 
- Uses deCODE protein QTL data
- Cross-validates eQTL findings at protein level
- Produces comparative analysis results

### 3. Single-cell Validation (`scripts/03_single_cell_eqtl_mr.R`)
- Brain cell type-specific eQTL analysis
- Validates bulk tissue findings
- Identifies cell-type-specific effects

### 4. Mediation Analysis (`scripts/04_mediation_analysis.R`)
- Tests Gene → Brain Imaging → MDD pathways
- Calculates direct and indirect effects
- Bootstrap confidence intervals

### 5. Sensitivity Analysis (`scripts/05_reverse_mr_sensitivity.R`)
- Reverse MR analysis
- Pleiotropy testing (MR-Egger intercept)
- Heterogeneity assessment (Cochran's Q)
- MR-PRESSO outlier detection

### 6. Figure Generation (`scripts/06_generate_figures.R`)
- Volcano plots
- Forest plots  
- Venn diagrams
- Composite publication figures

## Key Functions

### Core Analysis Functions
```r
# Primary MR analysis
perform_mr_analysis()          # Two-sample MR with multiple methods
run_gsmr_analysis()           # GSMR with HEIDI outlier detection  
calculate_mediation_effect()  # Mediation analysis with bootstrap

# Data processing
load_gwas_data()              # Standardized GWAS data loading
harmonize_exposure_outcome()  # SNP harmonization
perform_ld_clumping()         # LD-based SNP clumping

# Visualization
plot_volcano()                # Volcano plots
plot_forest()                 # Forest plots  
plot_mr_scatter()             # MR scatter plots
create_composite_figure()     # Multi-panel figures
```

### Utility Functions
```r
# Quality control
check_data_format()           # Validate input data
calculate_f_statistics()      # Instrument strength assessment
steiger_directional_test()    # Causal direction testing

# Results processing
export_results()              # Standardized result export
generate_summary_table()     # Publication-ready tables
```

## Results Structure

Results are organized by analysis type:
```
results/
├── eqtl_mr/                 # eQTL-MR results
│   ├── significant/         # Significant associations
│   └── sensitivity/         # Sensitivity analyses
├── pqtl_mr/                 # pQTL-MR results  
├── single_cell/             # Single-cell validation
├── mediation/               # Mediation analyses
└── figures/                 # All publication figures
    ├── volcano_plots/
    ├── forest_plots/
    └── composite/
```

## Citation

If you use this code in your research, please cite:

```bibtex
@article{hou2025multiomics,
  title={Multi-Omics Mendelian Randomization Identifies CKAP2 as a Causal Protective Factor in Major Depressive Disorder},
  author={Hou, Xuejiao and Ying Huang}
}
```

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup
```r
# Clone the repository
git clone https://github.com/yourusername/MDD-scMR-MultiOmics.git

# Install development dependencies
renv::restore()

# Run tests
source("tests/test_mr_functions.R")
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Xuejiao Hou**
- Email: [houxj5@mail.sysu.edu.cn]
- ORCID: [0000-0001-9562-7540]
- Institution: [Department of Psychiatry, The Fifth Affiliated Hospital of Sun Yat-sen University, 52 Meihua East Road, Xiangzhou District, Zhuhai, Guangdong 519000, China]

## Acknowledgments

- eQTLGen Consortium for expression QTL data
- deCODE genetics for protein QTL data  
- PGC-MDD working group for depression GWAS data
- UK Biobank for brain imaging data
- R/Bioconductor community for statistical packages

## Version History

- **v1.0.0** (2025-02-10): Initial release
  - Complete multi-omics MR pipeline
  - Comprehensive documentation
  - Example datasets and workflows

## Troubleshooting

### Common Issues

**Memory Issues with Large Datasets:**
```r
# Increase memory limits
options(java.parameters = "-Xmx8g")
gc()  # Force garbage collection
```

**Missing Dependencies:**
```r
# Check and install missing packages
renv::status()
renv::restore()
```

**Path Issues:**
- Ensure all file paths use forward slashes (/) or `file.path()`
- Check working directory with `getwd()`
- Use relative paths where possible

For more help, see the [issues page](https://github.com/yourusername/MDD-scMR-MultiOmics/issues) or contact the authors.
