# JM4QTN: Joint Mapping for Quantitative Trait Nucleotides

[![R-CMD-check](https://github.com/JunhuiLi1017/JM4QTN/workflows/R-CMD-check/badge.svg)](https://github.com/JunhuiLi1017/JM4QTN/actions)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Overview

**JM4QTN** is a comprehensive R package for joint mapping analysis of quantitative trait nucleotides (QTN) in genetic studies. The package implements sophisticated statistical methods for both Association Mapping (AM) and Linkage Mapping (LM) approaches to identify QTL affecting multiple traits simultaneously.

## Key Features

### 🔬 **Core Analysis Functions**
- **Joint Mapping Analysis**: Multivariate QTL analysis for multiple traits
- **Permutation Testing**: Empirical significance threshold determination
- **Linkage-Pleiotropy Analysis**: Distinguish between pleiotropic and linked QTL effects
- **Missing Genotype Imputation**: Using flanking marker information
- **Virtual Marker Creation**: For improved mapping resolution

### 📊 **Statistical Analysis**
- **Phenotype Analysis**: Normality tests, ANOVA, and least squares means
- **Stepwise Regression**: Multiple information criteria for cofactor selection
- **Parametric Bootstrap Testing**: For hypothesis testing
- **Population Support**: F2, RIL, DH, backcross, and advanced generations

### 🧬 **Population Types Supported**
- **F2**: F2 generation populations
- **RIL**: Recombinant Inbred Lines
- **DH**: Doubled Haploids
- **BCP1/BCP2**: Backcross to Parent 1/2
- **Fn**: Advanced generation populations (n > 2)

## Installation

```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("JunhuiLi1017/JM4QTN")

# Load the package
library(JM4QTN)
```

## Quick Start

### 1. Statistical Analysis of Phenotype Data

```r
# Analyze phenotype data with comprehensive statistics
results <- analyze_phenotype_statistics(pheno_data)

# View normality test results
results$Height$normality_test

# View ANOVA results
results$Height$ANOVA

# View least squares means
results$Height$lsmeans
```

### 2. Genotype Probability Calculations

```r
# Calculate genotype probabilities and impute missing data
result_am <- calculate_genotype_probabilities(genetic_map, geno_data, method = "AM")

# Linkage mapping with imputation for F2 population
result_lm <- calculate_genotype_probabilities(genetic_map, geno_data, method = "LM", 
                                            croType = "F2", steps = 0)

# Create virtual markers for improved resolution
result_vm <- calculate_genotype_probabilities(genetic_map, geno_data, method = "LM", 
                                            croType = "F2", steps = 5)
```

### 3. Permutation Testing

```r
# Permutation test for joint mapping
thresholds_am <- permutation_test_joint_mapping(c("Trait1", "Trait2"), pheno_data, geno_data, 
                                               method = "AM", npt = 100, alpha = 0.1)

# Permutation test for linkage-pleiotropy analysis
thresholds_lp <- permutation_test_linkage_pleiotropy(c("Trait1", "Trait2"), pheno_data, geno_data, 
                                                    npt = 100, alpha = 0.1)
```

### 4. Joint Mapping Analysis

```r
# Joint mapping analysis - Association mapping
results_am <- joint_mapping_analysis(c("Trait1", "Trait2"), h2, pheno_data, method = "AM",
                                    thresholds_am, geno_est)

# Joint mapping analysis - Linkage mapping
results_lm <- joint_mapping_analysis(c("Trait1", "Trait2"), h2, pheno_data, method = "LM",
                                    thresholds_lm, geno_est, geno_qtl)
```

### 5. Linkage-Pleiotropy Analysis

```r
# Linkage-pleiotropy analysis
results <- linkage_pleiotropy_analysis(c("Trait1", "Trait2"), pheno_data, geno_est, geno_qtl,
                                      thresholds, CChr = 1, Interval = c(0, 100), nPB = 100)
```

## Function Reference

### Core Analysis Functions

| Function | Description |
|----------|-------------|
| `joint_mapping_analysis()` | Main function for joint mapping analysis |
| `permutation_test_joint_mapping()` | Permutation testing for joint mapping |
| `permutation_test_linkage_pleiotropy()` | Permutation testing for linkage-pleiotropy |
| `linkage_pleiotropy_analysis()` | Linkage-pleiotropy analysis |

### Data Processing Functions

| Function | Description |
|----------|-------------|
| `calculate_genotype_probabilities()` | Genotype probability calculations and imputation |
| `calculate_genotype_frequencies()` | Genotype frequency calculations |
| `calculate_expected_genotype_distribution()` | Expected genotype distribution calculations |
| `haldane_mapping_function()` | Haldane mapping function for recombination |

### Statistical Analysis Functions

| Function | Description |
|----------|-------------|
| `analyze_phenotype_statistics()` | Comprehensive phenotype statistical analysis |

## Data Requirements

### Phenotype Data Structure
```r
pheno_data <- data.frame(
  Indi = paste0("Ind", 1:100),        # Individual identifiers
  Popu = rep(c("Pop1", "Pop2"), 50),  # Population identifiers
  MA = rep(1:2, 50),                  # Maternal allele information
  PA = rep(1:2, 50),                  # Paternal allele information
  Trait1 = rnorm(100, 100, 15),       # Trait measurements
  Trait2 = rnorm(100, 50, 8)          # Additional traits...
)
```

### Genetic Map Structure
```r
genetic_map <- data.frame(
  marker = c("M1", "M2", "M3"),       # Marker names
  chr = c(1, 1, 1),                   # Chromosome numbers
  pos = c(0, 10, 20)                  # Genetic positions (cM)
)
```

### Genotype Data Structure
```r
geno_data <- matrix(
  c(2, 0, 2, 1, 0, 1),               # Genotype codes: 0, 1, 2, 9/NA
  nrow = 3, ncol = 2,
  dimnames = list(c("Ind1", "Ind2", "Ind3"), c("M1", "M2"))
)
```

## Output Files

The package generates output files in organized directories:

- `OUTPUT_GenoProb/`: Genotype probability results
- `OUTPUT_ptJM/`: Permutation test results for joint mapping
- `OUTPUT_ptLP/`: Permutation test results for linkage-pleiotropy
- `OUTPUT_JM/`: Joint mapping analysis results
- `OUTPUT_LP/`: Linkage-pleiotropy analysis results

## Citation

If you use JM4QTN in your research, please cite:

```bibtex
@Manual{JM4QTN,
  title = {JM4QTN: Joint Mapping for Quantitative Trait Nucleotides},
  author = {Junhui Li and Wenxin Liu},
  year = {2025},
  note = {R package version 2.0.0},
  url = {https://github.com/JunhuiLi1017/JM4QTN}
}
```

## References

- Jiang, C. and Zeng, Z.B. (1995). Multiple trait analysis of genetic mapping for quantitative trait loci. *Genetics*, 140(3), 1111-1127.
- Churchill, G.A. and Doerge, R.W. (1994). Empirical threshold values for quantitative trait mapping. *Genetics*, 138(3), 963-971.
- Haldane, J.B.S. (1919). The combination of linkage values and the calculation of distances between the loci of linked factors. *Journal of Genetics*, 8(3), 299-309.

## Contributing

We welcome contributions! Please feel free to submit issues and pull requests.

## License

This project is licensed under the GNU General Public License v2.0 - see the [GNU GPLv2 license text](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html) for details.

## Contact

- **Maintainer**: Junhui Li <junhuili@cau.edu.cn>
- **GitHub**: [https://github.com/JunhuiLi1017/JM4QTN](https://github.com/JunhuiLi1017/JM4QTN)
