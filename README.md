# JM4QTN: Joint Mapping for Quantitative Trait Loci

[![R-CMD-check](https://github.com/JunhuiLi1017/JM4QTN/workflows/R-CMD-check/badge.svg)](https://github.com/JunhuiLi1017/JM4QTN/actions)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Overview

**JM4QTN** supports joint mapping and related analyses for quantitative trait loci (QTL). The package provides phenotype summaries, genotype probability and imputation (association vs linkage mapping), permutation-based thresholds with stepwise regression (via **StepReg**), and scan statistics comparing full vs reduced linear models (**joint_map**).

## Key Features

### Core analysis

- **Permutation testing** (`permutation_test`) for empirical p-value and LOD cutoffs under stepwise selection  
- **Skeleton model** (`skeleton_build` / `skeletion_build`) fitting the selected model using those cutoffs  
- **Joint mapping scan** (`joint_map`) over candidate terms vs the skeleton model  

### Genotypes and map

- **Genotype probabilities / imputation** (`genotype_prob`) for AM or LM, optional virtual markers  
- **Expected genotype distribution** (`expected_genotype_dist`) and **genotype frequencies** (`genotype_freq`)  
- **Haldane map** (`haldane_map`) for recombination  

### Phenotypes

- **Phenotype statistics** (`pheno_stats`): normality, ANOVA, least squares means (via **lsmeans**)

## Installation

```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("JunhuiLi1017/JM4QTN")

## or install from CRAN
install.pakcages("JM4QTN")

# Load the package
library(JM4QTN)
```

## Quick start

### 1. Phenotype analysis (`pheno_stats`)

Designed for multi-environment trials: columns **E**, **B**, **R**, **G**, then traits.

```r
pheno_data <- data.frame(
  E = rep(c("Env1", "Env2"), each = 60),
  B = rep(c("B1", "B2"), each = 30, times = 2),
  R = rep(1:5, 24),
  G = factor(rep(1:12, 10)),
  Height = rnorm(120, 175, 8),
  Weight = rnorm(120, 75, 12)
)

results <- pheno_stats(pheno_data)
results$Height$normality_test
results$Height$ANOVA
results$Height$lsmeans
```

### 2. Genotype probabilities (`genotype_prob`)

```r
genetic_map <- data.frame(
  marker = c("M1", "M2", "M3"),
  chr = c(1, 1, 1),
  pos = c(0, 10, 20)
)
geno_data <- matrix(sample(0:2, 15, TRUE), nrow = 5, ncol = 3,
                    dimnames = list(paste0("Ind", 1:5), c("M1", "M2", "M3")))

result_am <- genotype_prob(genetic_map, geno_data, method = "AM")
result_lm <- genotype_prob(genetic_map, geno_data, method = "LM",
                           croType = "F2", steps = 0)
```

### 3. Permutation thresholds, skeleton model, and joint mapping

Formula-based workflow with genotype columns bound into `data` (see package examples for `joint_map`).

```r
pheno_data <- data.frame(
  Trait1 = rnorm(100, 100, 15),
  Trait2 = rnorm(100, 50, 8),
  Popu = rep(c("Pop1", "Pop2"), each = 50)
)
geno_data <- matrix(sample(0:2, 100 * 50, TRUE), nrow = 100, ncol = 50)
colnames(geno_data) <- paste0("M", 1:50)
data1 <- cbind(pheno_data, geno_data)

terms <- c("Popu", paste0(colnames(geno_data), ":Popu"))
formula1 <- reformulate(terms, response = "Trait1")

cut_off_list <- permutation_test(formula1, data1, n = 100, alpha = 0.1)

skeleton <- skeleton_build(
  formula1, data1,
  strategy = "bidirection", metric = "SL",
  cut_off_list = cut_off_list
)

results <- joint_map(
  formula1, data1, skeleton,
  include = "Popu", cut_off_list = cut_off_list
)
results$p_value
results$lod
```

### 4. Supporting utilities

```r
haldane_map(0.1)   # genetic distance (Morgan) -> recombination fraction
genotype_freq("Fn", generation = 3, genotype_index = 1, recomb_aq = 0.1, recomb_qb = 0.2)
expected_genotype_dist("22", "F2", Gn = 2, x = 0.1, y = 0.2)
```

## Function reference (exported)

| Function | Description |
|----------|-------------|
| `joint_map()` | Compare candidate model terms to a fitted skeleton; returns p-values and LOD-like statistics |
| `permutation_test()` | Permutation distribution and empirical cutoffs for stepwise linear models |
| `skeleton_build()` | Alias for `skeletion_build()` (preferred spelling) |
| `skeletion_build()` | Fit stepwise skeleton using permutation p-value as StepReg entry/stay levels |
| `genotype_prob()` | Genotype probabilities / imputation (AM or LM) |
| `genotype_freq()` | Genotype class frequencies by cross type and generation |
| `expected_genotype_dist()` | Expected genotype distribution by marker pattern and cross |
| `haldane_map()` | Haldane mapping function |
| `pheno_stats()` | Normality tests, ANOVA, lsmeans for structured phenotype tables |

**Backward-compatible names** (same package, alternate spelling): `calculate_genotype_frequencies`, `calculate_expected_genotype_distribution`, `haldane_mapping_function` — see individual help pages.

## Data notes

- **`pheno_stats`**: expects **E**, **B**, **R**, **G** as the first four columns, then traits.  
- **`genotype_prob`**: rows = individuals, columns = markers; column order must align with `GeneticMap` rows.  
- **`joint_map` / `permutation_test`**: supply a single `data` frame containing phenotypes, population, and numeric genotype columns used in the formula.

## Citation

```bibtex
@Manual{JM4QTN,
  title = {JM4QTN: Joint Mapping for Quantitative Trait Loci},
  author = {Junhui Li and Wenxin Liu},
  year = {2026},
  note = {R package version 1.0.0},
  url = {https://github.com/JunhuiLi1017/JM4QTN}
}
```

## References

- Jiang, C. and Zeng, Z.B. (1995). Multiple trait analysis of genetic mapping for quantitative trait loci. *Genetics*, 140(3), 1111-1127.
- Churchill, G.A. and Doerge, R.W. (1994). Empirical threshold values for quantitative trait mapping. *Genetics*, 138(3), 963-971.
- Haldane, J.B.S. (1919). The combination of linkage values and the calculation of distances between the loci of linked factors. *Journal of Genetics*, 8(3), 299-309.

## Contributing

Issues and pull requests are welcome on [GitHub](https://github.com/JunhuiLi1017/JM4QTN).

## License

This project is licensed under the GNU General Public License v2.0 - see the [GNU GPLv2 license text](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html) for details.

## Contact

- **Maintainer**: Junhui Li <junhuili@cau.edu.cn>
- **GitHub**: [https://github.com/JunhuiLi1017/JM4QTN](https://github.com/JunhuiLi1017/JM4QTN)
