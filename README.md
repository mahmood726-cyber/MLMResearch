# MLMResearch: Multilevel Meta-Analysis Method Evaluation and Improvement

[![CRAN status](https://www.r-pkg.org/badges/version/MLMResearch)](https://CRAN.R-project.org/package=MLMResearch)

A research companion package to [MLM501](https://github.com/mahmood789/MLM501) that systematically evaluates modern multilevel meta-analysis (MLMA) methods across 501+ empirical Cochrane review datasets. 

While most MLMA methods (like REML and RVE) are evaluated via idealized simulations, `MLMResearch` provides a framework for **empirical stress-testing**. It identifies methodological flaws in real-world data, tests performance across different outcome types, and develops adaptive algorithms for small-sample corrections.

## Key Comparison

| Feature | `metafor` | `clubSandwich` | `robumeta` | `MLMResearch` |
| :--- | :---: | :---: | :---: | :---: |
| MLMA Model Fitting | ✓ | - | ✓ | ✓ (Wrappers) |
| Robust Variance Estimation | - | ✓ | ✓ | ✓ (Wrappers) |
| **Mass-Testing Framework** | - | - | - | **✓** |
| **Automated Issue Detection** | - | - | - | **✓** |
| **Adaptive RVE Selection** | - | - | - | **✓** |
| **Empirical Failure Pattern Analysis** | - | - | - | **✓** |

## Overview

This package provides tools to:

1. **Test** all modern multilevel meta-analysis methods systematically across large data cohorts.
2. **Detect** and catalog methodological issues (e.g., Variance Collapse, Convergence Failures).
3. **Analyze** patterns in method performance across different data topologies.
4. **Apply** improved, adaptive methods for small samples and sparse data.

## Installation

```r
# Install dependencies
install.packages(c("metafor", "clubSandwich"))

# Load the package
library(MLMResearch)
```

## Quick Start: The Evaluation Workflow

### 1. Load and Test on Real-World Data

```r
# Load dichotomous cohort (requires MLM501 data)
coh <- get_dichotomous_cohort()

# Run systematic tests (REML vs. RVE) on a sample
test_results <- run_dichotomous_tests(
  max_analyses = 50,
  methods = c("REML", "RVE-CR2", "RVE-CR1")
)
```

### 2. Automated Issue Detection

```r
# Detect all methodological failures (e.g., variance collapse)
all_issues <- batch_log_issues(test_results)
issue_report <- generate_issue_report()

print(issue_report$summary_table)
```

### 3. Apply Adaptive Methods

```r
# For datasets where standard methods fail due to small k, 
# use the adaptive RVE algorithm.
df_problematic <- coh[coh$analysis_id == "Analysis_005", ]
adaptive_result <- adaptive_rve(df_problematic)
```

## Scientific Contribution

`MLMResearch` shifts the evaluation of MLMA methods from theoretical simulation to massive, empirical validation. By exposing the high prevalence of issues like small-sample bias and variance component collapse in standard methods, the package provides a critical, automated safeguard for applied evidence synthesis.

## Dependencies

- **metafor**: Core meta-analysis functionality
- **clubSandwich**: Robust variance estimation

## Citation

If you use this package in research, please cite:

```
Ahmad M. MLMResearch: Multilevel Meta-Analysis Method Evaluation
and Improvement. R package version 0.1.0. 2024.
```

## License

MIT + file LICENSE
