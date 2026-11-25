Leak-Free Transcriptomic Machine Learning & Deep Learning Benchmarking Pipeline (TCGA-KIRC)

An end-to-end, reproducible bioinformatics and statistical machine learning pipeline in R for bulk RNA-sequencing biomarker discovery, model calibration, and clinical decision-support benchmarking in Clear Cell Renal Cell Carcinoma (**TCGA-KIRC**).

## Overview

High-dimensional transcriptomic classification ($p \gg n$) frequently suffers from data leakage, severe multicollinearity, and complete separation in unregularized generalized linear models. This pipeline enforces strict train/test isolation prior to normalization-dependent feature selection and evaluates seven distinct predictive architectures across discrimination, probability calibration, and clinical net benefit:

1. **Lasso Logistic Regression** (`glmnet`, $L_1$ penalty)
2. **Elastic-Net Logistic Regression** (`glmnet`, $L_1 + L_2$ penalty)
3. **Random Forest** (`ranger`, impurity-based feature importance)
4. **Gradient Boosted Trees** (`xgboost`)
5. **Support Vector Machine with Radial Basis Function Kernel** (`kernlab`)
6. **Stacked Meta-Ensemble** (`stacks`, penalized meta-learner)
7. **Deep Residual Multi-Layer Perceptron** (native R `torch` with BatchNorm, SiLU, Dropout, and AdamW)

---

## Key Methodological Features

* **Automated Open-Access Data Ingestion:** Streams harmonized GRCh38 protein-coding read counts for TCGA-KIRC (`Primary Tumor` vs. `Solid Tissue Normal`) via `recount3` and applies `edgeR` Trimmed Mean of M-values (TMM) normalization.
* **Strict Leak-Free Partitioning:** Performs a stratified 75/25 train/holdout split before differential expression (`limma-voom`), variance filtering (top 250 variable genes), and z-score standardization.
* **Inferential & Variance-Corrected Testing:**
  * Multicollinearity diagnostics via Variance Inflation Factors (VIF) and unregularized binomial GLM breakdown demonstration.
  * 100-replicate stratified bootstrap confidence intervals and selection frequencies for Lasso odds ratios.
  * Repeated stratified cross-validation ($3 \times 5$ folds) compared via **Nadeau–Bengio variance-corrected paired $t$-tests**.
  * Independent holdout ROC-AUC comparisons via **paired DeLong tests** (`pROC`).
* **Calibration & Clinical Utility:** Evaluates Brier score, Expected Calibration Error (ECE), Maximum Calibration Error (MCE), Spiegelhalter's $Z$-test, and **Decision Curve Analysis (DCA)** net clinical benefit across decision thresholds ($p_t \in [0.10, 0.90]$).

---

## Generated Outputs

Running the pipeline automatically creates `tables/` and `figures/` subdirectories containing **16 CSV tables** and **16 publication-ready figures** (300 DPI `.png` and vector `.pdf`):

### Tables (`tables/`)
| File | Description |
| :--- | :--- |
| `Table_01_Cohort_Sample_Summary.csv` | Diagnostic class sample counts, proportions, and median library sizes |
| `Table_02_Sequencing_QC_Metrics.csv` | Filtering counts, mapped read depth, and TMM normalization summaries |
| `Table_03_Differential_Expression_Top30.csv` | Top 30 differentially expressed protein-coding genes (`limma-voom`, training fold) |
| `Table_04_Top25_Biomarker_Summary.csv` | Class-stratified expression mean/SD and inferential statistics for top variable genes |
| `Table_05_Multicollinearity_VIF.csv` | Variance Inflation Factors and tolerance across top 20 biomarkers |
| `Table_06_PCA_Variance_Explained.csv` | Eigenvalue standard deviations and cumulative variance across PC1–PC10 |
| `Table_07_Unregularized_GLM_Instability.csv` | Coefficient and odds-ratio divergence under unregularized high-dimensional GLM |
| `Table_08_Optimal_Hyperparameters.csv` | Selected tuning hyperparameters across all parametric and non-parametric models |
| `Table_09_Lasso_Odds_Ratios_Bootstrap_CI.csv` | Standardized Lasso log-odds, odds ratios, 95% bootstrap CIs, and selection frequencies |
| `Table_10_Cross_Validation_Performance.csv` | Out-of-fold $3 \times 5$ CV mean, SD, and 95% CIs (ROC-AUC, Accuracy, Brier, MCC) |
| `Table_11_Nadeau_Bengio_Corrected_Tests.csv` | Variance-corrected paired $t$-tests comparing ML architectures against Lasso baseline |
| `Table_12_Holdout_Test_Performance.csv` | Independent holdout classification metrics with exact binomial and DeLong CIs |
| `Table_13_DeLong_ROC_Comparison_Tests.csv` | Paired DeLong $Z$-statistics and $p$-values on holdout ROC curves |
| `Table_14_Probability_Calibration_Metrics.csv` | Brier score, ECE, MCE, and Spiegelhalter calibration $Z$-test results |
| `Table_15_Threshold_Sensitivity_Net_Benefit.csv` | Sensitivity, specificity, error counts, and DCA net benefit across thresholds |
| `Table_16_Consensus_Biomarker_Importance.csv` | Consensus feature ranking combining Random Forest Gini impurity and Lasso $\beta$ |

### Figures (`figures/`)
* `Figure_01` – Sequencing library size and TMM normalization distributions
* `Figure_02` – Differential expression volcano plot (`limma-voom`)
* `Figure_03` – Pairwise Pearson correlation heatmap of top 20 biomarkers
* `Figure_04` – Principal Component Analysis (PCA) biplot with 95% confidence ellipses
* `Figure_05` – Non-linear UMAP manifold embedding
* `Figure_06` – Lasso cross-validation regularization penalty path
* `Figure_07` – Forest plot of top 15 Lasso odds ratios with 95% bootstrap CIs
* `Figure_08` – `torch` Residual MLP training/validation BCE loss and ROC-AUC convergence
* `Figure_09` – Out-of-fold cross-validation metric distributions across models
* `Figure_10` – Independent holdout ROC curve overlay
* `Figure_11` – Precision-Recall curve overlay
* `Figure_12` – Holdout confusion matrix grid across all 7 architectures
* `Figure_13` – Probability calibration reliability diagrams
* `Figure_14` – Decision Curve Analysis (net clinical benefit vs. threshold probability)
* `Figure_15` – Top 20 Random Forest mean decrease in Gini impurity
* `Figure_16` – Sensitivity vs. specificity threshold trade-off curves

---

## Requirements & Installation

* **R** ($\ge 4.3.0$)
* **CRAN Packages:** `tidymodels`, `stacks`, `glmnet`, `ranger`, `xgboost`, `kernlab`, `torch`, `pROC`, `uwot`, `dplyr`, `tidyr`, `tibble`, `purrr`, `ggplot2`, `patchwork`, `scales`, `matrixStats`, `coro`
* **Bioconductor Packages:** `recount3`, `SummarizedExperiment`, `edgeR`, `limma`

Missing dependencies and `LibTorch` binaries are installed automatically on first run.

---

## Usage

1. Update `out_dir` at the top of the script to your desired local working directory.
2. Execute the script cleanly from R or RStudio:
   ```r
   source("path/to/genomic_ml_pipeline.R")
