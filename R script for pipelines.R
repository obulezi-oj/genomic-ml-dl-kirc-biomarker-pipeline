local({
  # ============================================================================
  # PUBLICATION-GRADE BIOINFORMATICS & GENOMIC ML/DL BENCHMARKING PIPELINE
  # Output Directory: C:\Users\Dr. O. J. Obulezi\Documents\R projects\genomic
  # Generates: 16 Publication Tables (.csv) + 16 Clean Figures (.png & .pdf)
  # ============================================================================
  
  # 0. PATCH LIBRARY & LOCK CLEANUP
  out_dir   <- "C:/Users/Dr. O. J. Obulezi/Documents/R projects/genomic"
  patch_lib <- file.path(out_dir, "r_lib_patch")
  dir.create(patch_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(patch_lib, .libPaths()))
  
  for (lp in .libPaths()) {
    locks <- list.files(lp, pattern = "^00LOCK", full.names = TRUE)
    if (length(locks) > 0) unlink(locks, recursive = TRUE, force = TRUE)
  }
  
  if (packageVersion("rlang") < "1.3.0") {
    install.packages("rlang", lib = patch_lib, repos = "https://cloud.r-project.org", type = "binary")
    for (pkg in c("purrr", "tibble", "dplyr", "tidyr", "ggplot2", "pillar", "vctrs", "lifecycle", "cli", "glue", "rlang")) {
      if (pkg %in% loadedNamespaces()) try(unloadNamespace(pkg), silent = TRUE)
    }
    loadNamespace("rlang", lib.loc = patch_lib)
  }
  
  # 1. PACKAGE SETUP
  cran_pkgs <- c("tidymodels", "stacks", "glmnet", "ranger", "xgboost", "kernlab", "torch", "pROC", "uwot", "dplyr", "tidyr", "tibble", "purrr", "ggplot2", "patchwork", "scales", "matrixStats", "coro")
  new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
  if (length(new_cran)) install.packages(new_cran, lib = patch_lib, repos = "https://cloud.r-project.org", type = "binary")
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = patch_lib, repos = "https://cloud.r-project.org", type = "binary")
  bioc_pkgs <- c("recount3", "SummarizedExperiment", "edgeR", "limma")
  new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
  if (length(new_bioc)) BiocManager::install(new_bioc, lib = patch_lib, update = FALSE, ask = FALSE)
  
  suppressPackageStartupMessages({
    library(recount3); library(SummarizedExperiment); library(edgeR); library(limma)
    library(tidymodels); library(stacks); library(torch); library(pROC); library(uwot)
    library(dplyr); library(tidyr); library(tibble); library(purrr)
    library(ggplot2); library(patchwork); library(scales); library(matrixStats)
  })
  
  if (!torch_is_installed()) install_torch()
  set.seed(42)
  torch_manual_seed(42)
  
  tab_dir <- file.path(out_dir, "tables")
  fig_dir <- file.path(out_dir, "figures")
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  save_table <- function(df, filename) {
    write.csv(df, file.path(tab_dir, filename), row.names = FALSE)
    message(sprintf("Saved Table: %s", filename))
  }
  
  save_plot <- function(plot_obj, filename_base, width = 9, height = 6) {
    clean_plot <- plot_obj + labs(title = "", subtitle = "") + theme(plot.title = element_blank(), plot.subtitle = element_blank())
    ggsave(file.path(fig_dir, paste0(filename_base, ".png")), plot = clean_plot, width = width, height = height, dpi = 300, bg = "white")
    ggsave(file.path(fig_dir, paste0(filename_base, ".pdf")), plot = clean_plot, width = width, height = height, device = cairo_pdf)
    message(sprintf("Saved Figure: %s", filename_base))
  }
  
  theme_pub <- theme_minimal(base_size = 12) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(face = "bold"), legend.position = "bottom", strip.text = element_text(face = "bold"))
  
  # 2. DATA INGESTION & GENOMIC QC
  message("Streaming open-access TCGA-KIRC cohort via recount3...")
  human_projects <- available_projects()
  kirc_info <- subset(human_projects, project == "KIRC" & project_home == "data_sources/tcga")
  rse <- create_rse(kirc_info)
  
  assays(rse)$counts <- transform_counts(rse)
  rse <- rse[rowData(rse)$gene_type == "protein_coding", ]
  
  sample_type <- colData(rse)$tcga.cgc_sample_sample_type
  keep_idx <- sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
  rse <- rse[, keep_idx]
  
  outcome_vec <- factor(ifelse(colData(rse)$tcga.cgc_sample_sample_type == "Primary Tumor", "Tumor", "Normal"), levels = c("Tumor", "Normal"))
  
  dge <- DGEList(counts = assays(rse)$counts)
  keep_genes <- filterByExpr(dge, group = outcome_vec)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
  
  # FIX 1: Synchronize clean gene symbols across BOTH log_cpm and dge$counts
  gene_symbols <- make.unique(make.names(rowData(rse)[rownames(log_cpm), "gene_name"]))
  rownames(log_cpm)    <- gene_symbols
  rownames(dge$counts) <- gene_symbols
  
  # TABLE 1: Cohort Demographics & Sample Summary
  t1_cohort <- tibble(
    Diagnostic_Class = levels(outcome_vec),
    Sample_Count = as.integer(table(outcome_vec)),
    Proportion_Pct = round(as.numeric(prop.table(table(outcome_vec))) * 100, 2),
    Median_Library_Size = as.numeric(tapply(dge$samples$lib.size, outcome_vec, median)),
    IQR_Library_Size = as.numeric(tapply(dge$samples$lib.size, outcome_vec, IQR)),
    Median_Norm_Factor = round(as.numeric(tapply(dge$samples$norm.factors, outcome_vec, median)), 4)
  )
  save_table(t1_cohort, "Table_01_Cohort_Sample_Summary.csv")
  
  # TABLE 2: Sequencing Library & QC Summary Statistics
  t2_qc <- data.frame(
    Metric = c("Total Samples", "Initial Protein-Coding Genes", "Retained Expressed Genes", "Library Size (Millions)", "TMM Normalization Factor"),
    Value_or_Mean = c(ncol(log_cpm), sum(rowData(rse)$gene_type == "protein_coding"), nrow(log_cpm), round(mean(dge$samples$lib.size) / 1e6, 2), round(mean(dge$samples$norm.factors), 4)),
    SD_or_Notes = c("Stratified 75/25 split", "GRCh38 annotation", "filterByExpr threshold", as.character(round(sd(dge$samples$lib.size) / 1e6, 2)), as.character(round(sd(dge$samples$norm.factors), 4))),
    stringsAsFactors = FALSE
  )
  save_table(t2_qc, "Table_02_Sequencing_QC_Metrics.csv")
  
  # FIGURE 1: Library Size & Normalization Distribution
  qc_df <- tibble(Class = outcome_vec, LibSize_M = dge$samples$lib.size / 1e6, NormFactor = dge$samples$norm.factors)
  p1_a <- ggplot(qc_df, aes(x = Class, y = LibSize_M, fill = Class)) + geom_violin(alpha = 0.7, trim = FALSE) + geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) + scale_fill_manual(values = c("Tumor" = "#D55E00", "Normal" = "#0072B2")) + labs(title = "", subtitle = "", main = "", x = NULL, y = "Mapped Reads (Millions)") + theme_pub + theme(legend.position = "none")
  p1_b <- ggplot(qc_df, aes(x = LibSize_M, y = NormFactor, color = Class)) + geom_point(alpha = 0.65, size = 2) + scale_color_manual(values = c("Tumor" = "#D55E00", "Normal" = "#0072B2")) + labs(title = "", subtitle = "", main = "", x = "Mapped Reads (Millions)", y = "TMM Factor") + theme_pub
  save_plot(p1_a + p1_b, "Figure_01_Sequencing_QC_Distributions", width = 10, height = 5)
  
  # 3. LEAK-FREE DATA PARTITIONING & INFERENTIAL BIOMARKER DISCOVERY
  X_all <- as.data.frame(t(log_cpm))
  genomic_df <- bind_cols(outcome = outcome_vec, X_all)
  
  data_split <- initial_split(genomic_df, prop = 0.75, strata = outcome)
  train_raw  <- training(data_split)
  test_raw   <- testing(data_split)
  
  train_counts <- dge$counts[, rownames(train_raw)]
  design_mat   <- model.matrix(~ 0 + train_raw$outcome)
  colnames(design_mat) <- c("Tumor", "Normal")
  v_obj <- voom(train_counts, design_mat, plot = FALSE)
  fit_lm <- lmFit(v_obj, design_mat)
  contr_mat <- makeContrasts(TumorVsNormal = Tumor - Normal, levels = design_mat)
  fit_eb <- eBayes(contrasts.fit(fit_lm, contr_mat))
  de_res <- topTable(fit_eb, number = Inf, sort.by = "P") %>% rownames_to_column("Gene")
  
  # TABLE 3: Top 30 Differentially Expressed Genes
  save_table(head(de_res, 30), "Table_03_Differential_Expression_Top30.csv")
  
  # FIGURE 2: Volcano Plot
  de_res <- de_res %>% mutate(Regulation = case_when(adj.P.Val < 0.01 & logFC > 2 ~ "Up-regulated in Tumor", adj.P.Val < 0.01 & logFC < -2 ~ "Down-regulated in Tumor", TRUE ~ "Not Significant"))
  top_labels <- head(de_res, 12)
  p2_volcano <- ggplot(de_res, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) + geom_point(alpha = 0.5, size = 1.3) + geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey40") + geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") + geom_text(data = top_labels, aes(label = Gene), size = 3, vjust = -0.6, check_overlap = TRUE, show.legend = FALSE) + scale_color_manual(values = c("Up-regulated in Tumor" = "#D55E00", "Down-regulated in Tumor" = "#0072B2", "Not Significant" = "grey70")) + labs(title = "", subtitle = "", main = "", x = "log2 Fold Change (Tumor vs. Normal)", y = "-log10 Adjusted P-Value") + theme_pub
  save_plot(p2_volcano, "Figure_02_Volcano_Differential_Expression", width = 9, height = 6)
  
  train_mat <- as.matrix(train_raw[, -1])
  gene_vars <- colVars(train_mat)
  top_250_genes <- names(sort(gene_vars, decreasing = TRUE))[1:250]
  
  # TABLE 4: Top 25 Variable Biomarkers (now matches clean gene symbols with zero NAs)
  t4_biomarkers <- tibble(
    Gene = top_250_genes[1:25],
    Mean_Tumor = colMeans(train_raw[train_raw$outcome == "Tumor", top_250_genes[1:25]]),
    SD_Tumor = colSds(as.matrix(train_raw[train_raw$outcome == "Tumor", top_250_genes[1:25]])),
    Mean_Normal = colMeans(train_raw[train_raw$outcome == "Normal", top_250_genes[1:25]]),
    SD_Normal = colSds(as.matrix(train_raw[train_raw$outcome == "Normal", top_250_genes[1:25]]))
  ) %>% left_join(de_res %>% select(Gene, logFC, t, adj.P.Val), by = "Gene")
  save_table(t4_biomarkers, "Table_04_Top25_Biomarker_Summary.csv")
  
  # 4. MULTICOLLINEARITY, PCA, & UMAP
  top20_mat <- scale(train_raw[, top_250_genes[1:20]])
  cor_mat <- cor(top20_mat)
  vif_vals <- diag(solve(cor_mat + diag(1e-5, ncol(cor_mat))))
  
  # TABLE 5: VIF
  t5_vif <- tibble(Gene = colnames(top20_mat), VIF = round(vif_vals, 3), Tolerance = round(1 / vif_vals, 4)) %>% arrange(desc(VIF))
  save_table(t5_vif, "Table_05_Multicollinearity_VIF.csv")
  
  # FIGURE 3: Correlation Heatmap
  cor_long <- as.data.frame(as.table(cor_mat))
  colnames(cor_long) <- c("Gene1", "Gene2", "Correlation")
  p3_cor <- ggplot(cor_long, aes(x = Gene1, y = Gene2, fill = Correlation)) + geom_tile() + scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-1, 1)) + labs(title = "", subtitle = "", main = "", x = NULL, y = NULL) + theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))
  save_plot(p3_cor, "Figure_03_Correlation_Heatmap", width = 8, height = 7)
  
  pca_fit <- prcomp(scale(train_raw[, top_250_genes]), center = FALSE, scale. = FALSE)
  var_exp <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  
  # TABLE 6: PCA Variance Explained
  t6_pca <- tibble(Component = paste0("PC", 1:10), Standard_Deviation = round(pca_fit$sdev[1:10], 4), Proportion_Variance = round(var_exp[1:10], 4), Cumulative_Proportion = round(cumsum(var_exp[1:10]), 4))
  save_table(t6_pca, "Table_06_PCA_Variance_Explained.csv")
  
  # FIGURE 4: PCA Biplot
  pca_df <- tibble(PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2], Class = train_raw$outcome)
  p4_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) + geom_point(alpha = 0.75, size = 2.2) + stat_ellipse(level = 0.95, linetype = "dashed") + scale_color_manual(values = c("Tumor" = "#D55E00", "Normal" = "#0072B2")) + labs(title = "", subtitle = "", main = "", x = sprintf("PC1 (%.1f%% Variance)", var_exp[1] * 100), y = sprintf("PC2 (%.1f%% Variance)", var_exp[2] * 100)) + theme_pub
  save_plot(p4_pca, "Figure_04_PCA_Biplot", width = 8, height = 6)
  
  # FIGURE 5: UMAP
  umap_coords <- umap(scale(train_raw[, top_250_genes]), n_neighbors = 15, min_dist = 0.1)
  umap_df <- tibble(UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2], Class = train_raw$outcome)
  p5_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) + geom_point(alpha = 0.75, size = 2.2) + scale_color_manual(values = c("Tumor" = "#D55E00", "Normal" = "#0072B2")) + labs(title = "", subtitle = "", main = "", x = "UMAP Dimension 1", y = "UMAP Dimension 2") + theme_pub
  save_plot(p5_umap, "Figure_05_UMAP_Manifold", width = 8, height = 6)
  
  # 5. UNREGULARIZED GLM INSTABILITY
  train_data <- train_raw[, c("outcome", top_250_genes)]
  test_data  <- test_raw[, c("outcome", top_250_genes)]
  
  glm_sub_df <- train_data[, c("outcome", top_250_genes[1:40])]
  glm_sub_df$y <- ifelse(glm_sub_df$outcome == "Tumor", 1, 0)
  unreg_glm <- suppressWarnings(glm(y ~ . - outcome, data = glm_sub_df, family = binomial()))
  glm_sum <- as.data.frame(summary(unreg_glm)$coefficients) %>% rownames_to_column("Term") %>% mutate(Odds_Ratio = exp(pmin(Estimate, 50)))
  
  # TABLE 7: Unregularized GLM Instability
  save_table(head(glm_sum, 20), "Table_07_Unregularized_GLM_Instability.csv")
  
  # 6. REPEATED STRATIFIED CROSS-VALIDATION & ENSEMBLE
  cv_folds <- vfold_cv(train_data, v = 5, repeats = 3, strata = outcome)
  genomic_recipe <- recipe(outcome ~ ., data = train_data) %>% step_zv(all_predictors()) %>% step_normalize(all_predictors())
  
  lasso_spec <- logistic_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet") %>% set_mode("classification")
  enet_spec  <- logistic_reg(penalty = tune(), mixture = tune()) %>% set_engine("glmnet") %>% set_mode("classification")
  rf_spec    <- rand_forest(mtry = tune(), min_n = tune(), trees = 500) %>% set_engine("ranger", importance = "impurity") %>% set_mode("classification")
  xgb_spec   <- boost_tree(trees = 300, tree_depth = tune(), learn_rate = tune(), mtry = tune()) %>% set_engine("xgboost") %>% set_mode("classification")
  svm_spec   <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>% set_engine("kernlab") %>% set_mode("classification")
  
  stack_ctrl  <- control_stack_grid()
  metrics_set <- metric_set(roc_auc, brier_class, accuracy, mcc)
  
  message("Running 3x5 Repeated Stratified Cross-Validation across ML architectures...")
  lasso_res <- tune_grid(workflow(genomic_recipe, lasso_spec), resamples = cv_folds, grid = 10, metrics = metrics_set, control = stack_ctrl)
  enet_res  <- tune_grid(workflow(genomic_recipe, enet_spec),  resamples = cv_folds, grid = 12, metrics = metrics_set, control = stack_ctrl)
  rf_res    <- tune_grid(workflow(genomic_recipe, rf_spec),    resamples = cv_folds, grid = 10, metrics = metrics_set, control = stack_ctrl)
  xgb_res   <- tune_grid(workflow(genomic_recipe, xgb_spec),   resamples = cv_folds, grid = 12, metrics = metrics_set, control = stack_ctrl)
  svm_res   <- tune_grid(workflow(genomic_recipe, svm_spec),   resamples = cv_folds, grid = 10, metrics = metrics_set, control = stack_ctrl)
  
  stacked_ensemble <- stacks() %>%
    add_candidates(lasso_res) %>% add_candidates(enet_res) %>% add_candidates(rf_res) %>% add_candidates(xgb_res) %>% add_candidates(svm_res) %>%
    blend_predictions(penalty = 10^(-4:-1), metric = metric_set(roc_auc)) %>% fit_members()
  
  fit_best_wf <- function(wf_obj, tune_res) {
    best_p <- select_best(tune_res, metric = "roc_auc")
    finalize_workflow(wf_obj, best_p) %>% fit(data = train_data)
  }
  
  final_lasso <- fit_best_wf(workflow(genomic_recipe, lasso_spec), lasso_res)
  final_enet  <- fit_best_wf(workflow(genomic_recipe, enet_spec),  enet_res)
  final_rf    <- fit_best_wf(workflow(genomic_recipe, rf_spec),    rf_res)
  final_xgb   <- fit_best_wf(workflow(genomic_recipe, xgb_spec),   xgb_res)
  final_svm   <- fit_best_wf(workflow(genomic_recipe, svm_spec),   svm_res)
  
  # TABLE 8: Optimal Hyperparameters
  t8_hyper <- bind_rows(
    select_best(lasso_res, metric = "roc_auc") %>% mutate(Model = "Lasso LR"),
    select_best(enet_res,  metric = "roc_auc") %>% mutate(Model = "Elastic-Net LR"),
    select_best(rf_res,    metric = "roc_auc") %>% mutate(Model = "Random Forest"),
    select_best(xgb_res,   metric = "roc_auc") %>% mutate(Model = "XGBoost"),
    select_best(svm_res,   metric = "roc_auc") %>% mutate(Model = "SVM-RBF")
  ) %>% select(Model, everything(), -.config)
  save_table(t8_hyper, "Table_08_Optimal_Hyperparameters.csv")
  
  # TABLE 9: Lasso Odds Ratios & Bootstrap CIs
  lasso_coefs <- tidy(extract_fit_parsnip(final_lasso)) %>% filter(term != "(Intercept)", estimate != 0) %>% mutate(OR_per_SD = exp(estimate)) %>% arrange(desc(abs(estimate)))
  boot_res <- bootstraps(train_data, times = 100, strata = outcome)
  best_pen <- select_best(lasso_res, metric = "roc_auc")$penalty
  lasso_fixed_wf <- finalize_workflow(workflow(genomic_recipe, lasso_spec), list(penalty = best_pen))
  
  boot_coefs <- map_dfr(boot_res$splits, function(spl) {
    fit_b <- fit(lasso_fixed_wf, data = analysis(spl))
    tidy(extract_fit_parsnip(fit_b)) %>% filter(term %in% lasso_coefs$term)
  }) %>% group_by(term) %>% summarise(OR_CI_Low = exp(quantile(estimate, 0.025)), OR_CI_High = exp(quantile(estimate, 0.975)), Selection_Freq = mean(estimate != 0), .groups = "drop")
  
  t9_lasso_or <- lasso_coefs %>% left_join(boot_coefs, by = "term") %>% select(Gene = term, Beta_per_SD = estimate, OR_per_SD, OR_CI_Low, OR_CI_High, Selection_Freq)
  save_table(t9_lasso_or, "Table_09_Lasso_Odds_Ratios_Bootstrap_CI.csv")
  
  # FIGURE 6: Lasso Regularization Path
  p6_lasso_path <- autoplot(lasso_res, metric = "roc_auc") + labs(title = "", subtitle = "", main = "") + theme_pub
  save_plot(p6_lasso_path, "Figure_06_Lasso_Regularization_Path", width = 8, height = 5)
  
  # FIGURE 7: Forest Plot of Top 15 Lasso Odds Ratios
  top15_or <- head(t9_lasso_or, 15)
  p7_forest <- ggplot(top15_or, aes(x = OR_per_SD, y = reorder(Gene, OR_per_SD))) + geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") + geom_errorbarh(aes(xmin = OR_CI_Low, xmax = OR_CI_High), height = 0.25, color = "grey30") + geom_point(aes(color = OR_per_SD > 1), size = 3) + scale_x_log10() + scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), labels = c("Protective (Lower Tumor Odds)", "Risk (Higher Tumor Odds)")) + labs(title = "", subtitle = "", main = "", x = "Odds Ratio (Log Scale)", y = NULL, color = "Direction") + theme_pub
  save_plot(p7_forest, "Figure_07_Lasso_Forest_Plot", width = 9, height = 6)
  
  # 7. DEEP LEARNING (`torch` RESIDUAL MLP)
  message("Training Deep Learning Neural Network (`torch` Residual MLP)...")
  X_train_raw <- as.matrix(train_data[, top_250_genes])
  X_test_raw  <- as.matrix(test_data[, top_250_genes])
  tr_center   <- colMeans(X_train_raw)
  tr_scale    <- colSds(X_train_raw)
  tr_scale[tr_scale == 0] <- 1
  
  X_train_mat <- sweep(sweep(X_train_raw, 2, tr_center, "-"), 2, tr_scale, "/")
  X_test_mat  <- sweep(sweep(X_test_raw,  2, tr_center, "-"), 2, tr_scale, "/")
  y_train_vec <- ifelse(train_data$outcome == "Tumor", 1, 0)
  y_test_vec  <- ifelse(test_data$outcome == "Tumor", 1, 0)
  
  val_idx  <- sample(seq_len(nrow(X_train_mat)), size = floor(0.15 * nrow(X_train_mat)))
  X_tr_sub <- X_train_mat[-val_idx, , drop = FALSE]
  y_tr_sub <- y_train_vec[-val_idx]
  X_va_sub <- X_train_mat[val_idx, , drop = FALSE]
  y_va_sub <- y_train_vec[val_idx]
  
  genomic_dataset <- dataset(
    name = "GenomicDataset",
    initialize = function(x, y) { self$x <- torch_tensor(x, dtype = torch_float()); self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2) },
    .getitem = function(i) list(x = self$x[i, ], y = self$y[i, ]),
    .length  = function() self$x$size()[[1]]
  )
  train_dl <- dataloader(genomic_dataset(X_tr_sub, y_tr_sub), batch_size = 32, shuffle = TRUE)
  
  GenomicResMLP <- nn_module(
    "GenomicResMLP",
    initialize = function(in_features) {
      self$input_proj <- nn_sequential(nn_linear(in_features, 128), nn_batch_norm1d(128), nn_silu(), nn_dropout(0.4))
      self$res_block  <- nn_sequential(nn_linear(128, 128), nn_batch_norm1d(128), nn_silu(), nn_dropout(0.3), nn_linear(128, 128), nn_batch_norm1d(128))
      self$act <- nn_silu()
      self$head <- nn_linear(128, 1)
    },
    forward = function(x) { h <- self$input_proj(x); h <- self$act(h + self$res_block(h)); self$head(h) }
  )
  
  device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
  dl_model  <- GenomicResMLP(ncol(X_train_mat))$to(device = device)
  optimizer <- optim_adamw(dl_model$parameters, lr = 1e-3, weight_decay = 1e-2)
  loss_fn   <- nn_bce_with_logits_loss()
  
  epochs <- 40
  dl_history <- tibble(Epoch = integer(), Train_Loss = numeric(), Val_Loss = numeric(), Val_AUC = numeric())
  
  for (ep in seq_len(epochs)) {
    dl_model$train()
    batch_losses <- c()
    coro::loop(for (b in train_dl) {
      optimizer$zero_grad()
      out <- dl_model(b$x$to(device = device))
      loss <- loss_fn(out, b$y$to(device = device))
      loss$backward()
      optimizer$step()
      batch_losses <- c(batch_losses, as.numeric(loss$item()))
    })
    dl_model$eval()
    with_no_grad({
      vx <- torch_tensor(X_va_sub, dtype = torch_float())$to(device = device)
      vy <- torch_tensor(y_va_sub, dtype = torch_float())$unsqueeze(2)$to(device = device)
      v_logits <- dl_model(vx)
      v_loss   <- as.numeric(loss_fn(v_logits, vy)$item())
      v_probs  <- as.numeric(torch_sigmoid(v_logits)$cpu())
      v_auc    <- as.numeric(auc(roc(y_va_sub, v_probs, quiet = TRUE)))
    })
    dl_history <- add_row(dl_history, Epoch = ep, Train_Loss = mean(batch_losses), Val_Loss = v_loss, Val_AUC = v_auc)
  }
  
  dl_model$eval()
  with_no_grad({
    tx <- torch_tensor(X_test_mat, dtype = torch_float())$to(device = device)
    dl_probs <- as.numeric(torch_sigmoid(dl_model(tx))$cpu())
  })
  
  # FIGURE 8: Deep Learning Convergence
  p8_dl_a <- ggplot(dl_history, aes(x = Epoch)) + geom_line(aes(y = Train_Loss, color = "Training BCE Loss"), linewidth = 1) + geom_line(aes(y = Val_Loss, color = "Validation BCE Loss"), linewidth = 1, linetype = "dashed") + scale_color_manual(values = c("Training BCE Loss" = "#0072B2", "Validation BCE Loss" = "#D55E00")) + labs(title = "", subtitle = "", main = "", y = "Binary Cross-Entropy Loss", color = NULL) + theme_pub
  p8_dl_b <- ggplot(dl_history, aes(x = Epoch, y = Val_AUC)) + geom_line(color = "#009E73", linewidth = 1.1) + labs(title = "", subtitle = "", main = "", y = "Validation ROC-AUC") + theme_pub
  save_plot(p8_dl_a + p8_dl_b, "Figure_08_Deep_Learning_Convergence", width = 10, height = 5)
  
  # 8. CROSS-VALIDATION SUMMARY & NADEAU-BENGIO TESTS
  extract_best_folds <- function(tune_res, model_name) {
    best_id <- select_best(tune_res, metric = "roc_auc")$.config
    collect_metrics(tune_res, summarize = FALSE) %>% filter(.config == best_id) %>% mutate(Model = model_name)
  }
  
  cv_all_folds <- bind_rows(
    extract_best_folds(lasso_res, "Lasso LR"),
    extract_best_folds(enet_res,  "Elastic-Net LR"),
    extract_best_folds(rf_res,    "Random Forest"),
    extract_best_folds(xgb_res,   "XGBoost"),
    extract_best_folds(svm_res,   "SVM-RBF")
  )
  
  # TABLE 10: Out-of-Fold CV Performance
  t10_cv_summary <- cv_all_folds %>% group_by(Model, .metric) %>% summarise(Mean = round(mean(.estimate), 4), SD = round(sd(.estimate), 4), CI_Low = round(pmax(0, Mean - 1.96 * SD / sqrt(n())), 4), CI_High = round(pmin(1, Mean + 1.96 * SD / sqrt(n())), 4), .groups = "drop") %>% pivot_wider(names_from = .metric, values_from = c(Mean, SD, CI_Low, CI_High))
  save_table(t10_cv_summary, "Table_10_Cross_Validation_Performance.csv")
  
  # FIGURE 9: CV Boxplots
  p9_cv_box <- cv_all_folds %>% filter(.metric %in% c("roc_auc", "accuracy", "brier_class")) %>% ggplot(aes(x = reorder(Model, .estimate, FUN = median), y = .estimate, fill = Model)) + geom_boxplot(alpha = 0.8, outlier.size = 1.5) + coord_flip() + facet_wrap(~ .metric, scales = "free_x") + labs(title = "", subtitle = "", main = "", x = NULL, y = "Metric Score") + theme_pub + theme(legend.position = "none")
  save_plot(p9_cv_box, "Figure_09_CV_Performance_Boxplots", width = 11, height = 5)
  
  # TABLE 11: Nadeau-Bengio Corrected Paired t-Tests
  lasso_auc_folds <- cv_all_folds %>% filter(Model == "Lasso LR", .metric == "roc_auc") %>% arrange(id, id2) %>% pull(.estimate)
  nadeau_bengio_test <- function(comp_model) {
    comp_auc <- cv_all_folds %>% filter(Model == comp_model, .metric == "roc_auc") %>% arrange(id, id2) %>% pull(.estimate)
    diffs <- lasso_auc_folds - comp_auc
    K <- length(diffs)
    var_corr <- (1 / K + 0.25) * var(diffs)
    t_stat <- ifelse(var_corr == 0, 0, mean(diffs) / sqrt(var_corr))
    p_val  <- ifelse(var_corr == 0, 1, 2 * pt(-abs(t_stat), df = K - 1))
    tibble(Baseline = "Lasso LR", Comparator = comp_model, Mean_AUC_Diff = round(mean(diffs), 5), Corrected_t = round(t_stat, 4), P_Value = round(p_val, 4), Significant_005 = p_val < 0.05)
  }
  t11_nb <- map_dfr(c("Elastic-Net LR", "Random Forest", "XGBoost", "SVM-RBF"), nadeau_bengio_test)
  save_table(t11_nb, "Table_11_Nadeau_Bengio_Corrected_Tests.csv")
  
  # 9. HOLDOUT EVALUATION, DELONG TESTS, CALIBRATION & DCA
  holdout_probs <- tibble(
    truth = test_data$outcome,
    `Lasso LR`              = predict(final_lasso, test_data, type = "prob")$.pred_Tumor,
    `Elastic-Net LR`        = predict(final_enet,  test_data, type = "prob")$.pred_Tumor,
    `Random Forest`         = predict(final_rf,    test_data, type = "prob")$.pred_Tumor,
    `XGBoost`               = predict(final_xgb,   test_data, type = "prob")$.pred_Tumor,
    `SVM-RBF`               = predict(final_svm,   test_data, type = "prob")$.pred_Tumor,
    `Stacked Meta-Ensemble` = predict(stacked_ensemble, test_data, type = "prob")$.pred_Tumor,
    `Deep Learning (torch)` = dl_probs
  )
  model_names <- colnames(holdout_probs)[-1]
  
  # TABLE 12: Holdout Performance
  t12_holdout <- map_dfr(model_names, function(m) {
    p <- holdout_probs[[m]]
    pred_f <- factor(ifelse(p >= 0.5, "Tumor", "Normal"), levels = c("Tumor", "Normal"))
    n_corr <- sum(pred_f == holdout_probs$truth)
    n_tot  <- nrow(holdout_probs)
    bt     <- binom.test(n_corr, n_tot)
    roc_o  <- roc(holdout_probs$truth, p, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE)
    ci_auc <- as.numeric(ci.auc(roc_o, method = "delong"))
    eval_d <- tibble(truth = holdout_probs$truth, pred = pred_f, prob = p)
    tibble(Model = m, Accuracy = round(n_corr / n_tot, 4), Acc_CI_Low = round(bt$conf.int[1], 4), Acc_CI_High = round(bt$conf.int[2], 4), Sensitivity = round(sens_vec(eval_d$truth, eval_d$pred, event_level = "first"), 4), Specificity = round(spec_vec(eval_d$truth, eval_d$pred, event_level = "first"), 4), F1_Score = round(f_meas_vec(eval_d$truth, eval_d$pred, event_level = "first"), 4), MCC = round(mcc_vec(eval_d$truth, eval_d$pred), 4), ROC_AUC = round(ci_auc[2], 4), AUC_CI_Low = round(ci_auc[1], 4), AUC_CI_High = round(ci_auc[3], 4), Brier_Score = round(brier_class_vec(eval_d$truth, eval_d$prob, event_level = "first"), 4))
  })
  save_table(t12_holdout, "Table_12_Holdout_Test_Performance.csv")
  
  # TABLE 13: DeLong ROC Comparison Tests
  roc_lasso <- roc(holdout_probs$truth, holdout_probs[["Lasso LR"]], levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE)
  t13_delong <- map_dfr(model_names[-1], function(m) {
    roc_m <- roc(holdout_probs$truth, holdout_probs[[m]], levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE)
    dt <- tryCatch(roc.test(roc_lasso, roc_m, method = "delong"), error = function(e) NULL)
    tibble(Model_1 = "Lasso LR", Model_2 = m, AUC_1 = round(as.numeric(auc(roc_lasso)), 4), AUC_2 = round(as.numeric(auc(roc_m)), 4), DeLong_Z = ifelse(is.null(dt) || is.nan(dt$statistic), NA, round(as.numeric(dt$statistic), 4)), P_Value = ifelse(is.null(dt) || is.nan(dt$p.value), NA, round(as.numeric(dt$p.value), 4)))
  })
  save_table(t13_delong, "Table_13_DeLong_ROC_Comparison_Tests.csv")
  
  # FIGURE 10: Holdout ROC Curves
  roc_curve_df <- map_dfr(model_names, function(m) { tibble(truth = holdout_probs$truth, prob = holdout_probs[[m]]) %>% roc_curve(truth, prob, event_level = "first") %>% mutate(Model = m) })
  p10_roc <- ggplot(roc_curve_df, aes(x = 1 - specificity, y = sensitivity, color = Model)) + geom_abline(linetype = "dashed", color = "grey60") + geom_path(linewidth = 1) + labs(title = "", subtitle = "", main = "", x = "1 - Specificity (False Positive Rate)", y = "Sensitivity (True Positive Rate)") + theme_pub
  save_plot(p10_roc, "Figure_10_Holdout_ROC_Curves", width = 8, height = 6)
  
  # FIGURE 11: Precision-Recall Curves
  pr_curve_df <- map_dfr(model_names, function(m) { tibble(truth = holdout_probs$truth, prob = holdout_probs[[m]]) %>% pr_curve(truth, prob, event_level = "first") %>% mutate(Model = m) })
  p11_pr <- ggplot(pr_curve_df, aes(x = recall, y = precision, color = Model)) + geom_path(linewidth = 1) + labs(title = "", subtitle = "", main = "", x = "Recall (Sensitivity)", y = "Precision (Positive Predictive Value)") + theme_pub
  save_plot(p11_pr, "Figure_11_Precision_Recall_Curves", width = 8, height = 6)
  
  # FIGURE 12: Confusion Matrices Grid
  cm_df <- map_dfr(model_names, function(m) { pred_f <- factor(ifelse(holdout_probs[[m]] >= 0.5, "Tumor", "Normal"), levels = c("Tumor", "Normal")); as.data.frame(table(Predicted = pred_f, Observed = holdout_probs$truth)) %>% mutate(Model = m) })
  p12_cm <- ggplot(cm_df, aes(x = Observed, y = Predicted, fill = Freq)) + geom_tile(color = "white") + geom_text(aes(label = Freq), fontface = "bold", size = 4) + scale_fill_gradient(low = "#F0F4F8", high = "#0072B2") + facet_wrap(~ Model, ncol = 4) + labs(title = "", subtitle = "", main = "") + theme_pub + theme(legend.position = "none")
  save_plot(p12_cm, "Figure_12_Holdout_Confusion_Matrices", width = 11, height = 6)
  
  # TABLE 14: Probability Calibration Metrics
  y_num <- ifelse(holdout_probs$truth == "Tumor", 1, 0)
  t14_cal <- map_dfr(model_names, function(m) {
    p <- pmin(pmax(holdout_probs[[m]], 1e-6), 1 - 1e-6)
    bins <- cut(p, breaks = seq(0, 1, length.out = 6), include.lowest = TRUE)
    bin_stats <- tibble(p = p, y = y_num, bin = bins) %>% group_by(bin) %>% summarise(n = n(), mean_p = mean(p), obs_y = mean(y), .groups = "drop")
    ece <- sum((bin_stats$n / length(p)) * abs(bin_stats$mean_p - bin_stats$obs_y))
    mce <- max(abs(bin_stats$mean_p - bin_stats$obs_y))
    spieg_num <- sum((y_num - p) * (1 - 2 * p))
    spieg_den <- sqrt(sum(((1 - 2 * p)^2) * p * (1 - p)))
    spieg_z   <- spieg_num / spieg_den
    tibble(Model = m, Brier_Score = round(mean((p - y_num)^2), 4), Expected_Cal_Error = round(ece, 4), Max_Cal_Error = round(mce, 4), Spiegelhalter_Z = round(spieg_z, 4), Spiegelhalter_P = round(2 * pnorm(-abs(spieg_z)), 4))
  })
  save_table(t14_cal, "Table_14_Probability_Calibration_Metrics.csv")
  
  # FIGURE 13: Calibration Reliability Diagrams (FIX 2: fixed probability bins immune to tied quantiles)
  cal_plot_df <- map_dfr(model_names, function(m) {
    p <- pmin(pmax(holdout_probs[[m]], 0), 1)
    bins <- cut(p, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)
    tibble(p = p, y = y_num, bin = bins, Model = m) %>%
      group_by(Model, bin) %>%
      summarise(Mean_Pred = mean(p), Obs_Prop = mean(y), .groups = "drop")
  })
  p13_cal <- ggplot(cal_plot_df, aes(x = Mean_Pred, y = Obs_Prop, color = Model)) + geom_abline(linetype = "dashed", color = "grey50") + geom_line(linewidth = 0.9) + geom_point(size = 2) + labs(title = "", subtitle = "", main = "", x = "Mean Predicted Probability (Tumor)", y = "Observed Tumor Proportion") + theme_pub
  save_plot(p13_cal, "Figure_13_Calibration_Reliability_Diagrams", width = 8, height = 6)
  
  # TABLE 15: Threshold Sensitivity & Net Benefit
  thresholds <- seq(0.10, 0.90, by = 0.10)
  N_hold <- length(y_num)
  t15_thresh <- map_dfr(c("Lasso LR", "Stacked Meta-Ensemble", "Deep Learning (torch)"), function(m) {
    p <- holdout_probs[[m]]
    map_dfr(thresholds, function(pt) {
      tp <- sum(p >= pt & y_num == 1); fp <- sum(p >= pt & y_num == 0); fn <- sum(p < pt & y_num == 1); tn <- sum(p < pt & y_num == 0)
      net_benefit <- (tp / N_hold) - (fp / N_hold) * (pt / (1 - pt))
      tibble(Model = m, Threshold = pt, Sensitivity = round(tp / (tp + fn), 4), Specificity = round(tn / (tn + fp), 4), False_Negatives = fn, False_Positives = fp, Net_Benefit = round(net_benefit, 4))
    })
  })
  save_table(t15_thresh, "Table_15_Threshold_Sensitivity_Net_Benefit.csv")
  
  # FIGURE 14: Decision Curve Analysis
  prev <- mean(y_num)
  treat_all_df <- tibble(Threshold = thresholds, Model = "Treat All", Net_Benefit = prev - (1 - prev) * (thresholds / (1 - thresholds)))
  dca_df <- bind_rows(t15_thresh %>% select(Threshold, Model, Net_Benefit), treat_all_df)
  p14_dca <- ggplot(dca_df, aes(x = Threshold, y = Net_Benefit, color = Model)) + geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") + geom_line(linewidth = 1) + coord_cartesian(ylim = c(-0.05, max(dca_df$Net_Benefit) + 0.05)) + labs(title = "", subtitle = "", main = "", x = "Probability Threshold (pt)", y = "Net Clinical Benefit") + theme_pub
  save_plot(p14_dca, "Figure_14_Decision_Curve_Analysis", width = 8, height = 6)
  
  # 10. FEATURE IMPORTANCE & THRESHOLD SENSITIVITY
  rf_imp <- extract_fit_parsnip(final_rf)$fit$variable.importance
  t16_imp <- tibble(Gene = names(rf_imp), RF_Gini_Importance = as.numeric(rf_imp)) %>% arrange(desc(RF_Gini_Importance)) %>% left_join(t9_lasso_or %>% select(Gene, Lasso_Beta = Beta_per_SD, Lasso_OR = OR_per_SD), by = "Gene")
  
  # TABLE 16: Consensus Biomarker Importance
  save_table(head(t16_imp, 30), "Table_16_Consensus_Biomarker_Importance.csv")
  
  # FIGURE 15: Top 20 Random Forest Gini Importance
  p15_imp <- head(t16_imp, 20) %>% ggplot(aes(x = RF_Gini_Importance, y = reorder(Gene, RF_Gini_Importance))) + geom_col(fill = "#0072B2", alpha = 0.85) + labs(title = "", subtitle = "", main = "", x = "Mean Decrease in Gini Impurity", y = NULL) + theme_pub
  save_plot(p15_imp, "Figure_15_Random_Forest_Feature_Importance", width = 8, height = 6)
  
  # FIGURE 16: Threshold Sensitivity Curves
  p16_thresh <- t15_thresh %>% pivot_longer(cols = c(Sensitivity, Specificity), names_to = "Metric", values_to = "Value") %>% ggplot(aes(x = Threshold, y = Value, color = Metric, linetype = Model)) + geom_line(linewidth = 1) + labs(title = "", subtitle = "", main = "", x = "Decision Threshold", y = "Score") + theme_pub
  save_plot(p16_thresh, "Figure_16_Threshold_Sensitivity_Curves", width = 9, height = 6)
  
  message("======================================================================")
  message(sprintf("PIPELINE COMPLETE: 16 Tables and 16 Figures saved in '%s'", out_dir))
  message("======================================================================")
})