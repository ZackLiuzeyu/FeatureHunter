#' Universal Feature Hunter across Multiple ML Models
#'
#' @description
#' Automatically detect the selected model type (MLP / RF / GLMNET / XGBoost / SVM / LDA / QDA / Naive Bayes)
#' from a leaderboard, train the corresponding model, and compute three types of feature importance:
#' \itemize{
#'   \item Model-based importance (weights, coefficients, gain, etc.)
#'   \item Permutation importance (performance drop by shuffling a feature)
#'   \item SHAP importance (model-agnostic explanation)
#' }
#' The three importance scores are fused via robust z-scores (median/MAD) into a composite score to rank features.
#' To reduce variance, multiple seeds are trained, predictions are ensembled first, then permutation/SHAP are computed
#' on the ensemble predictor; error bars use bootstrap CI over per-run composite scores.
#'
#' @param train_exp Numeric matrix of predictors (samples x features).
#' @param train_labels Response vector (0/1, factor, or convertible).
#' @param nshow Number of models shown in leaderboard.
#' @param namesS Character vector of metric names shown in leaderboard
#'   (default: \code{c("Accuracy","Recall","F-score")}).
#' @param score_index Integer index of metric to rank by (default: 3 = F-score).
#' @param pick_index Row index of leaderboard to pick (default: 1). Works for any model type.
#' @param top_models_csv Path to leaderboard CSV (auto-inferred if \code{NULL}).
#' @param num_runs Number of repeat runs (seeds) (default: 10).
#' @param num_coregene Number of core genes to keep in stability/importance plots.
#' @param n_likes Integer or "auto". Number of top genes to use in downstream analysis (UMAP, logistic). When "auto", UMAP/logistic will use the count of genes whose 95% CI does not cross 0 within the composite bar plot candidate set.
#' @param n_interest Number of genes to show in the stability heatmap.
#' @param seed Random seed (default: 424).
#' @param out_dir Output directory for plots (default: working dir).
#'
#' @param learningratei,batch_sizei,epochselecti,dropoutratei,cutoffi
#' Optional overrides for MLP hyperparameters (if \code{NULL}, use internal defaults
#' or values parsed from the leaderboard string where applicable).
#'
#' @param perm_metric Performance metric for permutation importance (\code{"f1"} or \code{"prauc"}).
#' @param perm_nrep Number of repetitions for permutation (default: 3).
#' @param shap_nsim Number of Monte Carlo simulations for SHAP (default: 10).
#' @param shap_subsample Number of samples for SHAP subsampling (default: 100).
#'
#' @param svm_cost Cost parameter for SVM (default: 1).
#' @param svm_gamma Gamma for SVM-RBF (default: \code{1/p} if \code{NULL}).
#'
#' @param xgb_nrounds,xgb_eta,xgb_max_depth Hyperparameters for XGBoost
#'   (defaults: \code{200}, \code{0.05}, \code{4}).
#'
#' @param rf_num_trees,rf_mtry Hyperparameters for Random Forest
#'   (defaults: \code{500} trees, \code{sqrt(p)} mtry if \code{NULL}).
#'
#' @param shuffle Logical; whether to shuffle data each epoch in MLP training (default: \code{TRUE}).
#'
#' @param standardize logical; if TRUE, standardize using train-set mean/sd and apply to val/test.
#'   Default: NULL (auto).
#' @param hidden_units integer vector; number of units per hidden layer,
#'   e.g. c(32L,16L,8L) or c(64L,32L,16L,8L). Default: NULL (auto).
#' @param activation character; activation function for hidden layers,
#'   e.g. "relu". Default: NULL (auto).
#' @param use_batchnorm logical; whether to insert BatchNorm after each hidden Dense.
#'   Default: NULL (auto).
#' @param l2 numeric; L2 regularization strength (e.g. 1e-4). Default: NULL (auto).
#' @param gaussian_noise_sd numeric; stddev for input GaussianNoise layer.
#'   Default: NULL (disabled).
#' @param min_lr numeric; minimum learning rate for ReduceLROnPlateau.
#'   Default: NULL (auto).
#' @param plateau_factor numeric; factor for ReduceLROnPlateau (e.g. 0.5).
#'   Default: NULL (auto).
#' @param plateau_patience integer; patience (epochs) for ReduceLROnPlateau.
#'   Default: NULL (auto).
#' @param early_patience integer; patience (epochs) for EarlyStopping.
#'   Default: NULL (auto).
#' @param imbalance_thresh numeric in \eqn{[0,1]}; threshold for enabling imbalance handling.
#'   Default: NULL (auto).
#' @param auto_th_method character; "youden", "f1", or "auto" for thresholding.
#'   Default: NULL (auto).
#'
#' @param method_weights length-3 named numeric vector for composite weighting,
#'   names must be c("model","perm","shap"). Default: c(model=1, perm=1, shap=1).
#' @param ci_level numeric; CI level for composite bootstrap over runs (default 0.95).
#' @param strict_min_freq numeric; minimum cross-run hit frequency in per-run Top-K to keep a gene
#'   (e.g. 0.6 to 0.7). NULL disables.
#' @param strict_min_effect numeric; minimum composite mean to keep a gene on the composite scale.
#'   NULL disables.
#' @param strict_ci_gate logical; if TRUE, require composite CI lower bound > 0.
#'   Default: FALSE.
#' @param strict_knee logical; if TRUE, apply a knee cutoff on sorted composite means.
#'   Default: FALSE.
#' @param strict_nmax integer; optional hard cap after knee (e.g. 10). NULL disables.
#' @param apply_selection_to_plots logical; if TRUE, bar/UMAP use selected genes;
#'   default FALSE keeps original plotting behavior.
#'
#' @return A list with elements:
#' \item{params}{List of parsed parameters and metadata}
#' \item{importance_df}{Data frame with per-gene importance scores (3 methods + composite)}
#' \item{top_list}{List of top genes per run}
#' \item{final_top}{Vector of selected top genes after selection controls}
#' \item{composite_mat}{Matrix of per-run composite scores}
#' \item{glm_summary}{Summary of logistic regression on top genes (coefficients + formula)}
#'
#' @details
#' Output files in \code{out_dir}:
#' \itemize{
#'   \item FI_Boxplot_Top20.pdf: Top-20 feature importance distributions
#'   \item FI_Bar_TopComposite.pdf: Top-N composite feature importance bar plot (mean with 95% CI)
#'   \item FI_Density_AllMethods.pdf: Distribution comparison of all importance methods
#'   \item UMAP_TopSignatureGenes.pdf: UMAP projection of top genes
#'   \item Stability_TopGenes_Heatmap.pdf: Stability heatmap of top genes across runs
#' }
#'
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv write.csv
#' @importFrom stats predict coef cov median mad quantile reorder sd glm
#' @importFrom graphics plot
#' @importFrom magrittr %>%
#' @export
fh_hunter <- function(
    train_exp,
    train_labels,
    nshow,
    namesS = c("Accuracy", "Recall", "F-score"),
    score_index = 3,
    pick_index = 1,
    top_models_csv = NULL,
    num_runs = 10,
    num_coregene,
    n_likes,
    n_interest,
    seed = 424,
    out_dir = getwd(),
    # manual overrides for MLP (simple)
    learningratei = NULL,
    batch_sizei = NULL,
    epochselecti = NULL,
    dropoutratei = NULL,
    cutoffi = NULL,
    # permutation / SHAP
    perm_metric = "f1",
    perm_nrep = 3L,
    shap_nsim = 10L,
    shap_subsample = 100L,
    # SVM / XGB / RF defaults
    svm_cost = 1,
    svm_gamma = NULL,
    xgb_nrounds = 200,
    xgb_eta = 0.05,
    xgb_max_depth = 4,
    rf_num_trees = 500,
    rf_mtry = NULL,
    # MLP trainer toggles (advanced; all NULL by default)
    shuffle = TRUE,
    standardize = NULL,
    hidden_units = NULL,
    activation = NULL,
    use_batchnorm = NULL,
    l2 = NULL,
    gaussian_noise_sd = NULL,
    min_lr = NULL,
    plateau_factor = NULL,
    plateau_patience = NULL,
    early_patience = NULL,
    imbalance_thresh = NULL,
    auto_th_method = NULL,
    # >>> NEW: selection controls and weighting (defaults preserve old behavior) <<<
    method_weights = c(model = 1, perm = 1, shap = 1),
    ci_level = 0.95,
    strict_min_freq = NULL,
    strict_min_effect = NULL,
    strict_ci_gate = FALSE,
    strict_knee = FALSE,
    strict_nmax = NULL,
    apply_selection_to_plots = FALSE) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  ## >>> NEW: holders for MLP training curves (saved later) <<<
  p_train_loss <- NULL
  p_train_acc <- NULL

  stopifnot(!missing(train_exp), !missing(train_labels))
  stopifnot(is.numeric(num_runs) && num_runs >= 1)
  stopifnot(!missing(num_coregene), !missing(n_likes), !missing(n_interest))

  perm_metric <- perm_metric[1]

  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  plots_dir <- file.path(out_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(seed)
  train_exp <- as.matrix(train_exp)
  n <- nrow(train_exp)
  p <- ncol(train_exp)
  train_y <- as.numeric(train_labels)
  if (any(!train_y %in% c(0, 1))) {
    if (is.factor(train_labels) || is.character(train_labels)) {
      ty <- tolower(as.character(train_labels))
      train_y <- ifelse(ty %in% c("1", "pos", "positive", "yes", "true", "case"), 1,
        ifelse(ty %in% c("0", "neg", "negative", "no", "false", "control"), 0, NA)
      )
      if (any(is.na(train_y))) stop("train_labels cannot be converted to 0/1.")
    } else {
      stop("train_labels must be 0/1 (or convertible).")
    }
  }

  if (is.null(top_models_csv)) {
    top_models_csv <- file.path("heatmap", paste0(nshow, "_top_", namesS[score_index], "_models.csv"))
  }
  if (!file.exists(top_models_csv)) stop("Leaderboard CSV not found: ", top_models_csv)
  Top_models <- utils::read.csv(top_models_csv, stringsAsFactors = FALSE)
  if (!("Model" %in% names(Top_models))) stop("Leaderboard CSV is missing column 'Model'.")

  pick_row <- pick_index
  model_str <- Top_models$Model[pick_row]
  if (is.na(model_str)) stop("Invalid pick index: ", pick_row)

  model_type <- .detect_model_type(model_str)
  if (model_type == "unknown") stop("Unknown model type: ", model_str)

  cutoff_p <- .cutoff_from_modelstr(model_str)
  parsed <- list(model_type = model_type, cutoff = cutoff_p)

  # --- Refactored: Parse params and build fit_params list ---
  if (model_type == "mlp") {
    lr_p <- .num1("lr\\s*:\\s*[0-9\\.]+", model_str)
    bs_p <- .num1("bs\\s*:\\s*[0-9]+", model_str)
    ep_p <- .num1("ep\\s*:\\s*[0-9]+", model_str)
    dr_p <- .num1("dropout\\s*:\\s*[0-9\\.]+", model_str)
    learningratei <- learningratei %||% lr_p
    batch_sizei <- batch_sizei %||% bs_p
    epochselecti <- epochselecti %||% ep_p
    dropoutratei <- dropoutratei %||% dr_p
    cutoffi <- cutoffi %||% cutoff_p
    if (any(is.na(c(learningratei, batch_sizei, epochselecti, dropoutratei, cutoffi)))) {
      stop("[MLP] failed to parse params; please provide learningratei/batch_sizei/epochselecti/dropoutratei/cutoffi")
    }
    parsed <- c(parsed, list(lr = learningratei, bs = batch_sizei, ep = epochselecti, dropout = dropoutratei))
  } else if (model_type == "rf") {
    mt_p <- .num1("mtry\\s*=\\s*[0-9]+", model_str)
    if (is.null(rf_mtry)) {
      if (!is.na(mt_p)) rf_mtry <- mt_p else rf_mtry <- floor(sqrt(p))
    }

    nt_p <- .num1("(ntree|num\\.trees)\\s*=\\s*[0-9]+", model_str)
    if (is.null(rf_num_trees)) {
      if (!is.na(nt_p)) rf_num_trees <- nt_p else rf_num_trees <- 500
    }
    parsed <- c(parsed, list(mtry = rf_mtry, num.trees = rf_num_trees))
  } else if (model_type == "glmnet") {
    s_low <- tolower(model_str)
    alpha_p <- .num1("alpha\\s*[:=]\\s*[0-9\\.]+", model_str)
    if (is.na(alpha_p)) {
      if (grepl("lasso", s_low)) alpha_p <- 1 else if (grepl("ridge|\\brr\\b", s_low)) alpha_p <- 0 else if (grepl("elastic\\s*net|\\benr\\b", s_low)) alpha_p <- 0.5
    }
    parsed <- c(parsed, list(alpha = alpha_p))
  } else if (model_type == "xgb") {
    md_p <- .num1("max_depth\\s*=\\s*[0-9]+", model_str)
    if (is.finite(md_p)) xgb_max_depth <- md_p
    et_p <- .num1("(eta|lr)\\s*[:=]\\s*[0-9\\.]+", model_str)
    if (is.finite(et_p)) xgb_eta <- et_p
    nr_p <- .num1("(nrounds|trees)\\s*=\\s*[0-9]+", model_str)
    if (is.finite(nr_p)) xgb_nrounds <- nr_p
    parsed <- c(parsed, list(max_depth = xgb_max_depth, eta = xgb_eta, nrounds = xgb_nrounds))
  } else if (startsWith(model_type, "svm")) {
    c_p <- .num1("C\\s*=\\s*[0-9\\.]+", model_str)
    if (is.finite(c_p)) svm_cost <- c_p
    g_p <- .num1("gamma\\s*=\\s*[0-9\\.]+", model_str)
    if (is.finite(g_p)) svm_gamma <- g_p
    if (is.null(svm_gamma)) svm_gamma <- 1 / ncol(train_exp)
    parsed <- c(parsed, list(cost = svm_cost, gamma = svm_gamma))
  }

  fit_params <- c(
    parsed,
    list(
      use_batchnorm = use_batchnorm,
      l2 = l2,
      gaussian_noise_sd = gaussian_noise_sd,
      activation = activation,
      min_lr = min_lr,
      plateau_factor = plateau_factor,
      plateau_patience = plateau_patience,
      early_patience = early_patience,
      shuffle = shuffle,
      imbalance_thresh = imbalance_thresh,
      hidden_units = hidden_units,
      nrounds = xgb_nrounds,
      eta = xgb_eta,
      max_depth = xgb_max_depth,
      cost = svm_cost,
      gamma = svm_gamma,
      num.trees = rf_num_trees,
      mtry = rf_mtry,
      alpha = parsed$alpha
    )
  )

  message("[hunter] pick=", pick_index, " | type=", model_type, " | parsed=", paste(names(parsed), parsed, collapse = ", "))

  dimn <- list(colnames(train_exp), NULL)
  A_model <- matrix(NA_real_, nrow = p, ncol = num_runs, dimnames = dimn)
  B_perm <- matrix(NA_real_, nrow = p, ncol = num_runs, dimnames = dimn)
  C_shap <- matrix(0, nrow = p, ncol = num_runs, dimnames = dimn)
  fits <- vector("list", length = num_runs)
  top_list <- vector("list", length = num_runs)
  message("Step 1/3: Computing model-based importance (multi-run with seeds)")
  pb <- utils::txtProgressBar(min = 0, max = num_runs, style = 3)

  for (run in seq_len(num_runs)) {
    # --- Refactored: Call training wrapper ---
    train_res <- .train_single_model(model_type, train_exp, train_y, seed + run, fit_params, run_id = run)
    fit <- train_res$model

    if (run == 1 && !is.null(train_res$plots_data)) {
      # Reconstruct ggplot objects from returned data
      pdata <- train_res$plots_data
      p_train_loss <- ggplot2::ggplot(pdata, ggplot2::aes(x = Epoch, y = Loss)) +
        ggplot2::geom_line(color = "firebrick") +
        ggplot2::labs(title = "Training Loss Curve", y = "Loss") +
        ggplot2::theme_minimal()

      if (!all(is.na(pdata$Accuracy))) {
        p_train_acc <- ggplot2::ggplot(pdata, ggplot2::aes(x = Epoch, y = Accuracy)) +
          ggplot2::geom_line(color = "dodgerblue") +
          ggplot2::labs(title = "Training Accuracy Curve", y = "Accuracy") +
          ggplot2::theme_minimal()
      }
    }

    fits[[run]] <- fit
    W_model <- .inner_importance(fit, train_exp, model_type)
    A_model[, run] <- W_model

    composite_run_tmp <- .zrob(A_model[, run])
    top_idx <- order(composite_run_tmp, decreasing = TRUE)[seq_len(min(num_coregene, length(composite_run_tmp)))]
    top_list[[run]] <- colnames(train_exp)[top_idx]

    utils::setTxtProgressBar(pb, run)
  }
  close(pb)


  message("Step 2/3: Computing permutation importance (feature shuffling)")
  # Ensemble-first evaluation for permutation and SHAP
  predict_ens <- function(X) .pred_proba_ensemble(fits, X, model_type)

  .pb <- function(n) utils::txtProgressBar(min = 0, max = n, style = 3)
  .bump <- function(pb, i) utils::setTxtProgressBar(pb, i)

  steps <- c("perm_ensemble", "shap_ensemble", "broadcast", "composite", "aggregate", "bootstrap")
  pb_blk <- .pb(length(steps))
  st <- 0

  ## 1) perm per-run
  st <- st + 1
  .bump(pb_blk, st)
  pb_perm <- utils::txtProgressBar(min = 0, max = num_runs, style = 3)
  for (r in seq_len(num_runs)) {
    # 针对第 r 个模型的预测函数
    predict_fun_r <- function(X) .pred_proba(fits[[r]], X, model_type)
    # 逐 run 计算 permutation importance（不同 seed，避免完全一致）
    B_perm[, r] <- .perm_importance_fun(
      predict_fun = predict_fun_r,
      X = train_exp, y = train_y,
      metric = perm_metric, nrep = perm_nrep,
      stratified = TRUE, seed = seed + 1000 + r,
      progress = FALSE
    )
    utils::setTxtProgressBar(pb_perm, r)
  }
  close(pb_perm)

  message("Step 3/3: Computing SHAP importance (skipped for parametric models) (policy=", ifelse(.uses_shap_for_type(model_type), "enabled", "skipped"), ")")
  ## 2) shap per-run
  st <- st + 1
  .bump(pb_blk, st)
  pb_shap <- utils::txtProgressBar(min = 0, max = num_runs, style = 3)

  for (r in seq_len(num_runs)) {
    if (identical(model_type, "xgb")) {
      # XGBoost: use native TreeSHAP contributions
      contrib <- stats::predict(fits[[r]], newdata = train_exp, predcontrib = TRUE)
      if (ncol(contrib) == ncol(train_exp) + 1) {
        contrib <- contrib[, -ncol(contrib), drop = FALSE] # drop bias term
      }
      C_shap[, r] <- as.numeric(colMeans(abs(contrib), na.rm = TRUE))
      names(C_shap[, r]) <- colnames(train_exp)
    } else if (.uses_shap_for_type(model_type)) {
      # nonparametric/black-box models: approximate SHAP via fastshap
      predict_fun_r <- function(X) .pred_proba(fits[[r]], X, model_type)
      C_shap[, r] <- .shap_importance_fun(
        predict_fun = predict_fun_r,
        X = train_exp,
        shap_nsim = shap_nsim,
        shap_subsample = shap_subsample,
        seed = seed + 2000 + r,
        progress = FALSE
      )
    } else {
      # parametric models: skip SHAP to save time (fill zeros to keep shapes)
      C_shap[, r] <- rep(0, ncol(train_exp))
      names(C_shap[, r]) <- colnames(train_exp)
    }

    utils::setTxtProgressBar(pb_shap, r)
  }
  close(pb_shap)
  close(pb_blk)

  ## --- Robust z-score per run (by column) ---
  zA_mat <- apply(A_model, 2, .zrob) # p x num_runs
  zB_mat <- apply(B_perm, 2, .zrob) # p x num_runs
  zC_mat <- apply(C_shap, 2, .zrob) # p x num_runs

  if (is.vector(zA_mat)) zA_mat <- matrix(zA_mat, nrow = p)
  if (is.vector(zB_mat)) zB_mat <- matrix(zB_mat, nrow = p)
  if (is.vector(zC_mat)) zC_mat <- matrix(zC_mat, nrow = p)

  dimnames(zA_mat) <- dimnames(A_model)
  dimnames(zB_mat) <- dimnames(B_perm)
  dimnames(zC_mat) <- dimnames(C_shap)

  ## --- Helper: renormalize weights over enabled components ---
  .renorm_weights <- function(w, enabled = c("model", "perm", "shap")) {
    w2 <- w
    disable <- setdiff(names(w), enabled)
    if (length(disable)) w2[disable] <- 0
    s <- sum(w2)
    if (s > 0) w2 <- w2 / s
    w2
  }

  ## --- Check and normalize weights ---
  w <- method_weights
  if (is.null(names(w)) || !all(sort(names(w)) == c("model", "perm", "shap"))) {
    stop("method_weights must be a named numeric vector with names c('model','perm','shap').")
  }

  ## --- Detect whether SHAP is effectively available this run ---
  ## Policy-free detection: if C_shap has any finite non-zero entry, treat as effective.
  shap_has_signal <- isTRUE(any(is.finite(C_shap) & (abs(C_shap) > 0)))
  shap_effective <- shap_has_signal

  ## --- Finalize weights: drop SHAP if ineffective, then renormalize ---
  if (!shap_effective) {
    w <- .renorm_weights(w, enabled = c("model", "perm"))
  } else {
    w <- .renorm_weights(w, enabled = c("model", "perm", "shap"))
  }

  ## --- Fuse per run (align columns r for zA/zB/zC) ---
  composite_mat <- matrix(NA_real_,
    nrow = p, ncol = num_runs,
    dimnames = list(colnames(train_exp), NULL)
  )

  for (r in seq_len(num_runs)) {
    zA <- zA_mat[, r]
    zB <- zB_mat[, r]
    zC <- if (shap_effective) zC_mat[, r] else rep(0, p) # avoid 0 * NA

    ## Gentle clipping to avoid single-run outliers
    cap <- 5
    zA <- pmax(pmin(zA, cap), -cap)
    zB <- pmax(pmin(zB, cap), -cap)
    zC <- pmax(pmin(zC, cap), -cap)

    composite_mat[, r] <- w["model"] * zA + w["perm"] * zB + w["shap"] * zC
  }

  ## --- Per-method summaries ---
  savg <- rowMeans(A_model, na.rm = TRUE)
  ssd <- apply(A_model, 1, stats::sd, na.rm = TRUE)
  perm_avg <- rowMeans(B_perm, na.rm = TRUE)
  perm_sd <- apply(B_perm, 1, stats::sd, na.rm = TRUE)

  if (shap_effective) {
    shap_avg <- rowMeans(C_shap, na.rm = TRUE)
    shap_sd <- apply(C_shap, 1, stats::sd, na.rm = TRUE)
  } else {
    shap_avg <- rep(0, p)
    shap_sd <- rep(0, p)
  }

  ## --- Composite mean (final ranking basis) ---
  composite_avg <- rowMeans(composite_mat, na.rm = TRUE)

  ## --- Bootstrap CI over runs for the composite ---
  alpha <- 1 - ci_level
  comp_ci <- t(vapply(
    seq_len(p),
    function(j) .boot_ci(composite_mat[j, ], R = 1000, alpha = alpha),
    numeric(3L)
  ))
  colnames(comp_ci) <- c("Composite_mean", "Composite_lo", "Composite_hi")

  ## --- Assemble importance table ---
  importance_df <- data.frame(
    Gene = colnames(train_exp),
    Model_mean = savg, Model_sd = ssd,
    Perm_mean = perm_avg, Perm_sd = perm_sd,
    SHAP_mean = shap_avg, SHAP_sd = shap_sd,
    Composite = composite_avg,
    Composite_lo = comp_ci[, "Composite_lo"],
    Composite_hi = comp_ci[, "Composite_hi"],
    stringsAsFactors = FALSE
  )

  f_imp <- file.path(out_dir, "Importance_Table.csv")
  utils::write.csv(importance_df, f_imp, row.names = FALSE)
  message(sprintf("[%s] Saved importance table to: %s", .ts(), f_imp))

  # --- Selection filters: apply strict_* to ALL plots ---
  selected_genes <- NULL
  rank_df <- importance_df[order(importance_df$Composite, decreasing = TRUE), ]
  rank_df$FreqProp <- 0
  # build a temporary frequency map from top_list (per-run Top-K)
  tmp_tab <- as.data.frame(table(unlist(lapply(top_list, function(genelist) head(genelist, n_interest)))))
  if (nrow(tmp_tab)) {
    colnames(tmp_tab) <- c("Gene", "Frequency")
    freq_map <- setNames(tmp_tab$Frequency / num_runs, as.character(tmp_tab$Gene))
    rank_df$FreqProp <- ifelse(rank_df$Gene %in% names(freq_map), freq_map[rank_df$Gene], 0)
  }
  keep <- rep(TRUE, nrow(rank_df))
  if (is.numeric(strict_min_freq) && is.finite(strict_min_freq)) keep <- keep & (rank_df$FreqProp >= strict_min_freq)
  if (is.numeric(strict_min_effect) && is.finite(strict_min_effect)) keep <- keep & (rank_df$Composite >= strict_min_effect)
  if (isTRUE(strict_ci_gate)) keep <- keep & (rank_df$Composite_lo > 0)
  cand_df <- rank_df[keep, , drop = FALSE]
  if (isTRUE(strict_knee) && nrow(cand_df) >= 3) {
    k <- .find_knee(cand_df$Composite) # Refactored helper call
    cand_df <- cand_df[seq_len(k), , drop = FALSE]
  }
  if (is.numeric(strict_nmax) && is.finite(strict_nmax) && strict_nmax > 0) {
    cand_df <- cand_df[seq_len(min(nrow(cand_df), as.integer(strict_nmax))), , drop = FALSE]
  }
  if (!nrow(cand_df)) {
    cand_df <- rank_df[seq_len(min(nrow(rank_df), n_likes)), , drop = FALSE]
  }
  selected_genes <- cand_df$Gene
  plot_set <- if (isTRUE(apply_selection_to_plots)) selected_genes else colnames(train_exp)

  # Determine effective n_likes for downstream (UMAP/logistic)
  if (is.character(n_likes) && tolower(n_likes) == "auto") {
    n_likes_eff <- sum(importance_df$Gene %in% plot_set & importance_df$Composite_lo > 0, na.rm = TRUE)
    n_likes_eff <- max(2L, as.integer(n_likes_eff)) # keep at least 2 for UMAP
    message(sprintf(
      "[fh_hunter] n_likes='auto' -> using %d genes with %.0f%% CI > 0 within plot_set",
      n_likes_eff, ci_level * 100
    ))
  } else {
    n_likes_eff <- as.integer(n_likes)
  }

  # Top lists summary (respect plot_set)
  top_df <- as.data.frame(table(unlist(lapply(top_list, function(genelist) {
    head(intersect(genelist, plot_set), n_interest)
  }))))
  colnames(top_df) <- c("Gene", "Frequency")
  if (nrow(top_df)) {
    top_df$Gene <- factor(top_df$Gene, levels = top_df$Gene[order(top_df$Frequency, decreasing = TRUE)])
  }

  f_top <- file.path(out_dir, "FI_TopGeneFrequencies.csv")


  ## >>> NEW: save MLP training curves at the end (only if MLP and plots exist) <<<
  if (identical(model_type, "mlp") && (!is.null(p_train_loss) || !is.null(p_train_acc))) {
    f_loss <- file.path(plots_dir, "Training_Loss.pdf")
    f_acc <- file.path(plots_dir, "Training_Accuracy.pdf")

    if (!is.null(p_train_loss)) {
      ggplot2::ggsave(f_loss, p_train_loss, width = 6, height = 4)
      message(sprintf("[%s] Saved plot to: %s", .ts(), f_loss))
    }
    if (!is.null(p_train_acc)) {
      ggplot2::ggsave(f_acc, p_train_acc, width = 6, height = 4)
      message(sprintf("[%s] Saved plot to: %s", .ts(), f_acc))
    }
  }
  ## Decide whether SHAP was computed for this model_type
  include_shap <- isTRUE(.uses_shap_for_type(model_type)) # Refactored helper call

  ## Boxplot (Top20) across runs
  topN <- min(20, nrow(importance_df))
  topN_genes <- importance_df |>
    dplyr::arrange(dplyr::desc(Composite)) |>
    dplyr::pull(Gene)
  topN_genes <- intersect(topN_genes, plot_set)
  topN_genes <- head(topN_genes, topN)

  if (length(topN_genes) > 0) {
    df_list_box <- list(
      tibble::tibble(
        Gene = rep(topN_genes, each = num_runs),
        Method = "Model",
        Value = as.vector(A_model[topN_genes, , drop = FALSE])
      ),
      tibble::tibble(
        Gene = rep(topN_genes, each = num_runs),
        Method = "Permutation",
        Value = as.vector(B_perm[topN_genes, , drop = FALSE])
      )
    )
    if (.uses_shap_for_type(model_type)) {
      df_list_box[[length(df_list_box) + 1]] <- tibble::tibble(
        Gene = rep(topN_genes, each = num_runs),
        Method = "SHAP",
        Value = as.vector(C_shap[topN_genes, , drop = FALSE])
      )
    }
    long_df <- dplyr::bind_rows(df_list_box)

    p_box <- ggplot2::ggplot(long_df, ggplot2::aes(x = Gene, y = Value, fill = Method)) +
      ggplot2::geom_boxplot() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0("Top ", topN, " genes importance distributions"),
        y = "Importance (per-run)", x = "Gene"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    f_box <- file.path(plots_dir, "FI_Boxplot_Top20.pdf")
    ggplot2::ggsave(f_box, p_box, width = 9, height = 7)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_box))
  }

  # Bar (Top Composite) with 95% CI
  top_genes_plot <- importance_df |>
    dplyr::arrange(dplyr::desc(Composite))
  top_genes_plot <- dplyr::filter(top_genes_plot, Gene %in% plot_set)
  top_genes_plot <- dplyr::slice_head(top_genes_plot, n = num_coregene)
  if (nrow(top_genes_plot)) {
    p_bar <- ggplot2::ggplot(
      top_genes_plot,
      ggplot2::aes(x = stats::reorder(Gene, Composite), y = Composite)
    ) +
      ggplot2::geom_col(fill = "#0072B2", alpha = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = Composite_lo, ymax = Composite_hi),
        width = 0.4, color = "#D55E00"
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0(num_coregene, " Genes by Composite Importance"),
        x = "Gene", y = "Composite (mean with 95% CI)"
      ) +
      ggplot2::theme_minimal()
    f_bar <- file.path(plots_dir, "FI_Bar_TopComposite.pdf")
    ggplot2::ggsave(f_bar, p_bar, width = 8, height = 6)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_bar))
  }

  # Density (all)
  # Build combined importance long df, dropping SHAP if not used
  include_shap <- isTRUE(.uses_shap_for_type(model_type))
  A_use <- A_model[plot_set, , drop = FALSE]
  B_use <- B_perm[plot_set, , drop = FALSE]

  df_list <- list(
    data.frame(Method = "Model", Value = as.vector(A_use)),
    data.frame(Method = "Permutation", Value = as.vector(B_use))
  )
  if (include_shap) {
    C_use <- C_shap[plot_set, , drop = FALSE]
    df_list[[length(df_list) + 1]] <- data.frame(Method = "SHAP", Value = as.vector(C_use))
  }
  combined_imp <- dplyr::bind_rows(df_list)

  p_den <- ggplot2::ggplot(combined_imp, ggplot2::aes(x = Value, fill = Method)) +
    ggplot2::geom_density(alpha = 0.6, color = NA) +
    ggplot2::geom_vline(
      data = combined_imp %>%
        dplyr::group_by(Method) %>%
        dplyr::summarise(Median = stats::median(Value, na.rm = TRUE)),
      ggplot2::aes(xintercept = Median, color = Method), linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~Method, scales = "free") +
    ggplot2::labs(
      title = "Feature Importance Distribution Comparison",
      x = "Importance Score", y = "Density"
    ) +
    ggplot2::theme_bw()
  f_den <- file.path(plots_dir, "FI_Density_AllMethods.pdf")
  ggplot2::ggsave(f_den, p_den, width = 9, height = 6)
  message(sprintf("[%s] Saved plot to: %s", .ts(), f_den))

  # UMAP (consistent with plot_set)
  final_importance <- importance_df |>
    dplyr::arrange(dplyr::desc(Composite))
  final_top <- intersect(final_importance$Gene, plot_set)
  final_top <- head(final_top, n_likes_eff)
  set.seed(seed)
  if (length(final_top) >= 2) {
    um <- umap::umap(train_exp[, final_top, drop = FALSE])
    umap_df <- data.frame(
      UMAP1 = um$layout[, 1],
      UMAP2 = um$layout[, 2],
      Label = factor(train_y, labels = c("Control", "Case"))
    )
    p_umap <- ggplot2::ggplot(umap_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = Label)) +
      ggplot2::geom_point(alpha = 0.7, size = 2.5) +
      ggplot2::scale_color_manual(values = c("#377eb8", "#e41a1c")) +
      ggplot2::labs(title = "UMAP Projection using top Signature Genes", color = "Group") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    f_umap <- file.path(plots_dir, "UMAP_TopSignatureGenes.pdf")
    ggplot2::ggsave(f_umap, p_umap, width = 6, height = 5)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_umap))
  }

  # Stability heatmap (already aligned via plot_set-aware top_df)
  if (nrow(top_df)) {
    p_stab <- ggplot2::ggplot(top_df, ggplot2::aes(x = Gene, y = "Runs", fill = Frequency)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = Frequency), color = "black", size = 4) +
      ggplot2::scale_fill_gradient(low = "#E6F5FF", high = "#1F77B4", name = "Frequency") +
      ggplot2::labs(
        title = paste0("Top Gene Stability Across ", num_runs, " Runs"),
        x = "Gene", y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    f_stab <- file.path(plots_dir, "Stability_TopGenes_Heatmap.pdf")
    ggplot2::ggsave(f_stab, p_stab, width = 9, height = 6)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_stab))
  }


  ### formula construction for top n_likes_eff genes
  feat_names <- final_top[1:min(n_likes_eff, length(final_top))]
  x_all <- as.matrix(train_exp[, feat_names, drop = FALSE])
  y_all <- as.numeric(train_y)

  # --- Refactored: Formula Report Generation ---
  report_res <- .generate_formula_report(model_type, x_all, y_all, parsed, svm_cost)
  formula_source <- report_res$formula_source
  formula_str <- report_res$formula_str
  coef_summary <- report_res$coef_summary

  # ---- persist formula to text (snapshot + rolling log) ----
  ts_disp <- if (exists(".ts")) .ts() else format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ts_file <- gsub("[: ]", "-", ts_disp)

  file_ts <- file.path(plots_dir, sprintf("fh_formula_%s.txt", ts_file))
  log_file <- file.path(plots_dir, "fh_formula_log.txt")

  lines <- c(
    sprintf("[fh_hunter] timestamp: %s", ts_disp),
    sprintf("[fh_hunter] model_type: %s", if (exists("model_type")) as.character(model_type) else NA_character_),
    sprintf("[fh_hunter] formula_source: %s", if (exists("formula_source")) as.character(formula_source) else NA_character_),
    sprintf("[fh_hunter] n_features: %s", if (exists("coef_summary")) nrow(coef_summary) - 1L else NA_integer_),
    sprintf("[fh_hunter] formula: %s", if (exists("formula_str")) formula_str else NA_character_)
  )

  writeLines(lines, file_ts) # snapshot file
  write(paste(lines, collapse = "\n"), file = log_file, append = TRUE)
  write("\n---\n", file = log_file, append = TRUE)

  message(sprintf("[%s] Saved formula to: %s", if (exists(".ts")) .ts() else ts_disp, file_ts))
  message(sprintf("[%s] Appended formula to: %s", if (exists(".ts")) .ts() else ts_disp, log_file))

  # ---- wrap up and return ----
  invisible(list(
    params = list(
      top_models_csv = top_models_csv,
      score_name = namesS[score_index],
      picked_row = pick_index,
      parsed = parsed,
      seed = seed,
      out_dir = out_dir,
      method_weights = method_weights,
      ci_level = ci_level,
      strict_min_freq = strict_min_freq,
      strict_min_effect = strict_min_effect,
      strict_ci_gate = strict_ci_gate,
      strict_knee = strict_knee,
      strict_nmax = strict_nmax,
      apply_selection_to_plots = apply_selection_to_plots,
      formula_source = formula_source
    ),
    importance_df = importance_df,
    top_list = top_list,
    final_top = final_top,
    composite_mat = composite_mat,
    glm_summary = list(
      coefficients = coef_summary,
      formula = formula_str
    )
  ))
}
