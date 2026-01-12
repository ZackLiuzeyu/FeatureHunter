#' Gradient Boosting (GBM) – default settings with optional auto-thresholding
#'
#' Trains a Bernoulli GBM on the training set, selects decision thresholds on the
#' validation set B (fixed + optional auto threshold), evaluates on A/B/C/…,
#' and writes results to an optional collector environment as well as returning
#' local accumulators invisibly.
#'
#' @param train_exp data.frame/matrix. Training features (samples × features).
#' @param train_labels vector. Binary training labels (0/1 or convertible by `as_binary01`).
#' @param test_exp list. External datasets (first = validation B, then C/D/…).
#' @param labels_list list. Label data.frames for A=1, B=2, … (first column = ID, second = label).
#' @param cutoff numeric. Fixed thresholds merged with auto threshold from B. Default: `c(0.25, 0.5, 0.75)`.
#' @param auto_th_method character. Auto threshold method: `"youden"` / `"f1"` / `"auto"`.
#'   With `"auto"`, `decide_threshold_method()` is used based on class imbalance
#'   and PR-vs-ROC performance; otherwise the specified method is used. Default: `"auto"`.
#' @param auto_imbalance_thresh numeric. Imbalance threshold passed to `decide_threshold_method()`. Default: `0.35`.
#' @param auto_pr_vs_roc_gate numeric. PR-vs-ROC advantage ratio gate for `decide_threshold_method()`. Default: `0.5`.
#' @param collector environment (optional). If supplied, each cutoff’s results are written into it via `.fh_collect()`.
#'
#' @return
#' Invisibly returns a list with four elements:
#' \itemize{
#'   \item `acc`     – named vectors of accuracy for A/B/C/…
#'   \item `recall`  – named vectors of recall for A/B/C/…
#'   \item `fs`      – named vectors of F1 for A/B/C/…
#'   \item `summary` – per-dataset prediction data.frames
#' }
#' If `collector` is provided, results are also appended into it for downstream use.
#'
#' @examples
#' \dontrun{
#' fh_gbm_default(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   auto_th_method = "auto",
#'   auto_imbalance_thresh = 0.35,
#'   auto_pr_vs_roc_gate   = 0.5,
#'   collector    = collector
#' )
#' }
#' @export
fh_gbm_default <- function(train_exp, train_labels, test_exp, labels_list,
                           cutoff = c(0.25, 0.5, 0.75),
                           auto_th_method = "auto",
                           auto_imbalance_thresh = 0.35,
                           auto_pr_vs_roc_gate   = 0.5,
                           collector = collector) {
  
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # ---- train GBM (Bernoulli) ----
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- as.numeric(as_binary01(train_labels))  # ensure 0/1
  model <- gbm::gbm(
    labels ~ . - labels,
    data = train_expd,
    distribution = "bernoulli",
    bag.fraction = 0.8,
    n.minobsinnode = 10,
    verbose = FALSE                     # ← 静默训练
  )
  
  # ← 只计算一次最佳树数，并静默 gbm.perf 的输出
  best_trees <- tryCatch({
    suppressMessages(gbm::gbm.perf(model, method = "OOB", plot.it = FALSE))
  }, error = function(e) {
    if (!is.null(model$n.trees)) model$n.trees else 100L
  })
  if (!is.finite(best_trees) || best_trees <= 0) {
    best_trees <- if (!is.null(model$n.trees)) model$n.trees else 100L
  }
  
  # ---- auto threshold on B (merged with fixed cutoffs) ----
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA_real_
  method_used <- auto_th_method
  
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    prob_v <- as.numeric(gbm::predict.gbm(
      model, newdata = val_x, type = "response", n.trees = best_trees  # ← 指定 n.trees
    ))
    
    if (length(unique(prob_v[is.finite(prob_v)])) >= 2) {
      y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(
          probs = prob_v, y_true = y_val,
          imbalance_thresh = auto_imbalance_thresh,
          pr_vs_roc_gate   = auto_pr_vs_roc_gate
        )
      }
      
      th_auto <- choose_threshold(prob_v, y_val, method = method_used)
      
      if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
        message(sprintf("[GBM] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
      } else {
        message("[GBM] Auto threshold skipped; keep original cutoff grid.")
      }
    } else {
      message("[GBM] Insufficient probability variation on B; keep original cutoff grid.")
    }
  } else {
    message("[GBM] Validation overlap < 2; keep original cutoff grid.")
  }
  
  # ---- evaluate over thresholds ----
  for (cutoffi in cutoff_now) {
    # Train (A)
    prob_tr <- as.numeric(gbm::predict.gbm(
      model, newdata = as.data.frame(train_exp), type = "response", n.trees = best_trees  # ← 指定
    ))
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cutoffi, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
    
    # External sets B/C/D/…
    all_result <- lapply(seq_along(test_exp), function(i) {
      data_i   <- as.data.frame(test_exp[[i]])
      label_df <- labels_list[[i + 1]]
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) {
        return(data.frame(
          predict_p = numeric(0),
          predict_result = factor(character()),
          real_label = factor(character())
        ))
      }
      expdata <- data_i[comd, , drop = FALSE]
      p  <- as.numeric(gbm::predict.gbm(
        model, newdata = expdata, type = "response", n.trees = best_trees  # ← 指定
      ))
      pr <- factor(ifelse(p > cutoffi, "positive", "negative"))
      out <- data.frame(predict_p = p, predict_result = pr)
      out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1, "positive", "negative"))
      out
    })
    
    all_result <- c(list(Train = train_result), all_result)
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x) {
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 &&
      (abs(cutoffi - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("GBM-default (cutoff:auto(",
             formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
    } else {
      paste0("GBM-default (cutoff:", formatC(cutoffi, format = "f", digits = 4), ")")
    }
    
    if (is.environment(collector)) {
      .fh_collect(
        collector, key,
        acc     = result_acc,
        recall  = result_recall,
        fs      = result_FS,
        summary = all_result
      )
    }
    
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  invisible(list(
    acc     = acc_local,
    recall  = recall_local,
    fs      = fs_local,
    summary = summary_local
  ))
}






#' GBM (caret CV-tuned) interface
#'
#' @description
#' Cross-validated hyperparameter tuning for GBM via \code{caret}, then refit a
#' single GBM model with \code{bestTune}. On validation set B, automatically
#' chooses a probability cutoff (F1/Youden or auto), merges it with fixed
#' thresholds, evaluates on A/B/C/..., and writes results into the optional
#' \code{collector} environment and also returns them invisibly.
#'
#' @param train_exp data.frame; training expression matrix (samples x features).
#' @param train_labels vector; binary training labels (0/1 or convertible).
#' @param test_exp list; external datasets (first = validation B, then C/D/...).
#' @param labels_list list; list of label data.frames (A=1, B=2, …; col1=ID, col2=label).
#' @param fold integer; CV folds for \code{caret::train} (default 10).
#' @param cores integer; number of parallel workers. If \code{NULL}, use \code{detectCores()-1}.
#' @param cutoff numeric; fixed cutoffs merged with the auto-threshold from B (default \code{c(0.25,0.5,0.75)}).
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to pick the strategy
#'   based on class imbalance and PR vs ROC performance; otherwise the specified method is used.
#' @param auto_imbalance_thresh Imbalance threshold for \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector (optional) environment collector; if supplied, each cutoff’s results are written into it.
#'
#' @return
#' Invisibly returns a list with \code{acc}, \code{recall}, \code{fs}, \code{summary},
#' and \code{collector}. If \code{collector} is an environment, results are also
#' appended into it via \code{.fh_collect()}.
#'
#' @examples
#' \dontrun{
#' # prepare a collector env if you want to aggregate results
#' collector <- new.env(parent = emptyenv())
#' out <- fh_gbm_cv_best(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   fold         = 10,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   auto_th_method = "auto",
#'   auto_imbalance_thresh = 0.35,
#'   auto_pr_vs_roc_gate   = 0.5,
#'   collector    = collector
#' )
#' }
#' @export
fh_gbm_cv_best <- function(train_exp, train_labels, test_exp, labels_list,
                           fold = 10, cores = NULL,
                           cutoff = c(0.25, 0.5, 0.75),
                           auto_th_method = "auto",
                           auto_imbalance_thresh = 0.35,
                           auto_pr_vs_roc_gate   = 0.5,
                           collector = collector) {
  # local accumulators to return invisibly
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # --- parallel workers
  if (is.null(cores)) {
    cores <- tryCatch(max(1L, parallel::detectCores() - 1L), error = function(e) 1L)
  }
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)
  
  # 1) caret tuning (Bernoulli)
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- as.numeric(as_binary01(train_labels))
  train_expd$labels <- factor(train_expd$labels, levels = c(0, 1))
  
  gbmGrid <- expand.grid(
    interaction.depth = c(3, 4, 5),
    n.trees           = c(seq(50, 300, 50), seq(400, 1000, 200)),
    shrinkage         = c(0.005, 0.01, 0.02, 0.05),
    n.minobsinnode    = c(5, 10, 15)
  )
  
  fitControl <- caret::trainControl(
    method = "cv",
    number = fold,
    verboseIter = FALSE          # ← 关闭 caret 训练日志
  )
  
  GBMberT <- caret::train(
    labels ~ . - labels,
    data        = train_expd,
    distribution = "bernoulli",
    method       = "gbm",
    trControl    = fitControl,
    tuneGrid     = gbmGrid,
    bag.fraction = 0.8,
    verbose      = FALSE         # ← 传给底层 gbm，关闭迭代表格
  )
  
  # 2) refit a single GBM with bestTune
  besttune <- GBMberT$bestTune
  train_expd2 <- as.data.frame(train_exp)
  train_expd2$labels <- as.numeric(as_binary01(train_labels))
  model <- gbm::gbm(
    labels ~ . - labels,
    data              = train_expd2,
    distribution      = "bernoulli",
    bag.fraction      = 0.8,
    n.minobsinnode    = besttune$n.minobsinnode,
    shrinkage         = besttune$shrinkage,
    n.trees           = besttune$n.trees,
    interaction.depth = besttune$interaction.depth,
    verbose           = FALSE     # ← 重训也关闭输出
  )
  
  # 3) fixed + auto threshold on B
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])  # B
  val_lab <- labels_list[[2]]              # labels for B
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA_real_; method_used <- auto_th_method
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    prob_v <- as.numeric(gbm::predict.gbm(
      model, newdata = val_x, type = "response", n.trees = besttune$n.trees
    ))
    
    if (length(unique(prob_v[is.finite(prob_v)])) >= 2) {
      y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(
          probs = prob_v, y_true = y_val,
          imbalance_thresh = auto_imbalance_thresh,
          pr_vs_roc_gate   = auto_pr_vs_roc_gate
        )
      }
      th_auto <- choose_threshold(prob_v, y_val, method = method_used)
      
      if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
        message(sprintf("[GBM CV best] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
      } else {
        message("[GBM CV best] Auto threshold skipped; keep original cutoff grid.")
      }
    } else {
      message("[GBM CV best] B has insufficient probability variation; keep original cutoff grid.")
    }
  } else {
    message("[GBM CV best] Validation overlap < 2; keep original cutoff grid.")
  }
  
  # 4) evaluate for each cutoff
  for (cutoffi in cutoff_now) {
    prob_tr <- as.numeric(gbm::predict.gbm(
      model, newdata = as.data.frame(train_exp), type = "response", n.trees = besttune$n.trees
    ))
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cutoffi, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
    
    all_result <- lapply(seq_along(test_exp), function(i) {
      data_i   <- as.data.frame(test_exp[[i]])
      label_df <- labels_list[[i + 1]]
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) {
        return(data.frame(
          predict_p = numeric(0),
          predict_result = factor(character()),
          real_label = factor(character())
        ))
      }
      expdata <- data_i[comd, , drop = FALSE]
      p <- as.numeric(gbm::predict.gbm(
        model, newdata = expdata, type = "response", n.trees = besttune$n.trees
      ))
      pr <- factor(ifelse(p > cutoffi, "positive", "negative"))
      out <- data.frame(predict_p = p, predict_result = pr)
      out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1, "positive", "negative"))
      out
    })
    
    all_result <- c(list(Train = train_result), all_result)
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x) {
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 &&
      (abs(cutoffi - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("GBM-CV:", fold, " fold (cutoff:auto(",
             formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
    } else {
      paste0("GBM-CV:", fold, " fold (cutoff:",
             formatC(cutoffi, format = "f", digits = 4), ")")
    }
    
    if (is.environment(collector)) {
      .fh_collect(
        collector, key,
        acc     = result_acc,
        recall  = result_recall,
        fs      = result_FS,
        summary = all_result
      )
    }
    
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  invisible(list(
    acc     = acc_local,
    recall  = recall_local,
    fs      = fs_local,
    summary = summary_local,
    collector = collector
  ))
}


#' LASSO + GBM (Default Parameters)
#'
#' Fit a Gradient Boosting Machine (GBM, Bernoulli) on LASSO-selected features.
#' On validation set B, auto-pick a probability threshold (F1/Youden or auto),
#' merge it with fixed thresholds, then evaluate on A/B/C/… and record results.
#'
#' @param train_exp Training expression matrix (rows = samples, columns = genes).
#' @param train_labels Training labels (convertible to 0/1).
#' @param test_exp A list of external sets: 1st = validation set B, others = C/D/….
#' @param labels_list A list of label tables for datasets (A=1, B=2, …; col1 = ID, col2 = label).
#' @param lassogene Character vector of genes selected by LASSO.
#' @param fold Integer; used only in the result key for naming consistency (default 10).
#' @param cutoff Fixed cutoff thresholds; default \code{c(0.25, 0.5, 0.75)}.
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to pick the strategy
#'   based on class imbalance and PR-vs-ROC performance; otherwise the specified method is used.
#' @param auto_imbalance_thresh Imbalance threshold for \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector (optional) environment collector; if supplied, each cutoff’s
#'   results are written into it via \code{.fh_collect()}.
#'
#' @return \code{invisible(list(acc=..., recall=..., fs=..., summary=...))}.
#'   If \code{collector} is provided, results are also written into that environment.
#'
#' @seealso \code{\link{as_binary01}}, \code{\link{decide_threshold_method}},
#'   \code{\link{choose_threshold}}, \code{\link{f1_score}}
#'
#' @examples
#' \dontrun{
#' collector <- new.env(parent = emptyenv())
#' fh_gbm_lasso_default(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   fold         = 10,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   auto_th_method = "auto",
#'   collector    = collector
#' )
#' }
#' @export
fh_gbm_lasso_default <- function(
    train_exp,
    train_labels,
    test_exp,
    labels_list,
    lassogene,
    fold = 10,
    cutoff = c(0.25, 0.5, 0.75),
    auto_th_method = "auto",
    auto_imbalance_thresh = 0.35,
    auto_pr_vs_roc_gate   = 0.5,
    collector = collector
) {
  # local accumulators to return invisibly
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # 1) Fit GBM on LASSO-selected features
  train_expd <- as.data.frame(train_exp)[, lassogene, drop = FALSE]
  train_expd$labels <- as.numeric(as_binary01(train_labels))  # ensure 0/1
  model <- gbm::gbm(
    labels ~ . - labels,
    data = train_expd,
    distribution = "bernoulli",
    bag.fraction = 0.8,
    n.minobsinnode = 10
  )
  n_trees <- model$n.trees  # <<< 用于静默预测
  
  # 2) Auto threshold on B and merge with fixed thresholds
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA_real_
  method_used <- auto_th_method
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    prob_v <- as.numeric(gbm::predict.gbm(model, newdata = val_x, type = "response", n.trees = n_trees))
    if (length(unique(prob_v[is.finite(prob_v)])) >= 2) {
      y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(
          probs = prob_v, y_true = y_val,
          imbalance_thresh = auto_imbalance_thresh,
          pr_vs_roc_gate   = auto_pr_vs_roc_gate
        )
      }
      th_auto <- choose_threshold(prob_v, y_val, method = method_used)
      
      if (is.finite(th_auto) && length(th_auto) == 1L && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
        message(sprintf("[lasso+GBM] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
      } else {
        message("[lasso+GBM] Auto threshold skipped; keep original cutoff grid.")
      }
    } else {
      message("[lasso+GBM] B has insufficient probability variation; keep original cutoff grid.")
    }
  } else {
    message("[lasso+GBM] Validation overlap < 2; keep original cutoff grid.")
  }
  
  # 3) Evaluate across cutoffs
  for (cutoffi in cutoff_now) {
    # Train (A)
    prob_tr <- as.numeric(gbm::predict.gbm(
      model, newdata = as.data.frame(train_exp)[, lassogene, drop = FALSE],
      type = "response", n.trees = n_trees
    ))
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cutoffi, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
    
    # External sets: B/C/D/…
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])[, lassogene, drop = FALSE]
      label_df <- labels_list[[i + 1]]
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) {
        return(data.frame(
          predict_p = numeric(0),
          predict_result = factor(character()),
          real_label = factor(character())
        ))
      }
      expdata <- data_i[comd, , drop = FALSE]
      p  <- as.numeric(gbm::predict.gbm(model, newdata = expdata, type = "response", n.trees = n_trees))
      pr <- factor(ifelse(p > cutoffi, "positive", "negative"))
      out <- data.frame(predict_p = p, predict_result = pr)
      out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1, "positive", "negative"))
      out
    })
    
    # prepend Train
    all_result <- c(list(Train = train_result), all_result)
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    # metrics
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x) {
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    # key
    is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 &&
      (abs(cutoffi - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("GBM-default+Lasso-CV:", fold, " fold (cutoff:auto(",
             formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
    } else {
      paste0("GBM-default+Lasso-CV:", fold, " fold (cutoff:",
             formatC(cutoffi, format = "f", digits = 4), ")")
    }
    
    # write to collector if provided
    if (is.environment(collector)) {
      .fh_collect(
        collector, key,
        acc     = result_acc,
        recall  = result_recall,
        fs      = result_FS,
        summary = all_result
      )
    }
    
    # also accumulate locally
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  invisible(list(
    acc     = acc_local,
    recall  = recall_local,
    fs      = fs_local,
    summary = summary_local
  ))
}


#' LASSO + GBM (best via caret tuning)
#'
#' Fit a GBM (Bernoulli) on LASSO-selected features. Hyper-parameters are
#' tuned via caret CV; a probability cutoff is then chosen on validation set B
#' and merged with fixed cutoffs for evaluation on A/B/C/... datasets.
#'
#' @param train_exp Training expression matrix (rows = samples, columns = genes).
#' @param train_labels Training labels (convertible to 0/1).
#' @param test_exp List of test sets: the 1st element = validation set B, others = C/D/….
#' @param labels_list List of label tables (A=1, B=2, …; column 1 = sample ID, column 2 = label).
#' @param lassogene Character vector of genes selected by LASSO.
#' @param fold Integer; folds for caret cross-validation (default 10).
#' @param cutoff Fixed cutoff thresholds; default \code{c(0.25, 0.5, 0.75)}.
#' @param cores Number of CPU cores for parallel tuning with caret/doParallel.
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} will pick the strategy.
#' @param auto_imbalance_thresh Imbalance threshold for \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector Optional environment; if supplied, each cutoff’s results are written into it
#'   via \code{.fh_collect()}.
#'
#' @return An invisible list with four elements:
#' \itemize{
#'   \item \code{acc}     — named numeric vectors of accuracy per dataset
#'   \item \code{recall}  — named numeric vectors of recall  per dataset
#'   \item \code{fs}      — named numeric vectors of F1      per dataset
#'   \item \code{summary} — lists of data.frames of predictions per dataset
#' }
#' If \code{collector} is provided, results are also written there. Use
#' \code{fh_as_lists(collector)} to convert it into plain lists.
#'
#' @examples
#' \dontrun{
#' # assume 'collector <- new.env()'
#' fh_gbm_lasso_best(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   fold         = 10,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   auto_th_method = "auto",
#'   auto_imbalance_thresh = 0.35,
#'   auto_pr_vs_roc_gate   = 0.5,
#'   collector    = collector
#' )
#' }
#' @export
fh_gbm_lasso_best <- function(
    train_exp,
    train_labels,
    test_exp,
    labels_list,
    lassogene,
    fold = 10,
    cutoff = c(0.25, 0.5, 0.75),
    cores = NULL,
    auto_th_method = "auto",
    auto_imbalance_thresh = 0.35,
    auto_pr_vs_roc_gate   = 0.5,
    collector = collector
) {
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  if (is.null(cores)) {
    cores <- tryCatch(max(1L, parallel::detectCores() - 1L), error = function(e) 1L)
  }
  
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)
  
  train_expd <- as.data.frame(train_exp)[, lassogene, drop = FALSE]
  train_expd$labels <- factor(as_binary01(train_labels), levels = c(0, 1))
  
  gbmGrid <- expand.grid(
    interaction.depth = c(3, 4, 5),
    n.trees            = c(seq(50, 300, 50), seq(400, 1000, 200)),
    shrinkage          = c(0.005, 0.01, 0.02, 0.05),
    n.minobsinnode     = c(5, 10, 15)
  )
  
  fitControl <- caret::trainControl(method = "cv", number = fold)
  GBMberT <- caret::train(
    labels ~ . - labels,
    data         = train_expd,
    distribution = "bernoulli",
    method       = "gbm",
    trControl    = fitControl,
    tuneGrid     = gbmGrid,
    bag.fraction = 0.8,
    verbose      = FALSE
  )
  
  besttune <- GBMberT$bestTune
  train_expd_final <- as.data.frame(train_exp)[, lassogene, drop = FALSE]
  train_expd_final$labels <- as.numeric(as_binary01(train_labels))
  model <- gbm::gbm(
    labels ~ . - labels,
    data               = train_expd_final,
    distribution       = "bernoulli",
    bag.fraction       = 0.8,
    n.minobsinnode     = besttune$n.minobsinnode,
    shrinkage          = besttune$shrinkage,
    n.trees            = besttune$n.trees,
    interaction.depth  = besttune$interaction.depth,
    verbose            = FALSE
  )
  n_trees <- besttune$n.trees  # <<< 显式指定预测用树数，避免 “Using xxx trees...” 提示
  
  # ---- fixed cutoffs + auto cutoff on B (merge) ----
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA_real_
  method_used <- auto_th_method
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    prob_v <- as.numeric(gbm::predict.gbm(model, newdata = val_x, type = "response", n.trees = n_trees))
    if (length(unique(prob_v[is.finite(prob_v)])) >= 2) {
      y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(
          probs = prob_v, y_true = y_val,
          imbalance_thresh = auto_imbalance_thresh,
          pr_vs_roc_gate   = auto_pr_vs_roc_gate
        )
      }
      th_auto <- choose_threshold(prob_v, y_val, method = method_used)
      
      if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
        message(sprintf("[lasso+GBM CV best] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
      } else {
        message("[lasso+GBM CV best] Auto threshold skipped; keep fixed thresholds.")
      }
    } else {
      message("[lasso+GBM CV best] B has insufficient probability variation; keep fixed thresholds.")
    }
  } else {
    message("[lasso+GBM CV best] Validation overlap < 2; keep fixed thresholds.")
  }
  
  # ---- evaluate per cutoff ----
  for (cutoffi in cutoff_now) {
    prob_tr <- as.numeric(gbm::predict.gbm(
      model, newdata = as.data.frame(train_exp)[, lassogene, drop = FALSE],
      type = "response", n.trees = n_trees
    ))
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cutoffi, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
    
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])[, lassogene, drop = FALSE]
      label_df <- labels_list[[i + 1]]
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) {
        return(data.frame(predict_p = numeric(0),
                          predict_result = factor(character()),
                          real_label = factor(character())))
      }
      expdata <- data_i[comd, , drop = FALSE]
      p  <- as.numeric(gbm::predict.gbm(model, newdata = expdata, type = "response", n.trees = n_trees))
      pr <- factor(ifelse(p > cutoffi, "positive", "negative"))
      out <- data.frame(predict_p = p, predict_result = pr)
      out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1, "positive", "negative"))
      out
    })
    
    all_result <- c(list(Train = train_result), all_result)
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 &&
      (abs(cutoffi - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("GBM-CV:", fold, " fold + Lasso-CV:", fold,
             " fold (cutoff:auto(", formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
    } else {
      paste0("GBM-CV:", fold, " fold + Lasso-CV:", fold,
             " fold (cutoff:", formatC(cutoffi, format = "f", digits = 4), ")")
    }
    
    if (is.environment(collector)) {
      .fh_collect(
        collector, key,
        acc     = result_acc,
        recall  = result_recall,
        fs      = result_FS,
        summary = all_result
      )
    }
    
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  invisible(list(
    acc     = acc_local,
    recall  = recall_local,
    fs      = fs_local,
    summary = summary_local
  ))
}

