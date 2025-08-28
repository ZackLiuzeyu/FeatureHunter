#' Naive Bayes (with Auto Threshold on B) and Collector Output
#'
#' Trains a Naive Bayes classifier (Laplace smoothing), selects a threshold on
#' validation set B (optionally auto-deciding the strategy), merges it with a
#' fixed cutoff grid, and evaluates A/B/C/D… datasets to compute ACC/RECALL/F1.
#' Results can be written into a provided `collector` environment and are also
#' accumulated into the returned lists.
#'
#' @param train_exp    Training expression matrix (rows = samples, columns = features).
#' @param train_labels Training labels (0/1 or a binary factor/coercible to 0/1).
#' @param test_exp     A list of external sets: the 1st element is validation set B,
#'                     the rest are C/D/….
#' @param labels_list  A list of label data.frames aligned with datasets
#'                     (A = \code{[[1]]}, B = \code{[[2]]}, …; column 1 = ID, column 2 = label).
#' @param cutoff       Fixed cutoff grid; will be merged with the auto-selected
#'                     threshold on B. Default: \code{c(0.25, 0.5, 0.75)}.
#' @param auto_th_method Method for auto thresholding (passed to \code{choose_threshold}):
#'                     \code{"youden"} / \code{"f1"} / \code{"auto"} (default \code{"auto"}).
#'                     With \code{"auto"}, \code{decide_threshold_method()} is used to pick
#'                     the strategy based on class imbalance and PR vs ROC performance.
#' @param auto_imbalance_thresh Imbalance threshold for \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate   PR-vs-ROC advantage ratio gate for
#'                              \code{decide_threshold_method()} (default 0.5).
#' @param collector    (Optional) environment. If supplied, each cutoff’s results
#'                     are written into it via \code{.fh_collect()}.
#'
#' @return \code{invisible(list(acc=..., recall=..., fs=..., summary=..., collector=...))}
#'         Local copies of metrics/results (and the \code{collector} if provided).
#'
#' @examples
#' \dontrun{
#' # suppose you prepared train/test objects and a collector env()
#' collector <- new.env(parent = emptyenv())
#' fh_naive_bayes(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B/C/D
#'   labels_list  = labels_list,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   collector    = collector
#' )
#' }
#' @importFrom e1071 naiveBayes
#' @export
fh_naive_bayes <- function(
    train_exp,
    train_labels,
    test_exp,
    labels_list,
    cutoff = c(0.25, 0.5, 0.75),
    auto_th_method = "auto",
    auto_imbalance_thresh = 0.35,
    auto_pr_vs_roc_gate   = 0.5,
    collector = collector
){
  # local accumulators (also returned invisibly)
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # --- fit NB with Laplace smoothing ---
  y_fac <- factor(as_binary01(train_labels), levels = c(0,1), labels = c("0","1"))
  model <- e1071::naiveBayes(
    x = as.data.frame(train_exp),
    y = y_fac,
    laplace = 1
  )
  
  # --- auto threshold on B (merge with fixed) ---
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA_real_
  method_used <- auto_th_method
  
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    y_val  <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    prob_val_df <- tryCatch(predict(model, newdata = val_x, type = "raw"),
                            error = function(e) NULL)
    if (!is.null(prob_val_df)) {
      # helper from your utils: choose the positive column
      pos_col  <- .get_pos_col(prob_val_df)
      prob_val <- as.numeric(prob_val_df[, pos_col])
      
      if (!.is_degenerate_probs(prob_val)) {
        # auto method decision if requested
        if (identical(auto_th_method, "auto")) {
          method_used <- decide_threshold_method(
            probs = prob_val, y_true = y_val,
            imbalance_thresh = auto_imbalance_thresh,
            pr_vs_roc_gate   = auto_pr_vs_roc_gate
          )
        }
        th_auto <- choose_threshold(prob_val, y_val, method = method_used)
        
        if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
          cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
          message(sprintf("[NaiveBayes] Auto threshold on B = %.4f (method=%s)",
                          th_auto, method_used))
        } else {
          message("[NaiveBayes] Auto threshold skipped; keep original cutoff grid.")
        }
      } else {
        message("[NaiveBayes] Validation probabilities degenerate; keep original cutoff grid.")
      }
    } else {
      message("[NaiveBayes] Could not get raw probabilities on B; keep original cutoff grid.")
    }
  } else {
    message("[NaiveBayes] Validation overlap < 2; keep original cutoff grid.")
  }
  
  # --- evaluate for each cutoff (include B; prepend Train) ---
  for (cuti in cutoff_now) {
    # Train (with your fallback helper)
    tr_pred <- .predict_nb_with_fallback(model, as.data.frame(train_exp), cutoff = cuti)
    train_result <- data.frame(
      predict_p = tr_pred$predict_p,
      predict_result = tr_pred$predict_result
    )
    train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
    
    # External sets: B/C/D/…
    all_result <- lapply(seq_along(test_exp), function(i){
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
      prd <- .predict_nb_with_fallback(model, expdata, cutoff = cuti)
      out <- data.frame(predict_p = prd$predict_p, predict_result = prd$predict_result)
      out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1,
                                      "positive", "negative"))
      out
    })
    
    # put Train first and name sets
    all_result <- c(list(Train = train_result), all_result)
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    # metrics
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x)
      if (nrow(x) == 0) NA_real_ else mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      if (nrow(x) == 0) return(NA_real_)
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    })
    
    # key with auto tag if applicable
    is_auto <- is.finite(th_auto) && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("NaiveBayes (cutoff:auto(",
             formatC(th_auto, format="f", digits=4), ",", method_used, "))")
    } else {
      paste0("NaiveBayes (cutoff:", formatC(cuti, format="f", digits=4), ")")
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
  
  if (!is.environment(collector)) collector <- new.env(parent = emptyenv())
  invisible(list(
    acc      = acc_local,
    recall   = recall_local,
    fs       = fs_local,
    summary  = summary_local,
    collector = collector
  ))
}



#' Naive Bayes on LASSO Features with Auto Thresholding on B
#'
#' Fit a Naive Bayes model using only `lassogene` features, optionally pick an
#' additional decision threshold automatically on validation set B (merged with
#' fixed thresholds), then evaluate on Train (A) and all external sets (B/C/D…).
#'
#' @param train_exp data.frame or matrix. Training expression matrix (samples in rows).
#' @param train_labels vector. Training labels (binary 0/1 or equivalent).
#' @param test_exp list. External expression matrices in a list; `test_exp[[1]]` is B.
#' @param labels_list list. Label tables aligned with train/test; `labels_list[[1]]` is A, `[[2]]` is B, etc.
#' @param lassogene character vector. Feature names kept by LASSO.
#' @param fold integer. Fold count for display tags in metric keys (no re-splitting done here).
#' @param cutoff numeric. Fixed cutoff thresholds; default `c(0.25, 0.5, 0.75)`. Will be merged with the auto-selected one on B.
#' @param auto_th_method Auto threshold method passed to \code{choose_threshold}:
#'   \code{"youden"} / \code{"f1"} / \code{"auto"} (default \code{"auto"}).
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to pick the strategy
#'   based on class imbalance and PR vs ROC performance; otherwise the specified method is used.
#' @param auto_imbalance_thresh numeric. Imbalance threshold passed to \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate numeric. PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector optional environment. If supplied, each cutoff’s results are written into it
#'   via \code{.fh_collect(collector, key, acc=..., recall=..., fs=..., summary=...)}.
#'
#' @return An invisible list with four elements (\code{acc}, \code{recall}, \code{fs}, \code{summary})
#'   that mirror what is also written into \code{collector} when provided.
#'
#' @examples
#' \dontrun{
#' # suppose you already prepared: train_exp/train_labels/exp_list/labels_list/lassogene/collector
#' out <- fh_naive_bayes_lasso(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   fold         = 10,
#'   collector    = collector
#' )
#' names(out$acc); names(out$summary)
#' }
#' @export
fh_naive_bayes_lasso <- function(train_exp,
                                 train_labels,
                                 test_exp,
                                 labels_list,
                                 lassogene,
                                 fold = 10,
                                 cutoff = c(0.25, 0.5, 0.75),
                                 auto_th_method = "auto",
                                 auto_imbalance_thresh = 0.35,
                                 auto_pr_vs_roc_gate   = 0.5,
                                 collector = collector) {
  # ---- Train NB on LASSO features ----
  y_fac <- factor(train_labels, levels = c(0, 1), labels = c("0", "1"))
  model <- e1071::naiveBayes(
    x = as.data.frame(train_exp[, lassogene, drop = FALSE]),
    y = y_fac
  )
  
  # ---- Auto threshold on B + merge with fixed ----
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA_real_
  method_used <- auto_th_method
  
  if (length(comd_val) >= 2) {
    val_x <- val_exp[comd_val, lassogene, drop = FALSE]
    prob_val_df <- tryCatch(predict(model, newdata = val_x, type = "raw"),
                            error = function(e) NULL)
    if (!is.null(prob_val_df)) {
      pos_col  <- .get_pos_col(prob_val_df)
      prob_val <- as.numeric(prob_val_df[, pos_col])
      y_val    <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (!.is_degenerate_probs(prob_val)) {
        if (identical(auto_th_method, "auto")) {
          method_used <- decide_threshold_method(
            probs = prob_val, y_true = y_val,
            imbalance_thresh = auto_imbalance_thresh,
            pr_vs_roc_gate   = auto_pr_vs_roc_gate
          )
        }
        th_auto <- choose_threshold(prob_val, y_val, method = method_used)
        
        if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
          cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
          message(sprintf("[lasso+NaiveBayes] Auto threshold on B = %.4f (method=%s)",
                          th_auto, method_used))
        } else {
          message("[lasso+NaiveBayes] Auto threshold skipped; keep fixed thresholds.")
        }
      } else {
        message("[lasso+NaiveBayes] Validation probs degenerate; keep fixed thresholds.")
      }
    } else {
      message("[lasso+NaiveBayes] Could not get raw probs on B; keep fixed thresholds.")
    }
  } else {
    message("[lasso+NaiveBayes] Validation overlap < 2; keep fixed thresholds.")
  }
  
  # ---- Local accumulators for return ----
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # ---- Evaluate for each cutoff ----
  for (cuti in cutoff_now) {
    # Train predictions with robust helper
    prd_tr <- .predict_nb_with_fallback(
      model,
      as.data.frame(train_exp[, lassogene, drop = FALSE]),
      cutoff = cuti
    )
    train_result <- data.frame(
      predict_p      = prd_tr$predict_p,
      predict_result = prd_tr$predict_result
    )
    train_result$real_label <- factor(ifelse(train_labels == 1, "positive", "negative"))
    
    # External sets (include B)
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])
      label_df <- labels_list[[i + 1]]
      ids_i    <- as.character(label_df[[1]])
      comd     <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) {
        return(data.frame(predict_p = numeric(0),
                          predict_result = factor(character()),
                          real_label = factor(character())))
      }
      expdata <- data_i[comd, lassogene, drop = FALSE]
      
      prd <- .predict_nb_with_fallback(model, expdata, cutoff = cuti)
      
      out <- data.frame(
        predict_p      = prd$predict_p,
        predict_result = prd$predict_result
      )
      out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1,
                                      "positive", "negative"))
      out
    })
    
    # Put Train in front
    all_result <- c(list(Train = train_result), all_result)
    
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    # Metrics
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x)
      mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" &
                  as.character(x$real_label)    == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" &
                  as.character(x$real_label)    == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    # Key labeling (auto/fixed)
    is_auto <- is.finite(th_auto) && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("lasso+NaiveBayes-CV:", fold, " fold (cutoff:auto(",
             formatC(th_auto, format="f", digits=4), ",", method_used, "))")
    } else {
      paste0("lasso+NaiveBayes-CV:", fold, " fold (cutoff:",
             formatC(cuti, format="f", digits=4), ")")
    }
    
    # Write into collector if provided
    if (is.environment(collector)) {
      .fh_collect(
        collector, key,
        acc     = result_acc,
        recall  = result_recall,
        fs      = result_FS,
        summary = all_result
      )
    }
    
    # Also keep local copies for return
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  if (!is.environment(collector)) collector <- new.env(parent = emptyenv())
  invisible(list(
    acc      = acc_local,
    recall   = recall_local,
    fs       = fs_local,
    summary  = summary_local,
    collector = collector
  ))
}

