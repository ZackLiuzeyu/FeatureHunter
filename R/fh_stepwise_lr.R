#' Stepwise Logistic Regression with AIC
#'
#' Forward stepwise logistic regression using AIC for feature selection.
#' After training, a threshold is chosen on validation set B (either a fixed
#' method: "youden"/"f1", or auto via \code{decide_threshold_method()}),
#' merged with fixed cutoffs, and evaluated on A/B/C/… datasets to compute
#' F1/Accuracy/Recall. If a \code{collector} environment is supplied, results
#' are written into it via \code{.fh_collect()}; results are also accumulated
#' locally and returned invisibly.
#'
#' @param train_exp data.frame. Training expression matrix (rows = samples, columns = features).
#' @param train_labels vector. Binary labels for the training set (0/1 or convertible).
#' @param test_exp list. External expression matrices; the first is validation set B, then C/D/….
#' @param labels_list list. Label data.frames corresponding to datasets (A = 1, B = 2, …; col1 = ID, col2 = label).
#' @param cutoff numeric. Fixed cutoff thresholds (default \code{c(0.25, 0.5, 0.75)}); merged with the auto threshold from set B.
#' @param auto_th_method character. \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to pick a strategy
#'   based on class imbalance and PR-vs-ROC performance.
#' @param auto_imbalance_thresh numeric. Imbalance threshold passed to \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate numeric. PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector environment (optional). If supplied, each cutoff’s results are written
#'   into it via \code{.fh_collect()}.
#'
#' @return
#' Invisibly returns a list with four elements:
#' \itemize{
#'   \item \code{acc}: named vectors of accuracy across datasets per key
#'   \item \code{recall}: named vectors of recall across datasets per key
#'   \item \code{fs}: named vectors of F1 across datasets per key
#'   \item \code{summary}: lists of per-dataset prediction data.frames per key
#' }
#' If \code{collector} is provided, results are also written into it using the same keys.
#'
#' @examples
#' \dontrun{
#' # initialize a collector environment
#' collector <- new.env(parent = emptyenv())
#'
#' fh_stepwise_logistic(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]), # B/C/D
#'   labels_list  = labels_list,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   collector    = collector
#' )
#' }
#'
#' @importFrom stats glm step predict
#' @importFrom MASS stepAIC
#' @export
fh_stepwise_lr <- function(
    train_exp,
    train_labels,
    test_exp,
    labels_list,
    cutoff = c(0.25, 0.5, 0.75),
    auto_th_method = "auto",
    auto_imbalance_thresh = 0.35,
    auto_pr_vs_roc_gate   = 0.5,
    collector = collector
) {
  # ---- helper: silence only specific glm warnings ----
  .muffle_glm_warnings <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("glm.fit: algorithm did not converge", msg) ||
            grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", msg)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  
  # local accumulators (also returned invisibly)
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # 1) Fit forward stepwise logistic regression (AIC)
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- as_binary01(train_labels)
  
  fullModel <- .muffle_glm_warnings(
    stats::glm(labels ~ . - labels, family = "binomial", data = train_expd)
  )
  nullModel <- .muffle_glm_warnings(
    stats::glm(labels ~ 1,        family = "binomial", data = train_expd)
  )
  
  model <- .muffle_glm_warnings(
    MASS::stepAIC(
      nullModel,
      direction = "forward",
      scope     = list(upper = fullModel, lower = nullModel),
      trace     = 0
    )
  )
  
  # 2) Auto threshold on validation set B, merge with fixed grid
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
    
    prob_val <- as.numeric(stats::predict(model, newdata = val_x, type = "response"))
    
    if (length(unique(prob_val[is.finite(prob_val)])) >= 2) {
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
        message(sprintf("[Stepwise LR] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
      } else {
        message("[Stepwise LR] Auto threshold skipped; keep original cutoff grid.")
      }
    } else {
      message("[Stepwise LR] B has insufficient probability variation; keep original cutoff grid.")
    }
  } else {
    message("[Stepwise LR] Validation overlap < 2; keep original cutoff grid.")
  }
  
  # 3) Evaluate for each cutoff over A/B/C/…
  for (cuti in cutoff_now) {
    # Train predictions
    prob_tr <- as.numeric(stats::predict(model, type = "response"))
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
    
    # External datasets (B/C/D/…)
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
      
      p  <- as.numeric(stats::predict(model, newdata = expdata, type = "response"))
      pr <- factor(ifelse(p > cuti, "positive", "negative"))
      out <- data.frame(predict_p = p, predict_result = pr)
      out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1, "positive", "negative"))
      out
    })
    
    # prepend Train (keep consistent naming order)
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
      mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    # key: mark auto cutoff if used
    is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 &&
      (abs(cuti - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("StepWise-AIC+LR (cutoff:auto(",
             formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
    } else {
      paste0("StepWise-AIC+LR (cutoff:", formatC(cuti, format = "f", digits = 4), ")")
    }
    
    # write into collector if provided
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

