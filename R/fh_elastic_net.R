#' Elastic Net interface (cv.glmnet)
#'
#' @description
#' Fits Elastic Net models over a set of `alpha` values (default 0.1.9) via
#' `glmnet::cv.glmnet`, extracts `lambda.min`, selects an automatic threshold on
#' validation set B (F1/Youden or auto), merges it with fixed cutoffs, and evaluates
#' F1/ACC/RECALL on datasets A/B/C/D.
#'
#' @param train_exp    Training expression matrix (rows = samples, cols = features).
#' @param train_labels Training labels (0/1 or convertible to 0/1).
#' @param test_exp     A list: first element = validation set B, others = C/D/.
#' @param labels_list  List of label data.frames (A=1, B=2, ; first two columns = ID and label).
#' @param fold         Number of folds for CV (passed to `cv.glmnet`; also recorded in result keys). Default: 10.
#' @param alpha_all    Sequence of alpha values to search. Default: `seq(0.1, 0.9, 0.1)`.
#' @param cutoff       Fixed cutoff set merged with the auto threshold from B. Default: `c(0.25, 0.5, 0.75)`.
#' @param auto_th_method Auto threshold method: `"youden"` / `"f1"` / `"auto"`. With `"auto"`,
#'   `decide_threshold_method()` is used to pick between F1 and Youden based on class imbalance
#'   and PR-vs-ROC performance; otherwise the specified method is used. Default: `"auto"`.
#' @param auto_imbalance_thresh Imbalance threshold passed to `decide_threshold_method()` (default 0.35).
#' @param auto_pr_vs_roc_gate   PR-vs-ROC advantage ratio gate for `decide_threshold_method()` (default 0.5).
#' @param collector    (optional) environment used as a result collector. If supplied, each cutoff
#'   results are written into it via `.fh_collect()`.
#'
#' @return
#' A `collector` environment containing four fields (`all_result_summary`, `all_result_acc`,
#' `all_result_recall`, `all_result_FS`) **and** an invisible list with the same four
#' sub-lists (`acc`, `recall`, `fs`, `summary`) accumulated for this call. Use
#' `fh_as_lists(collector)` to convert the environment into plain R lists.
#'
#' @examples
#' \dontrun{
#' # Prepare a collector
#' collector <- new.env(parent = emptyenv())
#'
#' fh_elastic_net(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   fold         = 10,
#'   alpha_all    = c(0.1, 0.5, 0.9),
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   auto_th_method = "auto",
#'   auto_imbalance_thresh = 0.35,
#'   auto_pr_vs_roc_gate   = 0.5,
#'   collector    = collector
#' )
#'
#' # Convert to plain lists if needed:
#' lists <- fh_as_lists(collector)
#' }
#'
#' @importFrom glmnet cv.glmnet
#' @export
fh_elastic_net <- function(train_exp, train_labels, test_exp, labels_list,
                           fold = 10,
                           alpha_all = seq(0.1, 0.9, 0.1),
                           cutoff = c(0.25, 0.5, 0.75),
                           auto_th_method = "auto",
                           auto_imbalance_thresh = 0.35,
                           auto_pr_vs_roc_gate   = 0.5,
                           collector = NULL) {
  # local accumulators to be returned invisibly
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  for (alphai in alpha_all) {
    ## ---- train cv.glmnet ----
    train_expd <- as.matrix(train_exp)
    y_glmnet   <- as.numeric(as_binary01(train_labels))
    cvfit <- glmnet::cv.glmnet(
      x = train_expd, y = y_glmnet,
      family = "binomial",
      nlambda = 100, alpha = alphai,
      nfolds = fold
    )
    lambda_min <- cvfit$lambda.min
    # CV error at lambda.min
    cv_error <- {
      idx <- which(cvfit$lambda == lambda_min)
      if (length(idx) == 1L) cvfit$cvm[idx] else NA_real_
    }
    
    ## ---- baseline cutoffs ----
    cutoff_now <- cutoff
    
    ## ---- auto threshold on B, merge with fixed ----
    val_exp <- as.data.frame(test_exp[[1]])
    val_lab <- labels_list[[2]]
    ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    
    # default
    th_auto <- NA_real_
    method_used <- auto_th_method
    
    if (length(comd_val) >= 2) {
      val_x    <- as.matrix(val_exp[comd_val, , drop = FALSE])
      prob_val <- as.numeric(predict(cvfit, newx = val_x, s = lambda_min, type = "response"))
      y_val    <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(
          probs = prob_val, y_true = y_val,
          imbalance_thresh = auto_imbalance_thresh,
          pr_vs_roc_gate   = auto_pr_vs_roc_gate
        )
      }
      
      th_auto <- choose_threshold(prob_val, y_val, method = method_used)
      
      if (is.finite(th_auto) && length(th_auto) == 1L && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
        message(sprintf(
          "[ENR-CV alpha=%.2f] Auto threshold on B = %.4f (method=%s, CVerror=%.4f)",
          alphai, th_auto, method_used, cv_error
        ))
      } else {
        message(sprintf(
          "[ENR-CV alpha=%.2f] Auto threshold skipped; keep fixed thresholds. (CVerror=%.4f)",
          alphai, cv_error
        ))
      }
    } else {
      message(sprintf(
        "[ENR-CV alpha=%.2f] Validation overlap < 2; keep fixed thresholds. (CVerror=%.4f)",
        alphai, cv_error
      ))
    }
    
    ## ---- evaluate over cutoffs ----
    for (cuti in cutoff_now) {
      # Train (A)
      prob_tr <- as.numeric(predict(cvfit, newx = train_expd, s = lambda_min, type = "response"))
      train_result <- data.frame(
        predict_p = prob_tr,
        predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
      )
      train_result$real_label <- factor(ifelse(train_labels == 1, "positive", "negative"))
      
      # External sets: B/C/D/…
      all_result <- lapply(seq_along(test_exp), function(i) {
        data_i   <- as.data.frame(test_exp[[i]])
        label_df <- labels_list[[i + 1]]  # A=1 → B=2 → C=3 …
        
        ids_i <- as.character(label_df[[1]])
        comd  <- intersect(rownames(data_i), ids_i)
        if (length(comd) == 0) {
          return(data.frame(
            predict_p = numeric(0),
            predict_result = factor(character()),
            real_label = factor(character())
          ))
        }
        
        expdata    <- as.matrix(data_i[comd, , drop = FALSE])
        labelsdata <- label_df[match(comd, ids_i), 2]
        
        p  <- as.numeric(predict(cvfit, newx = expdata, s = lambda_min, type = "response"))
        pr <- factor(ifelse(p > cuti, "positive", "negative"))
        
        out <- data.frame(predict_p = p, predict_result = pr)
        out$real_label <- factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative"))
        out
      })
      
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
      
      # key with alpha + fold + cutoff (auto/fixed)
      is_auto <- is.finite(th_auto) && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
      key <- if (is_auto) {
        paste0("ENR-CV:", fold, " fold (alpha:", formatC(alphai, format = "f", digits = 2),
               ") (cutoff:auto(", formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
      } else {
        paste0("ENR-CV:", fold, " fold (alpha:", formatC(alphai, format = "f", digits = 2),
               ") (cutoff:", formatC(cuti, format = "f", digits = 4), ")")
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
      
      # also accumulate to local return value
      acc_local[[key]]     <- result_acc
      recall_local[[key]]  <- result_recall
      fs_local[[key]]      <- result_FS
      summary_local[[key]] <- all_result
    }
  }
  
  # return collector (and invisible lists)
  if (!is.environment(collector)) collector <- new.env(parent = emptyenv())
  invisible(list(
    acc       = acc_local,
    recall    = recall_local,
    fs        = fs_local,
    summary   = summary_local,
    collector = collector
  ))
}

