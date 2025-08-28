#' LASSO interface (cv.glmnet, alpha = 1)
#'
#' @description
#' Fits \code{glmnet::cv.glmnet(alpha = 1)} on the training set to obtain \code{lambda.min}.
#' On validation set B, selects a decision threshold and merges it with fixed cutoffs.
#' Evaluates F1 / Accuracy / Recall on datasets A/B/C/D… and records results.
#'
#' @param train_exp    Training expression matrix (rows = samples, cols = features).
#' @param train_labels Training labels (0/1 or convertible to 0/1).
#' @param test_exp     A list: first element = validation set B, others = C/D/….
#' @param labels_list  List of label data.frames (A=1, B=2, …; first two cols = ID, label).
#' @param fold         Number of CV folds for \code{cv.glmnet} and stored in result keys.
#' @param cutoff       Fixed cutoff grid merged with the auto-selected threshold
#'                     (default \code{c(0.25, 0.5, 0.75)}).
#' @param auto_th_method Threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"} (default).
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to adaptively choose based on
#'   class imbalance and PR vs ROC performance.
#' @param auto_imbalance_thresh Imbalance threshold for \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate   PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector    Optional environment. If supplied, each cutoff’s results are written
#'                     into \code{collector} via \code{.fh_collect()}.
#'
#' @return
#' Returns a \strong{collector environment} (created if not supplied) populated with
#' \code{all_result_summary}, \code{all_result_acc}, \code{all_result_recall}, \code{all_result_FS};
#' and (invisibly) a plain list of the current-call results:
#' \code{list(acc = ..., recall = ..., fs = ..., summary = ...)} for convenience.
#' Use \code{fh_as_lists(collector)} to convert the environment to plain lists.
#'
#' @examples
#' \dontrun{
#' # Initialize a collector (recommended)
#' collector <- new.env(parent = emptyenv())
#'
#' res <- fh_lasso(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   fold         = 10,
#'   collector   = collector
#' )
#'
#' # Convert the collector to plain R lists (if needed):
#' lst <- fh_as_lists(collector)
#' str(lst$all_result_acc)
#' }
#'
#' @importFrom glmnet cv.glmnet
#' @export
fh_lasso <- function(train_exp, train_labels, test_exp, labels_list,
                     fold = 10,
                     cutoff = c(0.25, 0.5, 0.75),
                     auto_th_method = "auto",
                     auto_imbalance_thresh = 0.35,
                     auto_pr_vs_roc_gate   = 0.5,
                     collector = NULL) {
  
  # ensure collector env
  if (!is.environment(collector)) {
    collector <- new.env(parent = emptyenv())
  }
  
  # local accumulators to also return as plain lists
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # ---- train cv.glmnet (alpha=1) ----
  x_mat   <- as.matrix(train_exp)
  y_glmnet <- as.numeric(as_binary01(train_labels))
  cvfit <- glmnet::cv.glmnet(
    x = x_mat, y = y_glmnet,
    family = "binomial",
    nlambda = 100, alpha = 1,
    nfolds = fold
  )
  lambda_min <- cvfit$lambda.min
  
  # baseline cutoff grid
  cutoff_now <- cutoff
  
  ## ---- auto-threshold on B and merge ----
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  if (length(comd_val) >= 2) {
    val_x    <- as.matrix(val_exp[comd_val, , drop = FALSE])
    prob_val <- as.numeric(predict(cvfit, newx = val_x, s = lambda_min, type = "response"))
    y_val    <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    method_used <- auto_th_method
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
      message(sprintf("[Lasso-CV] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[Lasso-CV] Auto threshold skipped; keep fixed thresholds.")
    }
  } else {
    message("[Lasso-CV] Validation overlap < 2; keep fixed thresholds.")
  }
  
  ## ---- evaluate across cutoffs ----
  for (cuti in cutoff_now) {
    
    # A/train
    prob_tr <- as.numeric(predict(cvfit, newx = x_mat, s = lambda_min, type = "response"))
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
    
    # external sets (B/C/D/…)
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])
      label_df <- labels_list[[i + 1]]  # A=1 -> B=2 -> C=3…
      
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0L) {
        return(data.frame(predict_p = numeric(0),
                          predict_result = factor(character()),
                          real_label = factor(character())))
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
             "(", c("Val", rep("Test", length(test_exp) - 1L)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    # metrics
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x)
      mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" &
                  as.character(x$real_label)   == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" &
                  as.character(x$real_label)   == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    # key (distinguish auto/fixed)
    is_auto <- exists("th_auto", inherits = FALSE) && is.finite(th_auto) &&
      (abs(cuti - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("Lasso.R-CV:", fold, " fold (cutoff:auto(",
             formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
    } else {
      paste0("Lasso.R-CV:", fold, " fold (cutoff:",
             formatC(cuti, format = "f", digits = 4), ")")
    }
    
    # write into collector if provided
    .fh_collect(
      collector, key,
      acc     = result_acc,
      recall  = result_recall,
      fs      = result_FS,
      summary = all_result
    )
    
    # also keep local copies to return
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  invisible(list(
    collector = collector,
    acc       = acc_local,
    recall    = recall_local,
    fs        = fs_local,
    summary   = summary_local
  ))
}

