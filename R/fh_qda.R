#' @title Quadratic Discriminant Analysis (QDA) Interface
#'
#' @description
#' Trains a \code{MASS::qda} model on the training set;  
#' automatically selects threshold on validation set B (via \code{"f1"}, \code{"youden"}, or \code{"auto"}),  
#' merges it with fixed cutoffs, and evaluates Accuracy / Recall / F1 across datasets A/B/C/D…  
#' Results are written into the \code{collector} environment.
#'
#' @param train_exp Training expression matrix (samples x features).
#' @param train_labels Binary vector or factor of training labels (0/1).
#' @param test_exp A list: first element is validation set B; remaining elements are test sets C/D/…
#' @param labels_list Label list used to align external datasets (A=1, B=2, C=3, …; first two columns = ID and label).
#' @param cutoff Numeric vector of thresholds (default = c(0.25, 0.5, 0.75)); merged with auto threshold from B.
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to adaptively select
#'   based on class imbalance and PR vs ROC performance; otherwise the specified method is used.
#' @param auto_imbalance_thresh Imbalance threshold passed to \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector (optional) environment collector; if supplied, each cutoff’s results are written into it.
#'
#' @return Returns a \code{collector} environment containing four fields:  
#' \itemize{
#'   \item \code{all_result_summary}: predicted result data.frames for each dataset  
#'   \item \code{all_result_acc}: accuracy values for each dataset  
#'   \item \code{all_result_recall}: recall values for each dataset  
#'   \item \code{all_result_FS}: F1 scores for each dataset  
#' }  
#' Use \code{fh_as_lists(collector)} to convert into plain R lists.
#'
#' @examples
#' \dontrun{
#' collector <- .fh_new_collector()
#' fh_qda(
#'   train_exp       = train_exp,
#'   train_labels    = train_labels,
#'   test_exp        = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list     = labels_list,
#'   cutoff          = c(0.25, 0.5, 0.75),
#'   auto_th_method  = "auto",
#'   collector       = collector
#' )
#' results <- fh_as_lists(collector)
#' str(results$all_result_FS)
#' }
#'
#' @importFrom MASS qda
#' @importFrom dplyr %>% mutate across
#' @export
fh_qda <- function(train_exp, train_labels, test_exp, labels_list,
                   cutoff = c(0.25, 0.5, 0.75),
                   auto_th_method = "auto",
                   auto_imbalance_thresh = 0.35,
                   auto_pr_vs_roc_gate   = 0.5,
                   collector = collector) {
  train_expd <- as.data.frame(train_exp)
  train_expd <- train_expd %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.na(.), stats::median(., na.rm=TRUE), .)))
  train_expd$labels <- train_labels
  
  model <- qda(labels ~ . - labels, data = train_expd)
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    post_v <- predict(model, newdata = val_x)$posterior
    pos_col <- intersect(colnames(post_v), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col) == 0) pos_col <- tail(colnames(post_v), 1)
    prob_val <- as.numeric(post_v[, pos_col[1]])
    y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    method_used <- auto_th_method
    if (identical(auto_th_method, "auto")) {
      method_used <- decide_threshold_method(
        probs = prob_val, y_true = y_val,
        imbalance_thresh = auto_imbalance_thresh,
        pr_vs_roc_gate   = auto_pr_vs_roc_gate
      )
    }
    th_auto <- choose_threshold(prob_val, y_val, method = method_used)
    if (is.finite(th_auto) && length(th_auto) == 1L  && th_auto > 0.01 && th_auto < 0.99) {
      cutoff <- sort(unique(c(cutoff, th_auto)))
      message(sprintf("[QDA] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[QDA] Auto threshold skipped; keep original cutoff grid.")
    }
  }else {
    message("[QDA] Validation overlap < 2; keep original cutoff grid.")
  }
  
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  
  for (cuti in cutoff) {
    post_tr <- predict(model, newdata = as.data.frame(train_exp))$posterior
    pos_col_tr <- intersect(colnames(post_tr), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col_tr) == 0) pos_col_tr <- tail(colnames(post_tr), 1)
    prob_tr <- as.numeric(post_tr[, pos_col_tr[1]])
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative")),
      real_label = factor(ifelse(train_labels == 1, "positive", "negative"))
    )
    
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])
      label_df <- labels_list[[i + 1]]
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) return(data.frame(predict_p=numeric(0),
                                               predict_result=factor(character()),
                                               real_label=factor(character())))
      expdata    <- data_i[comd, , drop = FALSE]
      labelsdata <- label_df[match(comd, ids_i), 2]
      post_i <- predict(model, newdata = expdata)$posterior
      pos_col_i <- intersect(colnames(post_i), c("1","positive","POS","Yes","TRUE"))
      if (length(pos_col_i) == 0) pos_col_i <- tail(colnames(post_i), 1)
      p  <- as.numeric(post_i[, pos_col_i[1]])
      pr <- factor(ifelse(p > cuti, "positive", "negative"))
      data.frame(predict_p = p, predict_result = pr,
                 real_label = factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative")))
    })
    all_result <- c(list(Train=train_result), all_result)
    nice_names <- c("DatasetA(Train)",
                    paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
                           "(", c("Val", rep("Test", length(test_exp)-1)), ")"))
    names(all_result) <- nice_names[seq_along(all_result)]
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result)=="positive" & as.character(x$real_label)=="positive")
      fn <- sum(as.character(x$predict_result)=="negative" & as.character(x$real_label)=="positive")
      if ((tp+fn)==0) 0 else tp/(tp+fn)
    })
    
    is_auto <- is.finite(th_auto) && (abs(cuti-th_auto)<.Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("QDA (cutoff:auto(", formatC(th_auto, format="f", digits=4), ",", method_used, "))")
    } else paste0("QDA (cutoff:", formatC(cuti, format="f", digits=4), ")")
    
    if (is.environment(collector)) {
      .fh_collect(collector, key, acc=result_acc, recall=result_recall, fs=result_FS, summary=all_result)
    }
    acc_local[[key]] <- result_acc
    recall_local[[key]] <- result_recall
    fs_local[[key]] <- result_FS
    summary_local[[key]] <- all_result
  }
  invisible(list(
    acc     = acc_local,
    recall  = recall_local,
    fs      = fs_local,
    summary = summary_local
  ))
}


#' @title LASSO + Quadratic Discriminant Analysis (QDA) Interface
#'
#' @description
#' Applies LASSO-selected features (\code{lassogene}) to filter the expression matrix,
#' then trains a \code{MASS::qda} model.  
#' The classification threshold is selected on validation set B (via \code{"f1"}, \code{"youden"}, or \code{"auto"}),
#' merged with fixed cutoffs, and used to evaluate F1 / Accuracy / Recall across datasets A/B/C/D…  
#' Results are written into the \code{collector} environment.
#'
#' @param train_exp Training expression matrix (samples x features).
#' @param train_labels Binary vector or factor of training labels (0/1).
#' @param test_exp A list: first element is validation set B; remaining elements are test sets C/D/…
#' @param labels_list A list of label data.frames (A=1, B=2, C=3, …) used for dataset alignment.
#' @param lassogene Character vector of LASSO-selected feature names.
#' @param fold Integer; only used for labeling in result keys (not used in computation).
#' @param cutoff Numeric vector of thresholds (default = c(0.25, 0.5, 0.75)); merged with auto threshold from B.
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is called to adaptively select
#'   based on class imbalance and PR vs ROC performance; otherwise the specified method is used.
#' @param auto_imbalance_thresh Imbalance threshold for \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
#' @param collector (optional) environment collector; if supplied, results for each cutoff are written into it.
#'
#' @return Returns a \code{collector} environment containing four fields:  
#' \itemize{
#'   \item \code{all_result_summary}: predicted result data.frames for each dataset  
#'   \item \code{all_result_acc}: accuracy values for each dataset  
#'   \item \code{all_result_recall}: recall values for each dataset  
#'   \item \code{all_result_FS}: F1 scores for each dataset  
#' }  
#' Use \code{fh_as_lists(collector)} to convert into plain R lists.
#'
#' @examples
#' \dontrun{
#' collector <- .fh_new_collector()
#' fh_lasso_qda(
#'   train_exp       = train_exp,
#'   train_labels    = train_labels,
#'   test_exp        = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list     = labels_list,
#'   lassogene       = lassogene,
#'   fold            = 10,
#'   cutoff          = c(0.25, 0.5, 0.75),
#'   auto_th_method  = "auto",
#'   collector       = collector
#' )
#' results <- fh_as_lists(collector)
#' str(results$all_result_acc)
#' }
#'
#' @importFrom MASS qda
#' @importFrom dplyr %>% mutate across
#' @export
fh_lasso_qda <- function(train_exp, train_labels, test_exp, labels_list,
                         lassogene, fold,
                         cutoff = c(0.25, 0.5, 0.75),
                         auto_th_method = "auto",
                         auto_imbalance_thresh = 0.35,
                         auto_pr_vs_roc_gate   = 0.5,
                         collector = collector) {
  train_expd <- as.data.frame(train_exp)[, lassogene, drop = FALSE]
  train_expd <- train_expd %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.na(.), stats::median(., na.rm=TRUE), .)))
  train_expd$labels <- train_labels
  model <- qda(labels ~ . - labels, data = train_expd)
  
  val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    post_v <- predict(model, newdata = val_x)$posterior
    pos_col <- intersect(colnames(post_v), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col) == 0) pos_col <- tail(colnames(post_v), 1)
    prob_val <- as.numeric(post_v[, pos_col[1]])
    y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    method_used <- auto_th_method
    if (identical(auto_th_method, "auto")) {
      method_used <- decide_threshold_method(
        probs = prob_val, y_true = y_val,
        imbalance_thresh = auto_imbalance_thresh,
        pr_vs_roc_gate   = auto_pr_vs_roc_gate
      )
    }
    th_auto <- choose_threshold(prob_val, y_val, method = method_used)
    if (is.finite(th_auto) && length(th_auto) == 1L  && th_auto > 0.01 && th_auto < 0.99) {
      cutoff <- sort(unique(c(cutoff, th_auto)))
      message(sprintf("[lasso+QDA] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[lasso+QDA] Auto threshold skipped; keep original cutoff grid.")
    }
  }else {
    message("[lasso+QDA] Validation overlap < 2; keep original cutoff grid.")
  }
  
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  
  for (cuti in cutoff) {
    post_tr <- predict(model, newdata = as.data.frame(train_exp)[, lassogene, drop = FALSE])$posterior
    pos_col_tr <- intersect(colnames(post_tr), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col_tr) == 0) pos_col_tr <- tail(colnames(post_tr), 1)
    prob_tr <- as.numeric(post_tr[, pos_col_tr[1]])
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative")),
      real_label = factor(ifelse(train_labels == 1, "positive", "negative"))
    )
    
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])[, lassogene, drop = FALSE]
      label_df <- labels_list[[i + 1]]
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) return(data.frame(predict_p=numeric(0),
                                               predict_result=factor(character()),
                                               real_label=factor(character())))
      expdata    <- data_i[comd, , drop = FALSE]
      labelsdata <- label_df[match(comd, ids_i), 2]
      post_i <- predict(model, newdata = expdata)$posterior
      pos_col_i <- intersect(colnames(post_i), c("1","positive","POS","Yes","TRUE"))
      if (length(pos_col_i) == 0) pos_col_i <- tail(colnames(post_i), 1)
      p  <- as.numeric(post_i[, pos_col_i[1]])
      pr <- factor(ifelse(p > cuti, "positive", "negative"))
      data.frame(predict_p = p, predict_result = pr,
                 real_label = factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative")))
    })
    all_result <- c(list(Train=train_result), all_result)
    nice_names <- c("DatasetA(Train)",
                    paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
                           "(", c("Val", rep("Test", length(test_exp)-1)), ")"))
    names(all_result) <- nice_names[seq_along(all_result)]
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result)=="positive" & as.character(x$real_label)=="positive")
      fn <- sum(as.character(x$predict_result)=="negative" & as.character(x$real_label)=="positive")
      if ((tp+fn)==0) 0 else tp/(tp+fn)
    })
    
    is_auto <- is.finite(th_auto) && (abs(cuti-th_auto)<.Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("QDA+Lasso-CV:", fold, " fold (cutoff:auto(", formatC(th_auto, format="f", digits=4), ",", method_used, "))")
    } else paste0("QDA+Lasso-CV:", fold, " fold (cutoff:", formatC(cuti, format="f", digits=4), ")")
    
    if (is.environment(collector)) {
      .fh_collect(collector, key, acc=result_acc, recall=result_recall, fs=result_FS, summary=all_result)
    }
    acc_local[[key]] <- result_acc
    recall_local[[key]] <- result_recall
    fs_local[[key]] <- result_FS
    summary_local[[key]] <- all_result
  }
  invisible(list(
    acc     = acc_local,
    recall  = recall_local,
    fs      = fs_local,
    summary = summary_local
  ))
}


