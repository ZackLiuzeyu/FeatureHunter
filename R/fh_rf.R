#' Random Forest Interface
#'
#' @description
#' Train a Random Forest model with candidate mtry values.
#' Steps:
#' 1) Train with 500 trees to get variable importance and pick the best number of trees via OOB error;
#' 2) Retrain with the best number of trees;
#' 3) On validation set B, determine an automatic threshold (Youden/F1/auto) and merge with fixed thresholds;
#' 4) Evaluate F1/Accuracy/Recall on datasets A/B/C/D...;
#' 5) Variable importance is stored into \code{collector$all_result_importance} if provided.
#'
#' @param train_exp Training expression matrix (samples x features).
#' @param train_labels Binary/factor labels for training set.
#' @param test_exp list of external datasets (first = validation B, then C/D/...).
#' @param labels_list list of label data.frames (A=1, B=2, ...; first col=ID, second col=label).
#' @param com_genes character vector of gene names, used to determine mtry candidates.
#' @param mnum candidate mtry values (default = 25%/50%/75% quantiles of \code{length(com_genes)}).
#' @param cutoff numeric thresholds (default = \code{c(0.25, 0.5, 0.75)}).
#' @param auto_th_method threshold method: "auto" (default), "f1", or "youden".
#' @param auto_imbalance_thresh imbalance threshold for auto mode (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage gate (default 0.5).
#' @param collector environment optional results collector.
#'
#' @return collector environment and an invisible list (acc, recall, fs, summary).
#'
#' @examples
#' \dontrun{
#' collector <- new.env()
#' fh_rf(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),
#'   labels_list  = labels_list,
#'   com_genes    = com_genes,
#'   collector    = collector
#' )
#' }
#'
#' @importFrom randomForest randomForest
#' @export
fh_rf <- function(train_exp, train_labels, test_exp, labels_list, com_genes,
                  mnum   = c(round(stats::quantile(seq_along(com_genes), probs = 0.25)),
                             round(stats::quantile(seq_along(com_genes), probs = 0.5)),
                             round(stats::quantile(seq_along(com_genes), probs = 0.75))),
                  cutoff = c(0.25, 0.5, 0.75),
                  auto_th_method = "auto",
                  auto_imbalance_thresh = 0.35,
                  auto_pr_vs_roc_gate = 0.5,
                  collector = collector) {
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  
  for (mnumberi in mnum) {
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(train_labels)
    
    model0 <- randomForest::randomForest(labels ~ . - labels, data = train_expd, ntree = 500, mtry = mnumberi)
    if (is.environment(collector)) {
      collector$all_result_importance[[paste0("RF (mtry=", mnumberi, ")")]] <-
        as.data.frame(randomForest::importance(model0))
    }
    optionTrees <- which.min(model0$err.rate[, 1])
    model <- randomForest::randomForest(labels ~ . - labels, data = train_expd, ntree = optionTrees, mtry = mnumberi)
    
    # --- auto threshold on B
    cutoff_grid <- cutoff
    th_auto <- NA_real_
    method_used <- auto_th_method
    
    val_exp <- as.data.frame(test_exp[[1]])
    val_lab <- labels_list[[2]]
    ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    
    if (length(comd_val) >= 2) {
      prob_mat <- predict(model, newdata = val_exp[comd_val, , drop = FALSE], type = "prob")
      pos_col <- intersect(colnames(prob_mat), c("1","positive","POS","Yes","TRUE"))
      if (!length(pos_col)) pos_col <- tail(colnames(prob_mat), 1)
      prob_val <- as.numeric(prob_mat[, pos_col[1]])
      y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(prob_val, y_val,
                                               imbalance_thresh = auto_imbalance_thresh,
                                               pr_vs_roc_gate   = auto_pr_vs_roc_gate)
      }
      th_auto <- choose_threshold(prob_val, y_val, method = method_used)
      
      if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_grid <- sort(unique(c(cutoff_grid, th_auto)))
        message(sprintf("[RF] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
      } else {
        message("[RF] Auto threshold skipped; keep original cutoff grid.")
      }
    } else {
      message("[RF] Validation overlap < 2; keep original cutoff grid.")
    }
    
    # --- evaluation
    for (cuti in cutoff_grid) {
      prob_tr_mat <- predict(model, newdata = as.data.frame(train_exp), type = "prob")
      pos_col_tr <- intersect(colnames(prob_tr_mat), c("1","positive","POS","Yes","TRUE"))
      if (!length(pos_col_tr)) pos_col_tr <- tail(colnames(prob_tr_mat), 1)
      prob_tr <- as.numeric(prob_tr_mat[, pos_col_tr[1]])
      
      train_result <- data.frame(
        predict_p     = prob_tr,
        predict_result= factor(ifelse(prob_tr > cuti, "positive", "negative")),
        real_label    = factor(ifelse(train_labels == 1, "positive", "negative"))
      )
      
      all_result <- lapply(seq_along(test_exp), function(i) {
        data_i   <- as.data.frame(test_exp[[i]])
        label_df <- labels_list[[i + 1]]
        ids_i    <- as.character(label_df[[1]])
        comd     <- intersect(rownames(data_i), ids_i)
        if (!length(comd)) {
          return(data.frame(predict_p = numeric(0),
                            predict_result = factor(character()),
                            real_label = factor(character())))
        }
        prob_mat_i <- predict(model, newdata = data_i[comd, , drop = FALSE], type = "prob")
        pos_col_i  <- intersect(colnames(prob_mat_i), c("1","positive","POS","Yes","TRUE"))
        if (!length(pos_col_i)) pos_col_i <- tail(colnames(prob_mat_i), 1)
        p  <- as.numeric(prob_mat_i[, pos_col_i[1]])
        pr <- factor(ifelse(p > cuti, "positive", "negative"))
        data.frame(predict_p = p, predict_result = pr,
                   real_label = factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1, "positive", "negative")))
      })
      
      all_result <- c(list(Train = train_result), all_result)
      
      # FIX: consistent names: A(Train), B(Val), C/D...(Test)
      n_tests <- length(test_exp)
      letters_idx <- LETTERS[2:(n_tests + 1)]
      lab_tags <- c("Val", rep("Test", max(0, n_tests - 1)))
      names(all_result) <- c("DatasetA(Train)", paste0("Dataset", letters_idx, "(", lab_tags, ")"))
      
      result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
      result_recall <- sapply(all_result, function(x) {
        prd <- as.character(x$predict_result); ref <- as.character(x$real_label)
        tp <- sum(prd == "positive" & ref == "positive")
        fn <- sum(prd == "negative" & ref == "positive")
        if ((tp + fn) == 0) 0 else tp / (tp + fn)
      })
      
      is_auto <- is.finite(th_auto) && abs(cuti - th_auto) < .Machine$double.eps^.5
      key <- if (is_auto) {
        paste0("RF (mtry=", mnumberi, ") (cutoff:auto(", formatC(th_auto, format="f", digits=4), ",", method_used, "))")
      } else {
        paste0("RF (mtry=", mnumberi, ") (cutoff:", formatC(cuti, format="f", digits=4), ")")
      }
      
      if (is.environment(collector)) .fh_collect(collector, key, result_acc, result_recall, result_FS, all_result)
      acc_local[[key]] <- result_acc
      recall_local[[key]] <- result_recall
      fs_local[[key]] <- result_FS
      summary_local[[key]] <- all_result
    }
  }
  invisible(list(acc = acc_local, recall = recall_local, fs = fs_local, summary = summary_local))
}


#' LASSO + Random Forest Interface
#'
#' @description
#' Train a Random Forest model using only features selected by LASSO.
#' Steps:
#' 1) Train with 500 trees and pick the best tree number via OOB error;
#' 2) Retrain with this best number of trees;
#' 3) On validation set B, determine an automatic threshold (Youden/F1/auto) and merge with fixed thresholds;
#' 4) Evaluate F1/Accuracy/Recall on datasets A/B/C/D...
#'
#' @param train_exp Training expression matrix (samples x features).
#' @param train_labels Binary/factor labels for training set.
#' @param test_exp list of external datasets (first = validation B, then C/D/...).
#' @param labels_list list of label data.frames (A=1, B=2, ...).
#' @param lassogene character vector of features selected by LASSO.
#' @param fold integer used for tagging result keys (not computation).
#' @param mnum candidate mtry values (default = 25%/50%/75% of \code{length(lassogene)}).
#' @param cutoff numeric thresholds (default = \code{c(0.25, 0.5, 0.75)}).
#' @param auto_th_method threshold method: "auto" (default), "f1", or "youden".
#' @param auto_imbalance_thresh imbalance threshold for auto mode (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage gate (default 0.5).
#' @param collector environment optional results collector.
#'
#' @return collector environment and an invisible list (acc, recall, fs, summary).
#'
#' @examples
#' \dontrun{
#' collector <- new.env()
#' fh_lasso_rf(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   fold         = 10,
#'   collector    = collector
#' )
#' }
#'
#' @importFrom randomForest randomForest
#' @export
fh_lasso_rf <- function(train_exp, train_labels, test_exp, labels_list,
                        lassogene, fold,
                        mnum   = c(round(stats::quantile(seq_len(length(lassogene)), probs = 0.25)),
                                   round(stats::quantile(seq_len(length(lassogene)), probs = 0.5)),
                                   round(stats::quantile(seq_len(length(lassogene)), probs = 0.75))),
                        cutoff = c(0.25, 0.5, 0.75),
                        auto_th_method = "auto",
                        auto_imbalance_thresh = 0.35,
                        auto_pr_vs_roc_gate = 0.5,
                        collector = collector) {
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  
  for (mnumberi in mnum) {
    train_expd <- as.data.frame(train_exp)[, lassogene, drop = FALSE]
    train_expd$labels <- factor(train_labels)
    
    model0 <- randomForest::randomForest(labels ~ . - labels, data = train_expd, ntree = 500, mtry = mnumberi)
    optionTrees <- which.min(model0$err.rate[, 1])
    model <- randomForest::randomForest(labels ~ . - labels, data = train_expd, ntree = optionTrees, mtry = mnumberi)
    
    # --- auto threshold on B
    cutoff_grid <- cutoff
    th_auto <- NA_real_
    method_used <- auto_th_method
    
    val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
    val_lab <- labels_list[[2]]
    ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    
    if (length(comd_val) >= 2) {
      prob_mat <- predict(model, newdata = val_exp[comd_val, , drop = FALSE], type = "prob")
      pos_col <- intersect(colnames(prob_mat), c("1","positive","POS","Yes","TRUE"))
      if (!length(pos_col)) pos_col <- tail(colnames(prob_mat), 1)
      prob_val <- as.numeric(prob_mat[, pos_col[1]])
      y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
      
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(prob_val, y_val,
                                               imbalance_thresh = auto_imbalance_thresh,
                                               pr_vs_roc_gate   = auto_pr_vs_roc_gate)
      }
      th_auto <- choose_threshold(prob_val, y_val, method = method_used)
      
      if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_grid <- sort(unique(c(cutoff_grid, th_auto)))
        message(sprintf("[lasso+RF] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
      } else {
        message("[lasso+RF] Auto threshold skipped; keep original cutoff grid.")
      }
    } else {
      message("[lasso+RF] Validation overlap < 2; keep original cutoff grid.")
    }
    
    # --- evaluation
    for (cuti in cutoff_grid) {
      prob_tr_mat <- predict(model, newdata = as.data.frame(train_exp)[, lassogene, drop = FALSE], type = "prob")
      pos_col_tr <- intersect(colnames(prob_tr_mat), c("1","positive","POS","Yes","TRUE"))
      if (!length(pos_col_tr)) pos_col_tr <- tail(colnames(prob_tr_mat), 1)
      prob_tr <- as.numeric(prob_tr_mat[, pos_col_tr[1]])
      
      train_result <- data.frame(
        predict_p     = prob_tr,
        predict_result= factor(ifelse(prob_tr > cuti, "positive", "negative")),
        real_label    = factor(ifelse(train_labels == 1, "positive", "negative"))
      )
      
      all_result <- lapply(seq_along(test_exp), function(i) {
        data_i   <- as.data.frame(test_exp[[i]])[, lassogene, drop = FALSE]
        label_df <- labels_list[[i + 1]]
        ids_i    <- as.character(label_df[[1]])
        comd     <- intersect(rownames(data_i), ids_i)
        if (!length(comd)) {
          return(data.frame(predict_p = numeric(0),
                            predict_result = factor(character()),
                            real_label = factor(character())))
        }
        prob_mat_i <- predict(model, newdata = data_i[comd, , drop = FALSE], type = "prob")
        pos_col_i  <- intersect(colnames(prob_mat_i), c("1","positive","POS","Yes","TRUE"))
        if (!length(pos_col_i)) pos_col_i <- tail(colnames(prob_mat_i), 1)
        p  <- as.numeric(prob_mat_i[, pos_col_i[1]])
        pr <- factor(ifelse(p > cuti, "positive", "negative"))
        data.frame(predict_p = p, predict_result = pr,
                   real_label = factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2]) == 1, "positive", "negative")))
      })
      
      all_result <- c(list(Train = train_result), all_result)
      
      # FIX: consistent names: A(Train), B(Val), C/D...(Test)
      n_tests <- length(test_exp)
      letters_idx <- LETTERS[2:(n_tests + 1)]
      lab_tags <- c("Val", rep("Test", max(0, n_tests - 1)))
      names(all_result) <- c("DatasetA(Train)", paste0("Dataset", letters_idx, "(", lab_tags, ")"))
      
      result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
      result_recall <- sapply(all_result, function(x) {
        prd <- as.character(x$predict_result); ref <- as.character(x$real_label)
        tp <- sum(prd == "positive" & ref == "positive")
        fn <- sum(prd == "negative" & ref == "positive")
        if ((tp + fn) == 0) 0 else tp / (tp + fn)
      })
      
      is_auto <- is.finite(th_auto) && abs(cuti - th_auto) < .Machine$double.eps^.5
      key <- if (is_auto) {
        paste0("RF+Lasso-CV:", fold, " fold (mtry=", mnumberi, ") (cutoff:auto(",
               formatC(th_auto, format="f", digits=4), ",", method_used, "))")
      } else {
        paste0("RF+Lasso-CV:", fold, " fold (mtry=", mnumberi, ") (cutoff:",
               formatC(cuti, format="f", digits=4), ")")
      }
      
      if (is.environment(collector)) .fh_collect(collector, key, result_acc, result_recall, result_FS, all_result)
      acc_local[[key]] <- result_acc
      recall_local[[key]] <- result_recall
      fs_local[[key]] <- result_FS
      summary_local[[key]] <- all_result
    }
  }
  invisible(list(acc = acc_local, recall = recall_local, fs = fs_local, summary = summary_local))
}