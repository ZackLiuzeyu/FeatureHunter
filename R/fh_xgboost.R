#' XGBoost Interface
#'
#' @description
#' Train a basic XGBoost model (objective = \code{"binary:logistic"}). Thresholds are
#' automatically selected on validation set B and merged with fixed thresholds.
#' Performance (F1/ACC/RECALL) is computed on A/B/C/D… datasets. Variable importance
#' is also extracted.
#'
#' @param train_exp Training expression matrix (samples × features).
#' @param train_labels Binary vector of training labels (0/1).
#' @param test_exp A list: first = validation set B, others = test sets.
#' @param labels_list A list of label data.frames (A=1, B=2, …; first two columns = ID, label).
#' @param cutoff Numeric vector of fixed thresholds (default = \code{c(0.25, 0.5, 0.75)}).
#' @param nround,max_depth,eta XGBoost hyperparameters (defaults: 100, 6, 0.5).
#' @param auto_th_method Thresholding method: \code{"youden"}/\code{"f1"}/\code{"auto"}.
#'   When \code{"auto"}, \code{decide_threshold_method()} is used.
#' @param auto_imbalance_thresh Imbalance threshold for auto mode (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC gate for auto mode (default 0.5).
#' @param collector Optional environment to collect results; if supplied, results are
#'   written into it on the fly via \code{.fh_collect()}.
#'
#' @return
#' Collector environment and an invisible list of results
#' (acc, recall, fs, summary, importance).
#'
#' @examples
#' \dontrun{
#' # create a collector
#' collector <- fh_new_collector()
#'
#' out <- fh_xgboost(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   auto_th_method = "auto",
#'   collector   = collector
#' )
#'
#' # convert to plain lists if needed
#' lists <- fh_as_lists(collector)
#' }
#'
#' @importFrom xgboost xgb.DMatrix xgboost xgb.importance 
#' @export
fh_xgboost <- function(train_exp, train_labels, test_exp, labels_list,
                       cutoff = c(0.25, 0.5, 0.75),
                       nround = 100, max_depth = 6, eta = 0.5,
                       auto_th_method = "auto",
                       auto_imbalance_thresh = 0.35,
                       auto_pr_vs_roc_gate   = 0.5,
                       collector = collector) {
  cutoff_now <- cutoff
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_exp), label = as_binary01(train_labels))
  model  <- xgboost::xgboost(data = dtrain,
                             objective = "binary:logistic",
                             nround    = nround,
                             max_depth = max_depth,
                             eta       = eta,
                             verbose   = 0)
  
  # auto threshold on B
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  if (length(comd_val) >= 2) {
    val_x    <- as.matrix(val_exp[comd_val, , drop = FALSE])
    prob_val <- as.numeric(predict(model, val_x))
    y_val    <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    method_used <- auto_th_method
    if (identical(auto_th_method, "auto")) {
      method_used <- decide_threshold_method(prob_val, y_val,
                                             imbalance_thresh = auto_imbalance_thresh,
                                             pr_vs_roc_gate   = auto_pr_vs_roc_gate)
    }
    th_auto <- choose_threshold(prob_val, y_val, method = method_used)
    if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
      cutoff_now <- unique(c(cutoff_now, th_auto))
      message(sprintf("[XGB] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    }else {
      message("[XGB] Auto threshold failed; keep original cutoff grid.")
    }
  }else {
    message("[XGB] Validation overlap < 2; keep original cutoff grid.")
  }
  
  acc_local <- list(); recall_local <- list(); fs_local <- list()
  summary_local <- list(); importance_local <- list()
  
  for (cuti in cutoff_now) {
    im <- as.data.frame(xgboost::xgb.importance(model = model))
    rownames(im) <- im$Feature
    
    prob_tr <- as.numeric(predict(model, as.matrix(train_exp)))
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative")),
      real_label = factor(ifelse(train_labels == 1, "positive", "negative"))
    )
    
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i <- as.data.frame(test_exp[[i]])
      label_df <- labels_list[[i+1]]
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (!length(comd)) return(data.frame(predict_p=numeric(),
                                           predict_result=factor(character()),
                                           real_label=factor(character())))
      expdata <- as.matrix(data_i[comd, , drop=FALSE])
      labelsdata <- label_df[match(comd, ids_i), 2]
      p  <- as.numeric(predict(model, expdata))
      pr <- factor(ifelse(p > cuti, "positive", "negative"))
      data.frame(
        predict_p = p,
        predict_result = pr,
        real_label = factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative"))
      )
    })
    all_result <- c(list(Train = train_result), all_result)
    names(all_result) <- c("DatasetA(Train)",
                           paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
                                  "(", c("Val", rep("Test", length(test_exp) - 1)), ")"))
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    key <- paste0("XGBoost (cutoff:", formatC(cuti, format="f", digits=4), ")")
    
    if (is.environment(collector)) {
      .fh_collect(collector, key,
                  acc = result_acc, recall = result_recall,
                  fs = result_FS, summary = all_result)
    }
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  invisible(list(acc = acc_local,
                 recall = recall_local,
                 fs = fs_local,
                 summary = summary_local))
}



#' LASSO (default) + XGBoost Interface
#'
#' @description
#' Train an XGBoost model (objective = \code{"binary:logistic"}) using only
#' features in \code{lassogene}. Thresholds are selected on validation set B and
#' merged with fixed thresholds. Performance (F1/ACC/RECALL) is computed on A/B/C/D…
#' datasets. Variable importance is also extracted.
#'
#' @param train_exp Training expression matrix (samples × features).
#' @param train_labels Binary training labels (0/1).
#' @param test_exp A list of external datasets; the first is validation set B.
#' @param labels_list Label data.frames list (A=1, B=2, …; first two cols = ID, label).
#' @param lassogene Character vector of feature names selected by LASSO.
#' @param cutoff Numeric vector of fixed thresholds (default \code{c(0.25, 0.5, 0.75)}).
#' @param nround,max_depth,eta XGBoost hyperparameters (defaults: 100, 6, 0.5).
#' @param auto_th_method Thresholding method: \code{"youden"}/\code{"f1"}/\code{"auto"}.
#' @param auto_imbalance_thresh Imbalance threshold for auto mode (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC gate for auto mode (default 0.5).
#' @param collector Optional environment to collect results.
#'
#' @return
#' Collector environment and an invisible list of results
#' (acc, recall, fs, summary, importance).
#'
#' @examples
#' \dontrun{
#' collector <- fh_new_collector()
#'
#' out <- fh_lasso_xgboost_default(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   auto_th_method = "auto",
#'   collector   = collector
#' )
#'
#' lists <- fh_as_lists(collector)
#' }
#'
#' @importFrom xgboost xgb.DMatrix xgboost xgb.importance 
#' @export
fh_lasso_xgboost_default <- function(train_exp, train_labels, test_exp, labels_list,
                                     lassogene,
                                     cutoff = c(0.25, 0.5, 0.75),
                                     nround = 100, max_depth = 6, eta = 0.5,
                                     auto_th_method = "auto",
                                     auto_imbalance_thresh = 0.35,
                                     auto_pr_vs_roc_gate   = 0.5,
                                     collector = collector) {
  dtrain <- xgboost::xgb.DMatrix(
    data  = as.matrix(train_exp[, lassogene, drop = FALSE]),
    label = as_binary01(train_labels)
  )
  model <- xgboost::xgboost(
    data      = dtrain,
    objective = "binary:logistic",
    nround    = nround,
    max_depth = max_depth,
    eta       = eta,
    verbose   = 0
  )
  
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  if (length(comd_val) >= 2) {
    val_x    <- as.matrix(val_exp[comd_val, , drop = FALSE])
    prob_val <- as.numeric(predict(model, val_x))
    y_val    <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    method_used <- auto_th_method
    if (identical(auto_th_method, "auto")) {
      method_used <- decide_threshold_method(prob_val, y_val,
                                             imbalance_thresh = auto_imbalance_thresh,
                                             pr_vs_roc_gate   = auto_pr_vs_roc_gate)
    }
    th_auto <- choose_threshold(prob_val, y_val, method = method_used)
    if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
      cutoff_now <- unique(c(cutoff_now, th_auto))
      message(sprintf("[lasso+XGBoost default] Auto threshold on B = %.4f (method=%s)",th_auto, method_used))
    }else {
      message("[lasso+XGBoost default] Auto threshold failed; keep original cutoff grid.")
    }
  } else {
    message("[lasso+XGBoost default] Validation overlap < 2; keep original cutoff grid.")
  }
  
  acc_local <- list(); recall_local <- list(); fs_local <- list()
  summary_local <- list(); importance_local <- list()
  
  for (cuti in cutoff_now) {
    im <- as.data.frame(xgboost::xgb.importance(model = model))
    rownames(im) <- im$Feature
    
    prob_tr <- as.numeric(predict(model, as.matrix(train_exp[, lassogene, drop = FALSE])))
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
      if (!length(comd)) return(data.frame(predict_p=numeric(),
                                           predict_result=factor(character()),
                                           real_label=factor(character())))
      expdata    <- as.matrix(data_i[comd, , drop = FALSE])
      labelsdata <- label_df[match(comd, ids_i), 2]
      p  <- as.numeric(predict(model, expdata))
      pr <- factor(ifelse(p > cuti, "positive", "negative"))
      data.frame(
        predict_p = p,
        predict_result = pr,
        real_label = factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative"))
      )
    })
    
    all_result <- c(list(Train = train_result), all_result)
    names(all_result) <- c("DatasetA(Train)",
                           paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
                                  "(", c("Val", rep("Test", length(test_exp) - 1)), ")"))
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    key <- paste0("lasso+XGBoost default (cutoff:", formatC(cuti, format="f", digits=4), ")")
    
    if (is.environment(collector)) {
      .fh_collect(collector, key,
                  acc = result_acc, recall = result_recall,
                  fs = result_FS, summary = all_result)
    }
    acc_local[[key]]     <- result_acc
    recall_local[[key]]  <- result_recall
    fs_local[[key]]      <- result_FS
    summary_local[[key]] <- all_result
  }
  
  invisible(list(acc = acc_local,
                 recall = recall_local,
                 fs = fs_local,
                 summary = summary_local))
}

#' LASSO (best) + XGBoost Interface (caret CV)
#'
#' @description
#' Perform hyperparameter tuning with \code{caret::train(method="xgbTree")} on the
#' \code{lassogene} subset, take the \code{finalModel}, select thresholds automatically
#' on validation set B (or fixed), and evaluate on A/B/C/D… datasets. Variable
#' importance is also extracted.
#'
#' @param train_exp Training expression matrix (samples × features).
#' @param train_labels Binary training labels (0/1).
#' @param test_exp A list of external datasets; the first is validation set B.
#' @param labels_list Label data.frames list (A=1, B=2, …; first two cols = ID, label).
#' @param lassogene Character vector of feature names selected by LASSO.
#' @param fold Number of CV folds used in caret (default 10).
#' @param nrounds Maximum boosting rounds for the tuning grid (default 1000).
#' @param cutoff Numeric vector of fixed thresholds (default \code{c(0.25, 0.5, 0.75)}).
#' @param auto_th_method Thresholding method: \code{"youden"}/\code{"f1"}/\code{"auto"}.
#' @param auto_imbalance_thresh Imbalance threshold for auto mode (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC gate for auto mode (default 0.5).
#' @param collector Optional environment to collect results.
#'
#' @return
#' Collector environment and an invisible list of results
#' (acc, recall, fs, summary, importance).
#'
#' @examples
#' \dontrun{
#' collector <- fh_new_collector()
#'
#' out <- fh_lasso_xgboost_best(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   fold         = 10,
#'   auto_th_method = "auto",
#'   collector   = collector
#' )
#'
#' lists <- fh_as_lists(collector)
#' }
#'
#' @importFrom caret train trainControl
#' @importFrom xgboost xgb.importance 
#' @export
fh_lasso_xgboost_best <- function(train_exp, train_labels, test_exp, labels_list,
                                  lassogene, fold = 10, nrounds = 1000,
                                  cutoff = c(0.25, 0.5, 0.75),
                                  auto_th_method = "auto",
                                  auto_imbalance_thresh = 0.35,
                                  auto_pr_vs_roc_gate   = 0.5,
                                  collector = collector) {
  
  train_expd <- as.data.frame(train_exp)[, lassogene, drop = FALSE]
  train_expd$labels <- factor(train_labels)
  
  tune_grid <- expand.grid(
    nrounds = seq(from = 50, to = nrounds, by = 50),
    eta = c(0.025, 0.05, 0.1, 0.3),
    max_depth = c(2, 3, 4, 5, 6),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
  )
  fitControl <- caret::trainControl(method = "cv",
                                    number = fold,
                                    verboseIter = FALSE,
                                    search = "grid")
  
  xgboost::xgb.set.config(verbosity = 0)
  caret_xgb <- caret::train(labels ~ . - labels,
                            data = train_expd,
                            method = "xgbTree",
                            trControl = fitControl,
                            tuneGrid = tune_grid)
  
  model <- caret_xgb$finalModel
  
  # importance snapshot (store once under a stable key)
  im <- as.data.frame(xgboost::xgb.importance(model = model))
  rownames(im) <- im$Feature
  
  cutoff_now <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  if (length(comd_val) >= 2) {
    val_x    <- as.matrix(val_exp[comd_val, , drop = FALSE])
    prob_val <- as.numeric(predict(model, val_x))
    y_val    <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    method_used <- auto_th_method
    if (identical(auto_th_method, "auto")) {
      method_used <- decide_threshold_method(prob_val, y_val,
                                             imbalance_thresh = auto_imbalance_thresh,
                                             pr_vs_roc_gate   = auto_pr_vs_roc_gate)
    }
    th_auto <- choose_threshold(prob_val, y_val, method = method_used)
    if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
      cutoff_now <- unique(c(cutoff_now, th_auto))
      message(sprintf("[lasso+XGBoost CV best] Auto threshold on B = %.4f (method=%s)",th_auto, method_used))
    }else {
      message("[lasso+XGBoost CV best] Auto threshold failed; keep original cutoff grid.")
    }
  } else {
    message("[lasso+XGBoost CV best] Validation overlap < 2; keep original cutoff grid.")
  }
  
  acc_local <- list(); recall_local <- list(); fs_local <- list()
  summary_local <- list(); importance_local <- list()
  
  for (cuti in cutoff_now) {
    prob_tr <- as.numeric(predict(model, as.matrix(as.data.frame(train_exp)[, lassogene, drop = FALSE])))
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
      if (!length(comd)) return(data.frame(predict_p=numeric(),
                                           predict_result=factor(character()),
                                           real_label=factor(character())))
      expdata    <- as.matrix(data_i[comd, , drop = FALSE])
      labelsdata <- label_df[match(comd, ids_i), 2]
      p  <- as.numeric(predict(model, expdata))
      pr <- factor(ifelse(p > cuti, "positive", "negative"))
      data.frame(
        predict_p = p,
        predict_result = pr,
        real_label = factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative"))
      )
    })
    
    all_result <- c(list(Train = train_result), all_result)
    names(all_result) <- c("DatasetA(Train)",
                           paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
                                  "(", c("Val", rep("Test", length(test_exp) - 1)), ")"))
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    key <- paste0("XGBoost-CV:", fold, " fold + Lasso(best) (cutoff:",
                  formatC(cuti, format = "f", digits = 4), ")")
    
    if (is.environment(collector)) {
      .fh_collect(collector, key,
                  acc = result_acc, recall = result_recall,
                  fs = result_FS, summary = all_result)
    }
    acc_local[[key]]        <- result_acc
    recall_local[[key]]     <- result_recall
    fs_local[[key]]         <- result_FS
    summary_local[[key]]    <- all_result
  }
  
  invisible(list(acc = acc_local,
                 recall = recall_local,
                 fs = fs_local,
                 summary = summary_local))
}

