#' @title Logistic Regression Interface
#'
#' @description
#' Fits a \code{glm(family = binomial)} model on the training set;  
#' automatically selects a threshold on validation set B (via \code{choose_threshold}, default "youden",
#' with optional adaptive strategy), merges it with fixed cutoffs, and evaluates
#' F1 / Accuracy / Recall across datasets A/B/C/D…  
#' Optionally writes results into a mutable \code{collector} environment 
#' (created by internal \code{.fh_new_collector()}).
#'
#' @param train_exp Training expression matrix (rows = samples, cols = features).
#' @param train_labels Binary vector or factor of training labels (0/1).
#' @param test_exp A list: first element = validation set B; remaining elements = test sets C/D/…
#' @param labels_list Label alignment list (for \code{get_labels_for_matrix()}; A=1, B=2, …).
#' @param all_labels Merged label data.frame (already built in your workflow).
#' @param cutoff Numeric vector of thresholds (default = c(0.25, 0.5, 0.75)); merged with auto threshold from B.
#' @param collector (optional) environment collector; if supplied, each cutoff’s results are written into it.
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to adaptively select
#'   based on class imbalance and PR vs ROC performance; otherwise the specified method is used.
#' @param auto_imbalance_thresh Imbalance threshold passed to \code{decide_threshold_method()} (default 0.35).
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default 0.5).
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
#' # Method 1: use a collector to accumulate results across runs
#' collector <- .fh_new_collector()
#' fh_logistic(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   all_labels   = all_labels,
#'   collector    = collector
#' )
#' results <- fh_as_lists(collector)
#' str(results$all_result_acc)
#'
#' # Method 2: without collector, directly use return value
#' out <- fh_logistic(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]]),
#'   labels_list  = labels_list,
#'   all_labels   = all_labels
#' )
#' names(out$acc)
#' }
#'
#' @importFrom dplyr %>% mutate across
#' @export
fh_logistic <- function(train_exp, train_labels,
                        test_exp,
                        labels_list,
                        all_labels,
                        cutoff = c(0.25, 0.5, 0.75),
                        auto_th_method = "auto",
                        auto_imbalance_thresh = 0.35, 
                        auto_pr_vs_roc_gate   = 0.5, 
                        collector = collector) {
  
  # —— 基础与护栏 —— #
  stopifnot(
    is.matrix(train_exp) || is.data.frame(train_exp),
    length(train_labels) == nrow(train_exp),
    is.list(test_exp),
    is.null(collector) || is.environment(collector)
  )
  cutoff <- as.numeric(cutoff)
  
  # —— 训练 GLM —— #
  train_expd <- as.data.frame(train_exp)
  
  # 缺失填补（数值列用中位数）
  invisible(sum(complete.cases(train_expd)))
  train_expd <- train_expd %>%
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.numeric),
      ~ ifelse(is.na(.), stats::median(., na.rm = TRUE), .)
    ))
  invisible(sum(complete.cases(train_expd)))
  
  train_expd$labels <- train_labels
  # 静默 GLM 拟合相关警告
  model <- suppressWarnings(
    stats::glm(labels ~ . - labels, family = "binomial", data = train_expd)
  )
  
  # —— 验证集B自动阈值 —— #
  val_exp <- test_exp[[1]]
  val_lab <- labels_list[[2]]
  comd_val <- intersect(rownames(val_exp), val_lab$V1)
  
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    y_val  <- as_binary01(get_labels_for_matrix(val_exp, labels_list))
    # 静默预测阶段警告
    prob_val <- as.numeric(suppressWarnings(
      stats::predict(model, newdata = as.data.frame(val_x), type = "response")
    ))
    
    # 可选：自动决定策略（如 F1 vs Youden）
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
      message(sprintf("[LR] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[LR] Auto threshold skipped; keep original cutoff grid.")
    }
  } else {
    message("[LR] Validation overlap < 2; keep original cutoff grid.")
  }
  
  # —— 评估：A（Train）+ test_exp（B/C/...） —— #
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  for (cuti in cutoff) {
    # A: Train（静默预测）
    train_prob <- as.data.frame(suppressWarnings(
      stats::predict(model, type = "response")
    ))
    train_res  <- data.frame(
      predict_p      = train_prob[, 1],
      predict_result = factor(ifelse(train_prob[, 1] > cuti, "positive", "negative")),
      real_label     = factor(ifelse(train_labels == 1, "positive", "negative"))
    )
    
    # B/C/... from test_exp（静默预测）
    all_result <- lapply(test_exp, function(data, model, labelsdata) {
      data <- as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labs <- labelsdata[comd, 2]
      expd <- data[comd, , drop = FALSE]
      prob <- as.data.frame(suppressWarnings(
        stats::predict(model, expd, type = "response")
      ))
      data.frame(
        predict_p      = prob[, 1],
        predict_result = factor(ifelse(prob[, 1] > cuti, "positive", "negative")),
        real_label     = factor(ifelse(labs == 1, "positive", "negative"))
      )
    }, model = model, labelsdata = all_labels)
    
    # ① 把 Train 放最前
    all_result <- c(list(train_res), all_result)
    
    # ② 统一命名
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    # ③ 指标
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) {
      mean(as.character(x$predict_result) == as.character(x$real_label))
    })
    result_recall <- sapply(all_result, function(x) {
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    key <- sprintf("LR (cutoff:%.3f)", cuti)
    
    # —— 若提供 collector，就即时写入 —— #
    if (is.environment(collector)) {
      .fh_collect(
        collector, key,
        acc     = result_acc,
        recall  = result_recall,
        fs      = result_FS,
        summary = all_result
      )
    }
    
    # —— 同时把本轮结果累积到返回值 —— #
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

