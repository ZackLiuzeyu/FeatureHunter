#' Linear Discriminant Analysis (LDA) interface
#'
#' @description
#' Fits a \code{MASS::lda} model on the training set; automatically selects threshold
#' (F1 / Youden) on validation set B and merges it with fixed cutoffs; evaluates
#' F1 / ACC / RECALL on A/B/C/D… datasets. Results are written into global 
#' containers: \code{all_result_FS}, \code{all_result_acc}, 
#' \code{all_result_recall}, \code{all_result_summary}.
#'
#' @param train_exp Training expression matrix
#' @param train_labels Training labels (0/1 or binary factor)
#' @param test_exp List: first element = validation set B, others = test sets C/D/…
#' @param labels_list Label list for dataset alignment (A=1, B=2, C=3, …)
#' @param cutoff Fixed cutoff set (default \code{c(0.25, 0.5, 0.75)}); merged with auto thresholds from B
#' @param collector (optional) environment collector; if supplied, each cutoff’s results are written into it
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to pick the strategy
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
#' fh_lda(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B/C/D
#'   labels_list  = labels_list,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   collector    = collector
#' )
#' }
#'
#' @importFrom MASS lda
#' @importFrom dplyr %>% mutate across
#' @export
fh_lda <- function(train_exp, 
                   train_labels, 
                   test_exp, 
                   labels_list,
                   cutoff = c(0.25, 0.5, 0.75),
                   auto_th_method = "auto",
                   auto_imbalance_thresh = 0.35,
                   auto_pr_vs_roc_gate   = 0.5,
                   collector = collector) {
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  ################3. 线性判别分析 Linear discriminant analysis #####################
  # 创建LDA模型
  train_expd <- as.data.frame(train_exp)
  ###################注意
  sum(complete.cases(train_expd)) #此处应该显示可处理的样本数，若与总样本数不同则执行以下代码
  train_expd <- train_expd %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ifelse(is.na(.), median(., na.rm=TRUE), .)))
  sum(complete.cases(train_expd)) #此处应该显示正确样本数
  ##################注意解除
  train_expd$labels <- train_labels
  model <- lda(labels ~ . - labels, data = train_expd)
  
  # ===================== ① 在 B 上自动选阈值（与固定阈值合并） ===================== #
  cutoff <- cutoff                                      # 备选阈值
  # ★ B = test_exp[[1]]，对应标签 labels_list[[2]]
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    # ★ LDA 概率用 posterior
    post_val <- predict(model, newdata = val_x)$posterior
    # ★ 取“正类”列：先按常见正类别名匹配；匹配不到就取最后一列
    pos_alias <- c("1","positive","POS","Yes","TRUE")
    pos_col <- intersect(colnames(post_val), pos_alias)
    if (length(pos_col) == 0) pos_col <- tail(colnames(post_val), 1)
    prob_val <- as.numeric(post_val[, pos_col[1]])
    y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
    
    # ★ 方法可自动判别：F1 or Youden（你已有 decide_threshold_method/choose_threshold）
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
    
    if (is.finite(th_auto) && length(th_auto) == 1   && th_auto > 0.01 && th_auto < 0.99) {
      cutoff <- unique(c(cutoff, th_auto))                          # ★ 合并为 4 个阈值
      message(sprintf("[LDA] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[LDA] Auto threshold skipped; keep fixed thresholds.")
    }
  } else {
    message("[LDA] Validation overlap < 2; keep fixed thresholds.")
  }
  
  # ===================== ② 评估：B/C/D… + A(Train)，按 cutoff 遍历 ===================== #
  for (cuti in cutoff) {
    # ---- 训练集：先拿 posterior 概率，再按阈值二分类 ----
    post_tr <- predict(model, newdata = as.data.frame(train_exp))$posterior
    pos_col_tr <- intersect(colnames(post_tr), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col_tr) == 0) pos_col_tr <- tail(colnames(post_tr), 1)
    prob_tr <- as.numeric(post_tr[, pos_col_tr[1]])
    
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(train_labels == 1, "positive", "negative"))
    
    # ---- 外部集：包含 B/C/D…（后续排序你再排除 B 就好）----
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])      # i=1 -> B, i=2 -> C, ...
      label_df <- labels_list[[i + 1]]              # A=1, B=2, C=3 …
      
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) {
        return(data.frame(predict_p = numeric(0),
                          predict_result = factor(character()),
                          real_label = factor(character())))
      }
      
      expdata    <- data_i[comd, , drop = FALSE]
      labelsdata <- label_df[match(comd, ids_i), 2]
      
      post_i <- predict(model, newdata = expdata)$posterior
      pos_col_i <- intersect(colnames(post_i), c("1","positive","POS","Yes","TRUE"))
      if (length(pos_col_i) == 0) pos_col_i <- tail(colnames(post_i), 1)
      p  <- as.numeric(post_i[, pos_col_i[1]])
      pr <- factor(ifelse(p > cuti, "positive", "negative"))
      
      out <- data.frame(predict_p = p, predict_result = pr)
      out$real_label <- factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative"))
      out
    })
    
    # ---- 把训练集结果放在最后一列（顺序：B/C/D/…/A）----
    all_result <- c(list(Train = train_result), all_result)
    nice_names <- c(
      "DatasetA(Train)",
      paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
             "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
    )
    names(all_result) <- nice_names[seq_along(all_result)]
    
    # ---- 指标 ----
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    # ---- 写出（区分 auto / fixed 便于追踪）----
    is_auto <- exists("th_auto") && is.finite(th_auto) && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("LDA (cutoff:auto(", formatC(th_auto, format = "f", digits = 4), ",", auto_th_method, "))")
    } else {
      paste0("LDA (cutoff:", formatC(cuti, format = "f", digits = 4), ")")
    }
    
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








#' @title LASSO + LDA Interface
#'
#' @description
#' This function applies LASSO-selected features (\code{lassogene}) 
#' to filter columns before training a Linear Discriminant Analysis (LDA) model.  
#' Automatic thresholds are determined on validation set B and merged with fixed thresholds.  
#' Performance metrics (F1, Accuracy, Recall) are then computed on datasets A/B/C/D…  
#' Results are written back to the result collector environment.
#'
#' @param train_exp Training expression matrix (rows = samples, cols = features).
#' @param train_labels Binary vector or factor of training labels (0/1).
#' @param test_exp A list tidyselect::where the first element is validation set B, followed by test sets C/D/…
#' @param labels_list A list of label data frames (A=1, B=2, C=3, …), used for aligning external datasets.
#' @param lassogene Character vector of LASSO-selected feature names.
#' @param fold Integer; only used for labeling in result keys (not used in computation).
#' @param cutoff Numeric vector of thresholds (default = c(0.25, 0.5, 0.75)); merged with auto thresholds from B.
#' @param collector (optional) environment collector; if supplied, each cutoff’s results are written into it.
#' @param auto_th_method Auto threshold method: \code{"youden"} / \code{"f1"} / \code{"auto"}.
#'   With \code{"auto"}, \code{decide_threshold_method()} is used to pick the strategy
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
#' fh_lasso_lda(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   fold         = 10,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   collector    = NULL
#' )
#' }
#'
#' @importFrom MASS lda
#' @importFrom dplyr %>% mutate across
#' @export
fh_lasso_lda <- function(train_exp, 
                         train_labels, 
                         test_exp, 
                         labels_list,
                         lassogene, 
                         fold, 
                         cutoff = c(0.25, 0.5, 0.75),
                         auto_th_method = "auto",
                         auto_imbalance_thresh = 0.35, 
                         auto_pr_vs_roc_gate   = 0.5, 
                         collector = collector
                         ) {
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  ################4. lasso + LDA ################
  train_expd <- as.data.frame(train_exp)[, lassogene]
  ###################注意
  sum(complete.cases(train_expd))
  train_expd <- train_expd %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),
                                ~ ifelse(is.na(.), stats::median(., na.rm = TRUE), .)))
  sum(complete.cases(train_expd))
  ##################注意解除
  
  train_expd$labels <- train_labels
  model <- lda(labels ~ . - labels, data = train_expd)
  
  cutoff <- cutoff
  # B = test_exp[[1]]；标签表用 labels_list[[2]]
  val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  if (length(comd_val) >= 2) {
    val_x  <- val_exp[comd_val, , drop = FALSE]
    post_v <- predict(model, newdata = val_x)$posterior
    # 取“正类”列名；没有匹配就取最后一列
    pos_alias <- c("1","Disease","positive","POS","Yes","TRUE")
    pos_col <- intersect(colnames(post_v), pos_alias)
    if (length(pos_col) == 0) pos_col <- tail(colnames(post_v), 1)
    prob_val <- as.numeric(post_v[, pos_col[1]])
    
    y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
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
    
    if (is.finite(th_auto) && length(th_auto) == 1 && th_auto > 0.01 && th_auto < 0.99) {
      cutoff <- unique(c(cutoff, th_auto))                        # 合并成 4 个阈值
      message(sprintf("[Lasso+LDA] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[Lasso+LDA] Auto threshold skipped; keep fixed thresholds.")
    }
  } else {
    message("[Lasso+LDA] Validation overlap < 2; keep fixed thresholds.")
  }
  
  for (cuti in cutoff) {
    # ---- 训练集：先拿 posterior 概率，再按阈值二分类 ----
    post_tr <- predict(model, newdata = as.data.frame(train_exp)[, lassogene, drop = FALSE])$posterior
    pos_col_tr <- intersect(colnames(post_tr), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col_tr) == 0) pos_col_tr <- tail(colnames(post_tr), 1)
    prob_tr <- as.numeric(post_tr[, pos_col_tr[1]])
    
    train_result <- data.frame(
      predict_p = prob_tr,
      predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
    )
    train_result$real_label <- factor(ifelse(train_labels == 1, "positive", "negative"))
    
    # ---- 外部各集：包含 B/C/D…（后续排序可排除 B）----
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i   <- as.data.frame(test_exp[[i]])[, lassogene, drop = FALSE]  # i=1 -> B
      label_df <- labels_list[[i + 1]]                                     # A=1，所以 B=2、C=3…
      
      ids_i <- as.character(label_df[[1]])
      comd  <- intersect(rownames(data_i), ids_i)
      if (length(comd) == 0) {
        return(data.frame(predict_p = numeric(0),
                          predict_result = factor(character()),
                          real_label = factor(character())))
      }
      
      expdata    <- data_i[comd, , drop = FALSE]
      labelsdata <- label_df[match(comd, ids_i), 2]
      
      post_i <- predict(model, newdata = expdata)$posterior
      pos_col_i <- intersect(colnames(post_i), c("1","positive","POS","Yes","TRUE"))
      if (length(pos_col_i) == 0) pos_col_i <- tail(colnames(post_i), 1)
      p  <- as.numeric(post_i[, pos_col_i[1]])
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
    
    # ---- 指标 ----
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){
      tp <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      fn <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      if ((tp + fn) == 0) 0 else tp / (tp + fn)
    })
    
    # ---- key：带上 fold 与 cutoff（避免覆盖）----
    is_auto <- exists("th_auto") && is.finite(th_auto) && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
    key <- if (is_auto) {
      paste0("LDA+Lasso-CV:", fold, " fold (cutoff:auto(", formatC(th_auto, format="f", digits=4), ",", auto_th_method, "))")
    } else {
      paste0("LDA+Lasso-CV:", fold, " fold (cutoff:", formatC(cuti, format="f", digits=4), ")")
    }
    
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














