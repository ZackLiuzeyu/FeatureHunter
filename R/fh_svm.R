#' Support Vector Machine (e1071) Interface
#'
#' @description
#' Train SVM classifiers with specified kernels (default: "linear", "polynomial", "radial")
#' using \code{e1071::svm(probability = TRUE)}. On validation set B, choose a threshold
#' (F1/Youden; supports an automatic policy) and merge it with fixed thresholds
#' (\code{c(0.25, 0.5, 0.75)}). Evaluate F1 / Accuracy / Recall on datasets A/B/C/D…
#' Results can be written into a \emph{collector} environment and are also returned
#' as plain lists.
#'
#' @param train_exp data.frame; training expression matrix (samples × features).
#' @param train_labels vector; training labels (0/1 or convertible to 0/1; coerced to factor for training).
#' @param test_exp list; external datasets (first = validation B, then C/D/…).
#' @param labels_list list; label data.frames (A=1, B=2, …; first column = ID, second = label).
#' @param kernel_all character; kernels to use (default \code{c("linear","polynomial","radial")}).
#' @param cutoff numeric; fixed thresholds (default \code{c(0.25, 0.5, 0.75)}), merged with the auto threshold from B.
#' @param auto_th_method character; \code{"youden"} / \code{"f1"} / \code{"auto"} (default \code{"auto"}).
#'   With \code{"auto"}, \code{decide_threshold_method()} is used based on class imbalance and PR-vs-ROC gain.
#' @param auto_imbalance_thresh numeric; imbalance threshold for \code{decide_threshold_method()} (default \code{0.35}).
#' @param auto_pr_vs_roc_gate numeric; PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default \code{0.5}).
#' @param collector environment (optional); if supplied, each cutoff’s results are appended via \code{.fh_collect()}.
#'
#' @return
#' Collector environment and an invisible list of results (acc, recall, fs, summary).
#' Specifically, returns \code{invisible(list(acc=..., recall=..., fs=..., summary=...))}.
#' When \code{collector} is provided, results are also written into it.
#'
#' @examples
#' \dontrun{
#' # Prepare a collector once:
#' collector <- new.env(parent = emptyenv())
#'
#' fh_svm(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   kernel_all   = c("linear","polynomial","radial"),
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   auto_th_method = "auto",
#'   collector    = collector
#' )
#'
#' # Convert the collector to simple lists if needed:
#' res_lists <- fh_as_lists(collector)
#' }
#'
#' @importFrom e1071 svm
#' @importFrom dplyr mutate across
#' @export
fh_svm <- function(train_exp, train_labels, test_exp, labels_list,
                   kernel_all = c("linear", "polynomial", "radial"),
                   cutoff = c(0.25, 0.5, 0.75),
                   auto_th_method = "auto",
                   auto_imbalance_thresh = 0.35,
                   auto_pr_vs_roc_gate   = 0.5,
                   collector = collector) {
  
  # helper: extract P(positive) from e1071::predict(..., probability=TRUE)
  .get_pos_prob <- function(pred_obj) {
    probs <- tryCatch(attr(pred_obj, "probabilities"), error = function(e) NULL)
    if (is.null(probs)) return(NULL)
    cn <- colnames(probs)
    pick <- which(cn %in% c("1","pos","positive","YES","True","TRUE"))
    if (length(pick) == 0) pick <- length(cn)
    as.numeric(probs[, pick[1]])
  }
  
  # local accumulators (also returned invisibly)
  acc_local        <- list()
  recall_local     <- list()
  fs_local         <- list()
  summary_local    <- list()
  importance_local <- list()   # ★ 攒在本地返回，不再写全局
  
  for (kerneli in kernel_all) {
    # impute training NA by column median
    train_exp_fix <- as.data.frame(train_exp)
    train_exp_fix <- train_exp_fix |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),
                                  ~ ifelse(is.na(.), stats::median(., na.rm = TRUE), .)))
    train_exp_fix <- as.matrix(train_exp_fix)
    
    # train SVM with probability=TRUE
    train_expd <- as.data.frame(train_exp_fix)
    train_expd$labels <- factor(train_labels)
    model <- e1071::svm(
      labels ~ . - labels,
      data = train_expd,
      kernel = kerneli,
      probability = TRUE
    )
    
    # optional: a simple importance proxy（保留你的写法，但不再写全局）
    im <- tryCatch({
      M <- t(model$coefs) %*% model$SV
      M <- t(as.data.frame(abs(M)))
      colnames(M) <- "importance"
      M
    }, error = function(e) NULL)
    if (!is.null(im)) {
      importance_local[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- im
    }
    
    cutoff_now <- cutoff
    th_auto <- NA_real_
    method_used <- auto_th_method
    
    # ===== B: auto threshold (merge with fixed) =====
    val_exp <- as.data.frame(test_exp[[1]])
    val_exp <- val_exp |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),
                                  ~ ifelse(is.na(.), stats::median(., na.rm = TRUE), .)))
    val_lab <- labels_list[[2]]
    ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    
    if (length(comd_val) >= 2) {
      val_x  <- val_exp[comd_val, , drop = FALSE]
      pred_val <- stats::predict(model, as.data.frame(val_x), probability = TRUE)
      prob_val <- .get_pos_prob(pred_val)
      
      if (!is.null(prob_val) && length(unique(prob_val[is.finite(prob_val)])) >= 2) {
        y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
        
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
          message(sprintf("[SVM %s] Auto threshold on B = %.4f (method=%s)",
                          kerneli, th_auto, method_used))
        } else {
          message(sprintf("[SVM %s] Auto threshold skipped; keep fixed thresholds.", kerneli))
        }
      } else {
        message(sprintf("[SVM %s] B has insufficient probability variation; keep fixed thresholds.", kerneli))
      }
    } else {
      message(sprintf("[SVM %s] Validation overlap < 2; keep fixed thresholds.", kerneli))
    }
    
    # ===== evaluate for each cutoff =====
    for (cuti in cutoff_now) {
      # train-set probs
      pred_tr <- stats::predict(model, as.data.frame(train_exp_fix), probability = TRUE)
      prob_tr <- .get_pos_prob(pred_tr)
      if (is.null(prob_tr)) prob_tr <- as.numeric(as.character(pred_tr))
      
      train_result <- data.frame(
        predict_p = prob_tr,
        predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
      )
      train_result$real_label <- factor(ifelse(as_binary01(train_labels) == 1, "positive", "negative"))
      
      # external sets B/C/D/...
      all_result <- lapply(seq_along(test_exp), function(i) {
        data_i <- as.data.frame(test_exp[[i]])
        data_i <- data_i |>
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),
                                      ~ ifelse(is.na(.), stats::median(., na.rm = TRUE), .)))
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
        expdata    <- data_i[comd, , drop = FALSE]
        labelsdata <- label_df[match(comd, ids_i), 2]
        
        pred_i <- stats::predict(model, as.data.frame(expdata), probability = TRUE)
        p <- .get_pos_prob(pred_i)
        if (is.null(p)) p <- as.numeric(as.character(pred_i))
        pr <- factor(ifelse(p > cuti, "positive", "negative"))
        
        out <- data.frame(predict_p = p, predict_result = pr)
        out$real_label <- factor(ifelse(as_binary01(labelsdata) == 1, "positive", "negative"))
        out
      })
      
      # naming
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
      is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
      key <- if (is_auto) {
        paste0("SVM-default (kernel: ", kerneli, ") (cutoff:auto(",
               formatC(th_auto, format = "f", digits = 4), ",", method_used, "))")
      } else {
        paste0("SVM-default (kernel: ", kerneli, ") (cutoff:",
               formatC(cuti, format = "f", digits = 4), ")")
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
      
      # also accumulate locally (returned invisibly)
      acc_local[[key]]     <- result_acc
      recall_local[[key]]  <- result_recall
      fs_local[[key]]      <- result_FS
      summary_local[[key]] <- all_result
    }
  }
  
  invisible(list(
    acc        = acc_local,
    recall     = recall_local,
    fs         = fs_local,
    summary    = summary_local,
    importance = importance_local   # ★ 本地返回的“重要性”
  ))
}


#' Optimal Support Vector Machine (SVM, CV best) Interface
#'
#' @description
#' Implements two parts:
#' - **Part A**: Train `e1071::svm(probability=TRUE)` on **LASSO-selected features**
#'   for kernels in `kernel_all` ("linear", "polynomial", "radial"). On validation set B,
#'   an automatic threshold is chosen (F1/Youden via `auto_th_method` decided by
#'   `decide_threshold_method`) and merged with fixed thresholds (`c(0.25, 0.5, 0.75)`),
#'   then evaluate on A/B/C/D… datasets.
#' - **Part B**: Use `caret::train` to perform CV-tuning for
#'   `svmLinear` / `svmRadial` / `svmPoly`, following the same
#'   “auto threshold on B + evaluation” workflow.
#'
#' Metrics can be written into a provided collector environment and are also
#' accumulated locally and returned invisibly.
#'
#' @param train_exp data.frame. Training expression matrix (samples × features).
#' @param train_labels vector. Training labels (0/1 or convertible to 0/1).
#' @param test_exp list. External datasets (first = validation B, then C/D/…).
#' @param labels_list list. Label data.frames (A=1, B=2, …; first col = ID, second col = label).
#' @param fold integer. Number of CV folds (used in caret and stored in result keys).
#' @param lassogene character. LASSO-selected feature names (used in Part A).
#' @param kernel_all character. Kernels to iterate (default `c("linear","polynomial","radial")`).
#' @param cutoff numeric. Fixed thresholds (default `c(0.25, 0.5, 0.75)`), merged with auto thresholds from B.
#' @param collector environment (optional). If provided, results are written into it via `.fh_collect(...)`.
#' @param auto_th_method character; \code{"youden"} / \code{"f1"} / \code{"auto"} (default \code{"auto"}).
#'   With \code{"auto"}, \code{decide_threshold_method()} is used based on class imbalance and PR-vs-ROC gain.
#' @param auto_imbalance_thresh numeric; imbalance threshold for \code{decide_threshold_method()} (default \code{0.35}).
#' @param auto_pr_vs_roc_gate numeric; PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default \code{0.5}).
#'
#'
#' @return
#' Writes to the `collector` environment (if provided) and invisibly returns a
#' list of results with four elements: `acc`, `recall`, `fs`, and `summary`.
#' Each element is a named list keyed by the constructed model tag.
#'
#' @examples
#' \dontrun{
#' # Prepare containers (collector is optional)
#' collector <- new.env(parent = emptyenv())
#'
#' out <- fh_svm_best(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   fold         = 10,
#'   lassogene    = lassogene,
#'   collector    = collector
#' )
#'
#' # Inspect returned results
#' names(out$acc)
#' }
#'
#' @importFrom e1071 svm
#' @importFrom caret train trainControl
#' @importFrom dplyr mutate across
#' @export
fh_svm_best <- function(train_exp, train_labels, test_exp, labels_list,
                        fold = 10, lassogene,
                        kernel_all = c("linear", "polynomial", "radial"),
                        cutoff = c(0.25, 0.5, 0.75),
                        auto_th_method = "auto",
                        auto_imbalance_thresh = 0.35,
                        auto_pr_vs_roc_gate   = 0.5,
                        collector = collector) {
  
  # local accumulators (returned invisibly)
  acc_local        <- list()
  recall_local     <- list()
  fs_local         <- list()
  summary_local    <- list()
  importance_local <- list()   # ← 存放本地的重要性，不再写全局
  
  # helpers: pick positive-class probability
  .get_pos_prob_mat <- function(probs) {
    if (is.null(probs)) return(NULL)
    cn <- colnames(probs)
    pick <- which(cn %in% c("1","Disease","positive","POS","Yes","TRUE","True"))
    if (length(pick) == 0) pick <- length(cn)
    as.numeric(probs[, pick[1]])
  }
  .get_pos_prob_df <- function(df) {
    if (is.null(df) || !ncol(df)) return(NULL)
    cn <- colnames(df)
    pick <- which(cn %in% c("1","Disease","positive","POS","Yes","TRUE","True","pos"))
    if (length(pick) == 0) pick <- ncol(df)
    as.numeric(df[[pick[1]]])
  }
  
  ################ Part A: e1071::svm on LASSO subset ################
  for (kerneli in kernel_all) {
    train_expd <- as.data.frame(train_exp)[, lassogene, drop = FALSE]
    train_expd$labels <- factor(train_labels)
    model <- e1071::svm(labels ~ . - labels,
                        data = train_expd,
                        kernel = kerneli,
                        probability = TRUE)
    
    # importance (unchanged logic) → 只存到本地 importance_local
    im <- tryCatch({
      M <- t(model$coefs) %*% model$SV
      M <- t(as.data.frame(abs(M))); colnames(M) <- "importance"; M
    }, error = function(e) NULL)
    if (!is.null(im)) {
      importance_local[[paste0("SVM-default+Lasso-CV:", fold, " fold (kernel: ", kerneli, ")")]] <- im
    }
    
    # auto threshold on B
    cutoff_now <- cutoff
    val_exp <- as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE]
    val_lab <- labels_list[[2]]
    ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    th_auto <- NA_real_; auto_th_method <- "youden"
    if (length(comd_val) >= 2) {
      val_x <- val_exp[comd_val, , drop = FALSE]
      pred_val <- predict(model, newdata = val_x, probability = TRUE)
      prob_val <- .get_pos_prob_mat(attr(pred_val, "probabilities"))
      if (!is.null(prob_val) && length(unique(prob_val[is.finite(prob_val)])) >= 2) {
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
        if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
          cutoff_now <- unique(c(cutoff_now, th_auto))
          message(sprintf("[SVM+CV best %s] Auto threshold on B = %.4f (%s)",kerneli, th_auto, auto_th_method))
        }else {
          message(sprintf("[SVM+CV best %s] Auto threshold skipped; keep fixed thresholds.", kerneli))
        }
      }else {
        message(sprintf("[SVM+CV best %s] B has insufficient probability variation; keep fixed thresholds.", kerneli))
      }
    }else {
      message(sprintf("[SVM %s] Validation overlap < 2; keep fixed thresholds.", kerneli))
    }
    
    # iterate thresholds: evaluate A/B/C/D…
    for (cuti in cutoff_now) {
      pred_tr <- predict(model, newdata = as.data.frame(train_exp)[, lassogene, drop = FALSE],
                         probability = TRUE)
      prob_tr <- .get_pos_prob_mat(attr(pred_tr, "probabilities"))
      if (is.null(prob_tr)) prob_tr <- as.numeric(as.character(pred_tr))
      
      train_result <- data.frame(
        predict_p = prob_tr,
        predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
      )
      train_result$real_label <- factor(ifelse(train_labels == 1, "positive", "negative"))
      
      all_result <- lapply(seq_along(test_exp), function(i){
        data_i <- as.data.frame(test_exp[[i]])[, lassogene, drop = FALSE]
        label_df <- labels_list[[i + 1]]
        ids_i <- as.character(label_df[[1]])
        comd  <- intersect(rownames(data_i), ids_i)
        if (length(comd) == 0) {
          return(data.frame(predict_p = numeric(0),
                            predict_result = factor(character()),
                            real_label = factor(character())))
        }
        expdata <- data_i[comd, , drop = FALSE]
        pred_i  <- predict(model, newdata = expdata, probability = TRUE)
        p       <- .get_pos_prob_mat(attr(pred_i, "probabilities"))
        if (is.null(p)) p <- as.numeric(as.character(pred_i))
        pr      <- factor(ifelse(p > cuti, "positive", "negative"))
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
      
      result_FS  <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
      result_recall <- sapply(all_result, function(x){
        tp <- sum(as.character(x$predict_result)=="positive" & as.character(x$real_label)=="positive")
        fn <- sum(as.character(x$predict_result)=="negative" & as.character(x$real_label)=="positive")
        if ((tp+fn)==0) 0 else tp/(tp+fn)
      })
      
      is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
      key <- if (is_auto) {
        paste0("SVM-default+Lasso-CV:", fold, " fold (kernel: ", kerneli, ") (cutoff:auto(",
               formatC(th_auto, format="f", digits=4), ",", auto_th_method, "))")
      } else {
        paste0("SVM-default+Lasso-CV:", fold, " fold (kernel: ", kerneli, ") (cutoff:",
               formatC(cuti, format="f", digits=4), ")")
      }
      
      # write to collector if provided + accumulate locally
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
  }
  
  ################ Part B: caret::train SVMs ################
  
  # ---- 1) svmLinear ----
  {
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(ifelse(as_binary01(train_labels) == 1, "pos", "neg"),
                                levels = c("neg","pos"))
    ctrl <- caret::trainControl(method = "cv", number = fold, classProbs = TRUE)
    model <- caret::train(labels ~ . - labels, data = train_expd,
                          method = "svmLinear", trControl = ctrl)
    
    cutoff_now <- cutoff
    val_exp <- as.data.frame(test_exp[[1]])
    val_lab <- labels_list[[2]]; ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    th_auto <- NA_real_; auto_th_method <- "youden"
    if (length(comd_val) >= 2) {
      val_x <- val_exp[comd_val, , drop = FALSE]
      prob_val_df <- predict(model, newdata = val_x, type = "prob")
      prob_val <- .get_pos_prob_df(prob_val_df)
      if (!is.null(prob_val) && length(unique(prob_val[is.finite(prob_val)])) >= 2) {
        y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
        auto_th_method <- decide_threshold_method(prob_val, y_val)
        th_auto <- choose_threshold(prob_val, y_val, method = auto_th_method)
        if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) cutoff_now <- unique(c(cutoff_now, th_auto))
      }
    }
    
    for (cuti in cutoff_now) {
      prob_tr_df <- predict(model, newdata = as.data.frame(train_exp), type = "prob")
      prob_tr <- .get_pos_prob_df(prob_tr_df)
      train_result <- data.frame(
        predict_p = prob_tr,
        predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
      )
      train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
      
      all_result <- lapply(seq_along(test_exp), function(i){
        data_i <- as.data.frame(test_exp[[i]])
        label_df <- labels_list[[i+1]]
        ids_i <- as.character(label_df[[1]])
        comd <- intersect(rownames(data_i), ids_i)
        if (length(comd)==0) return(data.frame(predict_p=numeric(0),
                                               predict_result=factor(character()),
                                               real_label=factor(character())))
        expdata <- data_i[comd, , drop=FALSE]
        prob_df <- predict(model, newdata = expdata, type = "prob")
        p <- .get_pos_prob_df(prob_df)
        pr <- factor(ifelse(p > cuti, "positive", "negative"))
        out <- data.frame(predict_p=p, predict_result=pr)
        out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2])==1, "positive", "negative"))
        out
      })
      
      all_result <- c(list(Train = train_result), all_result)
      nice_names <- c(
        "DatasetA(Train)",
        paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
               "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
      )
      names(all_result) <- nice_names[seq_along(all_result)]
      result_FS  <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
      result_recall <- sapply(all_result, function(x){
        tp <- sum(as.character(x$predict_result)=="positive" & as.character(x$real_label)=="positive")
        fn <- sum(as.character(x$predict_result)=="negative" & as.character(x$real_label)=="positive")
        if ((tp+fn)==0) 0 else tp/(tp+fn)
      })
      is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
      key <- if (is_auto) {
        paste0("SVM-CV:", fold, " fold (kernel: linear) (cutoff:auto(",
               formatC(th_auto, format="f", digits=4), ",", auto_th_method, "))")
      } else {
        paste0("SVM-CV:", fold, " fold (kernel: linear) (cutoff:",
               formatC(cuti, format="f", digits=4), ")")
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
  }
  
  # ---- 2) svmRadial ----
  {
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(ifelse(as_binary01(train_labels) == 1, "pos", "neg"),
                                levels = c("neg","pos"))
    ctrl <- caret::trainControl(method = "cv", number = fold, classProbs = TRUE)
    model <- caret::train(labels ~ . - labels, data = train_expd,
                          method = "svmRadial", trControl = ctrl)
    
    cutoff_now <- cutoff
    val_exp <- as.data.frame(test_exp[[1]])
    val_lab <- labels_list[[2]]; ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    th_auto <- NA_real_; auto_th_method <- "youden"
    if (length(comd_val) >= 2) {
      val_x <- val_exp[comd_val, , drop = FALSE]
      prob_val_df <- predict(model, newdata = val_x, type = "prob")
      prob_val <- .get_pos_prob_df(prob_val_df)
      if (!is.null(prob_val) && length(unique(prob_val[is.finite(prob_val)])) >= 2) {
        y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
        auto_th_method <- decide_threshold_method(prob_val, y_val)
        th_auto <- choose_threshold(prob_val, y_val, method = auto_th_method)
        if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) cutoff_now <- unique(c(cutoff_now, th_auto))
      }
    }
    
    for (cuti in cutoff_now) {
      prob_tr_df <- predict(model, newdata = as.data.frame(train_exp), type = "prob")
      prob_tr <- .get_pos_prob_df(prob_tr_df)
      train_result <- data.frame(
        predict_p = prob_tr,
        predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
      )
      train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
      
      all_result <- lapply(seq_along(test_exp), function(i){
        data_i <- as.data.frame(test_exp[[i]])
        label_df <- labels_list[[i+1]]
        ids_i <- as.character(label_df[[1]])
        comd <- intersect(rownames(data_i), ids_i)
        if (length(comd)==0) return(data.frame(predict_p=numeric(0),
                                               predict_result=factor(character()),
                                               real_label=factor(character())))
        expdata <- data_i[comd, , drop=FALSE]
        prob_df <- predict(model, newdata = expdata, type = "prob")
        p <- .get_pos_prob_df(prob_df)
        pr <- factor(ifelse(p > cuti, "positive", "negative"))
        out <- data.frame(predict_p=p, predict_result=pr)
        out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2])==1, "positive", "negative"))
        out
      })
      
      all_result <- c(list(Train = train_result), all_result)
      nice_names <- c(
        "DatasetA(Train)",
        paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
               "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
      )
      names(all_result) <- nice_names[seq_along(all_result)]
      result_FS  <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
      result_recall <- sapply(all_result, function(x){
        tp <- sum(as.character(x$predict_result)=="positive" & as.character(x$real_label)=="positive")
        fn <- sum(as.character(x$predict_result)=="negative" & as.character(x$real_label)=="positive")
        if ((tp+fn)==0) 0 else tp/(tp+fn)
      })
      is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
      key <- if (is_auto) {
        paste0("SVM-CV:", fold, " fold (kernel: radial) (cutoff:auto(",
               formatC(th_auto, format="f", digits=4), ",", auto_th_method, "))")
      } else {
        paste0("SVM-CV:", fold, " fold (kernel: radial) (cutoff:",
               formatC(cuti, format="f", digits=4), ")")
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
  }
  
  # ---- 3) svmPoly ----
  {
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(ifelse(as_binary01(train_labels) == 1, "pos", "neg"),
                                levels = c("neg","pos"))
    ctrl <- caret::trainControl(method = "cv", number = fold, classProbs = TRUE)
    model <- caret::train(labels ~ . - labels, data = train_expd,
                          method = "svmPoly", trControl = ctrl)
    
    cutoff_now <- cutoff
    val_exp <- as.data.frame(test_exp[[1]])
    val_lab <- labels_list[[2]]; ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    th_auto <- NA_real_; auto_th_method <- "youden"
    if (length(comd_val) >= 2) {
      val_x <- val_exp[comd_val, , drop = FALSE]
      prob_val_df <- predict(model, newdata = val_x, type = "prob")
      prob_val <- .get_pos_prob_df(prob_val_df)
      if (!is.null(prob_val) && length(unique(prob_val[is.finite(prob_val)])) >= 2) {
        y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
        auto_th_method <- decide_threshold_method(prob_val, y_val)
        th_auto <- choose_threshold(prob_val, y_val, method = auto_th_method)
        if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) cutoff_now <- unique(c(cutoff_now, th_auto))
      }
    }
    
    for (cuti in cutoff_now) {
      prob_tr_df <- predict(model, newdata = as.data.frame(train_exp), type = "prob")
      prob_tr <- .get_pos_prob_df(prob_tr_df)
      train_result <- data.frame(
        predict_p = prob_tr,
        predict_result = factor(ifelse(prob_tr > cuti, "positive", "negative"))
      )
      train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
      
      all_result <- lapply(seq_along(test_exp), function(i){
        data_i <- as.data.frame(test_exp[[i]])
        label_df <- labels_list[[i+1]]
        ids_i <- as.character(label_df[[1]])
        comd <- intersect(rownames(data_i), ids_i)
        if (length(comd)==0) return(data.frame(predict_p=numeric(0),
                                               predict_result=factor(character()),
                                               real_label=factor(character())))
        expdata <- data_i[comd, , drop=FALSE]
        prob_df <- predict(model, newdata = expdata, type = "prob")
        p <- .get_pos_prob_df(prob_df)
        pr <- factor(ifelse(p > cuti, "positive", "negative"))
        out <- data.frame(predict_p=p, predict_result=pr)
        out$real_label <- factor(ifelse(as_binary01(label_df[match(comd, ids_i), 2])==1, "positive", "negative"))
        out
      })
      
      all_result <- c(list(Train = train_result), all_result)
      nice_names <- c(
        "DatasetA(Train)",
        paste0("Dataset", LETTERS[1 + seq_along(test_exp)],
               "(", c("Val", rep("Test", length(test_exp) - 1)), ")")
      )
      names(all_result) <- nice_names[seq_along(all_result)]
      result_FS  <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) == as.character(x$real_label)))
      result_recall <- sapply(all_result, function(x){
        tp <- sum(as.character(x$predict_result)=="positive" & as.character(x$real_label)=="positive")
        fn <- sum(as.character(x$predict_result)=="negative" & as.character(x$real_label)=="positive")
        if ((tp+fn)==0) 0 else tp/(tp+fn)
      })
      is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
      key <- if (is_auto) {
        paste0("SVM-CV:", fold, " fold (kernel: polynomial) (cutoff:auto(",
               formatC(th_auto, format="f", digits=4), ",", auto_th_method, "))")
      } else {
        paste0("SVM-CV:", fold, " fold (kernel: polynomial) (cutoff:",
               formatC(cuti, format="f", digits=4), ")")
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
  }
  
  invisible(list(
    acc        = acc_local,
    recall     = recall_local,
    fs         = fs_local,
    summary    = summary_local,
    importance = importance_local   # ← 本地返回的重要性
  ))
}





#' LASSO + Optimal SVM Interface (caret CV, threshold on B)
#'
#' @description
#' Train optimal SVM models on the LASSO-selected feature subset (`lassogene`).
#' Uses caret cross-validation (svmLinear / svmRadial / svmPoly), chooses an
#' automatic threshold on validation set B (F1 or Youden), merges it with fixed
#' thresholds, evaluates on datasets A/B/C/D..., and records results.
#'
#' @param train_exp data.frame; training expression matrix (samples x features).
#' @param train_labels vector; training labels (0/1 or convertible to 0/1).
#' @param test_exp list; external datasets (first = validation set B, then C/D/...).
#' @param labels_list list; label data.frames (A=1, B=2, ...; first col = ID, second = label).
#' @param fold integer; number of CV folds for caret (tagged in result keys, default 10).
#' @param lassogene character; feature names selected by LASSO.
#' @param cutoff numeric; fixed thresholds merged with the auto threshold from B
#'   (default `c(0.25, 0.5, 0.75)`).
#' @param collector environment or NULL; optional results collector. If supplied,
#'   `.fh_collect(collector, key, acc=..., recall=..., fs=..., summary=...)` is called.
##' @param auto_th_method character; \code{"youden"} / \code{"f1"} / \code{"auto"} (default \code{"auto"}).
#'   With \code{"auto"}, \code{decide_threshold_method()} is used based on class imbalance and PR-vs-ROC gain.
#' @param auto_imbalance_thresh numeric; imbalance threshold for \code{decide_threshold_method()} (default \code{0.35}).
#' @param auto_pr_vs_roc_gate numeric; PR-vs-ROC advantage ratio gate for \code{decide_threshold_method()} (default \code{0.5}).
#' 
#' 
#' @return
#' Collector environment (if provided) is updated, and the function also returns
#' an **invisible list** of per-key results:
#' `list(acc = <named numeric>, recall = <named numeric>, fs = <named numeric>, summary = <named lists>)`.
#'
#' @examples
#' \dontrun{
#' # Minimal example (with a collector):
#' collector <- new.env(parent = emptyenv())
#' out <- fh_lasso_svm_best(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]], exp_list[[4]]),  # B / C / D
#'   labels_list  = labels_list,
#'   fold         = 10,
#'   lassogene    = lassogene,
#'   collector    = collector
#' )
#' }
#'
#' @importFrom e1071 svm
#' @importFrom caret train trainControl
#' @importFrom dplyr mutate across
#' @export
fh_lasso_svm_best <- function(train_exp, 
                              train_labels, 
                              test_exp, 
                              labels_list,
                              fold = 10, 
                              lassogene,
                              cutoff = c(0.25, 0.5, 0.75),
                              auto_th_method = "auto",
                              auto_imbalance_thresh = 0.35,
                              auto_pr_vs_roc_gate   = 0.5,
                              collector = collector) {
  # local accumulators (returned invisibly)
  acc_local     <- list()
  recall_local  <- list()
  fs_local      <- list()
  summary_local <- list()
  
  # helpers
  .get_pos_prob_df <- function(df) {
    if (is.null(df) || !ncol(df)) return(NULL)
    cn <- colnames(df)
    pick <- which(cn %in% c("1","positive","POS","Yes","TRUE","True","pos"))
    if (length(pick) == 0) pick <- ncol(df)
    as.numeric(df[[pick[1]]])
  }
  .make_pred_df_prob <- function(prob_valec, cuti, y_true = NULL) {
    out <- data.frame(
      predict_p = as.numeric(prob_valec),
      predict_result = factor(ifelse(prob_valec > cuti, "positive", "negative"))
    )
    if (!is.null(y_true)) {
      out$real_label <- factor(ifelse(as_binary01(y_true) == 1, "positive", "negative"))
    }
    out
  }
  .impute_median_df <- function(df) {
    df <- as.data.frame(df)
    dplyr::mutate(df, dplyr::across(tidyselect::where(is.numeric),
                                    ~ ifelse(is.na(.), stats::median(., na.rm = TRUE), .)))
  }
  
  # --- common objects
  X_train_all <- .impute_median_df(train_exp)
  y_train     <- train_labels
  
  make_names <- function(n_test) {
    nm <- c("DatasetA(Train)",
            paste0("Dataset", LETTERS[1:n_test + 1],
                   "(", c("Val", rep("Test", n_test - 1)), ")"))
    nm
  }
  
  # inner worker for one caret SVM method on LASSO subset
  .run_caret_svm <- function(method_name, kernel_label) {
    train_expd <- X_train_all[, lassogene, drop = FALSE]
    train_expd <- .impute_median_df(train_expd)
    train_expd$labels <- factor(ifelse(as_binary01(y_train)==1, "pos","neg"),
                                levels = c("neg","pos"))
    ctrl  <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE)
    model <- caret::train(labels ~ . - labels, data = train_expd,
                          method = method_name, trControl = ctrl)
    
    cutoff_now <- cutoff
    # ---- auto threshold on B
    val_exp <- .impute_median_df(as.data.frame(test_exp[[1]])[, lassogene, drop = FALSE])
    val_lab <- labels_list[[2]]
    ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    
    th_auto <- NA_real_
    method_used <- "youden"  # default placeholder
    if (length(comd_val) >= 2) {
      prob_val_df <- predict(model, newdata = val_exp[comd_val, , drop = FALSE], type = "prob")
      prob_val <- .get_pos_prob_df(prob_val_df)
      if (!is.null(prob_val) && length(unique(prob_val[is.finite(prob_val)])) >= 2) {
        y_val <- as_binary01(val_lab[match(comd_val, ids_val), 2])
        method_used <- decide_threshold_method(prob_val, y_val)  # 自动决定
        th_auto <- choose_threshold(prob_val, y_val, method = method_used)
        if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 && length(th_auto) == 1L && th_auto > 0.01 && th_auto < 0.99) {
          cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
          message(sprintf("[lasso+SVM+CV best %s] Auto threshold on B = %.4f (method=%s)",
                          kernel_label, th_auto, method_used))
        } else {
          message(sprintf("[lasso+SVM+CV best %s] Auto threshold skipped; keep fixed thresholds.", kernel_label))
        }
      } else {
        message(sprintf("[lasso+SVM+CV best %s] B has insufficient probability variation; keep fixed thresholds.", kernel_label))
      }
    } else {
      message(sprintf("[lasso+SVM+CV best %s] Validation overlap < 2; keep fixed thresholds.", kernel_label))
    }
    
    # ---- evaluate each cutoff
    for (cuti in cutoff_now) {
      prob_tr_df <- predict(model, newdata = X_train_all[, lassogene, drop = FALSE], type = "prob")
      prob_tr <- .get_pos_prob_df(prob_tr_df)
      train_result <- .make_pred_df_prob(prob_tr, cuti, y_true = y_train)
      
      all_result <- lapply(seq_along(test_exp), function(i){
        Xi <- .impute_median_df(as.data.frame(test_exp[[i]])[, lassogene, drop = FALSE])
        Li <- labels_list[[i + 1]]
        ids_i <- as.character(Li[[1]])
        comd  <- intersect(rownames(Xi), ids_i)
        if (length(comd) == 0) {
          return(data.frame(predict_p = numeric(0),
                            predict_result = factor(character()),
                            real_label = factor(character())))
        }
        prob_df <- predict(model, newdata = Xi[comd, , drop = FALSE], type = "prob")
        p <- .get_pos_prob_df(prob_df)
        .make_pred_df_prob(p, cuti, y_true = Li[match(comd, ids_i), 2])
      })
      
      all_result <- c(list(Train = train_result), all_result)
      names(all_result) <- make_names(length(test_exp))[seq_along(all_result)]
      
      result_FS <- sapply(all_result, function(x) f1_score(x$predict_result, x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result) ==
                                                          as.character(x$real_label)))
      result_recall <- sapply(all_result, function(x){
        tp <- sum(as.character(x$predict_result)=="positive" & as.character(x$real_label)=="positive")
        fn <- sum(as.character(x$predict_result)=="negative" & as.character(x$real_label)=="positive")
        if ((tp+fn)==0) 0 else tp/(tp+fn)
      })
      
      is_auto <- is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99 && (abs(cuti - th_auto) < .Machine$double.eps^0.5)
      key <- if (is_auto) {
        paste0("SVM-CV:", fold, " fold + Lasso-CV:", fold,
               " fold (kernel: ", kernel_label, ") (cutoff:auto(",
               formatC(th_auto, format="f", digits=4), ",", method_used, "))")
      } else {
        paste0("SVM-CV:", fold, " fold + Lasso-CV:", fold,
               " fold (kernel: ", kernel_label, ") (cutoff:",
               formatC(cuti, format="f", digits=4), ")")
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
  }
  
  # ========== run three kernels via caret ==========
  .run_caret_svm("svmLinear",   "linear")
  .run_caret_svm("svmPoly",     "polynomial")
  .run_caret_svm("svmRadial",   "radial")
  
  
  invisible(list(
    acc     = acc_local,
    recall  = recall_local,
    fs      = fs_local,
    summary = summary_local
  ))
}

