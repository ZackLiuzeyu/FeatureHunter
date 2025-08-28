#' Decision Tree Interface
#'
#' @description
#' Train a decision tree model using \code{tree::tree}.  
#' Thresholds are selected automatically on validation set B (F1/Youden/auto),  
#' merged with fixed thresholds (\code{c(0.25, 0.5, 0.75)}),  
#' and performance metrics (F1, ACC, RECALL) are evaluated on A/B/C/D… datasets.  
#' Results are written into the \code{collector} environment if supplied.
#'
#' @param train_exp data.frame Training expression matrix (samples × features).  
#' @param train_labels vector Binary or factor labels for the training set.  
#' @param test_exp list External datasets (first = validation B, then C/D/…).  
#' @param labels_list list A list of label data.frames (A=1, B=2, …; first col=ID, second col=label).  
#' @param cutoff numeric Fixed thresholds (default \code{c(0.25, 0.5, 0.75)}).  
#' @param auto_th_method Threshold method: "auto" (default), "f1", or "youden".  
#' @param auto_imbalance_thresh imbalance threshold for auto mode (default 0.35).  
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage gate (default 0.5).  
#' @param collector environment Optional collector environment for results.  
#'
#' @return collector environment and invisible list of results (acc, recall, fs, summary).
#'
#' @examples
#' \dontrun{
#' collector <- new.env()
#' fh_tree(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]]),
#'   labels_list  = labels_list,
#'   auto_th_method = "auto",
#'   collector    = collector
#' )
#' }
#'
#' @importFrom tree tree
#' @export
fh_tree <- function(train_exp, train_labels, test_exp, labels_list,
                    cutoff = c(0.25, 0.5, 0.75),
                    auto_th_method = "auto",
                    auto_imbalance_thresh = 0.35,
                    auto_pr_vs_roc_gate   = 0.5,
                    collector = collector) {
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- factor(train_labels)
  model <- tree::tree(labels ~ . - labels, data = train_expd)
  
  # Validation set B
  cutoff <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA; method_used <- auto_th_method
  if (length(comd_val) >= 2) {
    val_x <- val_exp[comd_val,,drop=FALSE]
    prob_mat <- predict(model, newdata=val_x, type="vector")
    pos_col <- intersect(colnames(prob_mat), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col) == 0) pos_col <- tail(colnames(prob_mat),1)
    prob_val <- as.numeric(prob_mat[, pos_col[1]])
    y_val <- as_binary01(val_lab[match(comd_val,ids_val),2])
    
    if (identical(auto_th_method,"auto")) {
      method_used <- decide_threshold_method(prob_val, y_val,
                                             imbalance_thresh=auto_imbalance_thresh,
                                             pr_vs_roc_gate=auto_pr_vs_roc_gate)
    }
    th_auto <- choose_threshold(prob_val, y_val, method=method_used)
    if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
      cutoff <- sort(unique(c(cutoff, th_auto)))
      message(sprintf("[DT] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[DT] Auto threshold skipped; keep original cutoff grid.")
    }
  } else {
    message("[DT] Validation overlap < 2; keep original cutoff grid.")
  }
  
  for (cuti in cutoff) {
    # Training set
    prob_tr_mat <- predict(model, newdata=as.data.frame(train_exp), type="vector")
    pos_col_tr <- intersect(colnames(prob_tr_mat), c("1","positive","POS","Yes","TRUE"))
    if (length(pos_col_tr) == 0) pos_col_tr <- tail(colnames(prob_tr_mat),1)
    prob_tr <- as.numeric(prob_tr_mat[,pos_col_tr[1]])
    
    train_result <- data.frame(
      predict_p=prob_tr,
      predict_result=factor(ifelse(prob_tr>cuti,"positive","negative")),
      real_label=factor(ifelse(train_labels==1,"positive","negative"))
    )
    
    # External sets
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i <- as.data.frame(test_exp[[i]])
      label_df <- labels_list[[i+1]]
      ids_i <- as.character(label_df[[1]])
      comd <- intersect(rownames(data_i), ids_i)
      if (!length(comd)) return(data.frame(predict_p=numeric(0),predict_result=factor(character()),real_label=factor(character())))
      expdata <- data_i[comd,,drop=FALSE]
      labelsdata <- label_df[match(comd,ids_i),2]
      prob_mat_i <- predict(model, newdata=expdata, type="vector")
      pos_col_i <- intersect(colnames(prob_mat_i), c("1","positive","POS","Yes","TRUE"))
      if (length(pos_col_i) == 0) pos_col_i <- tail(colnames(prob_mat_i),1)
      p <- as.numeric(prob_mat_i[,pos_col_i[1]])
      pr <- factor(ifelse(p>cuti,"positive","negative"))
      data.frame(predict_p=p,predict_result=pr,
                 real_label=factor(ifelse(as_binary01(labelsdata)==1,"positive","negative")))
    })
    all_result <- c(list(Train=train_result), all_result)
    names(all_result) <- c("DatasetA(Train)", paste0("Dataset", LETTERS[1:length(test_exp)], "(", c("Val", rep("Test", length(test_exp)-1)), ")"))
    
    result_FS <- sapply(all_result, function(x) f1_score(x$predict_result,x$real_label))
    result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result)==as.character(x$real_label)))
    result_recall <- sapply(all_result, function(x){tp<-sum(as.character(x$predict_result)=="positive"&as.character(x$real_label)=="positive"); fn<-sum(as.character(x$predict_result)=="negative"&as.character(x$real_label)=="positive"); if((tp+fn)==0)0 else tp/(tp+fn)})
    
    is_auto <- is.finite(th_auto) && abs(cuti-th_auto)<.Machine$double.eps^.5
    key <- if(is_auto) paste0("DT (cutoff:auto(",formatC(th_auto,format="f",digits=4),",",method_used,"))")
    else paste0("DT (cutoff:",formatC(cuti,format="f",digits=4),")")
    
    if (is.environment(collector)) .fh_collect(collector,key,result_acc,result_recall,result_FS,all_result)
    acc_local[[key]]<-result_acc; recall_local[[key]]<-result_recall; fs_local[[key]]<-result_FS; summary_local[[key]]<-all_result
  }
  invisible(list(acc=acc_local, recall=recall_local, fs=fs_local, summary=summary_local))
}


#' LASSO + Decision Tree Interface
#'
#' @description
#' Train a decision tree (\code{tree::tree}) using only features selected by LASSO.  
#' Thresholds are selected automatically on validation set B (F1/Youden/auto),  
#' merged with fixed thresholds, and performance metrics are computed on A/B/C/D… datasets.  
#' Results are written into the \code{collector} environment if supplied.
#'
#' @param train_exp data.frame Training expression matrix.  
#' @param train_labels vector Binary labels.  
#' @param test_exp list External datasets (first = validation B, then C/D/…).  
#' @param labels_list list Label data.frames (A=1, B=2, …).  
#' @param lassogene character Vector of features selected by LASSO.  
#' @param fold integer Used in result key annotation.  
#' @param cutoff numeric Fixed thresholds (default \code{c(0.25, 0.5, 0.75)}).  
#' @param auto_th_method Threshold method: "auto" (default), "f1", or "youden".  
#' @param auto_imbalance_thresh imbalance threshold for auto mode (default 0.35).  
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage gate (default 0.5).  
#' @param collector environment Optional collector environment for results.  
#'
#' @return collector environment and invisible list of results (acc, recall, fs, summary).  
#'
#' @importFrom tree tree
#' @export
fh_lasso_tree <- function(train_exp, train_labels, test_exp, labels_list,
                          lassogene, fold,
                          cutoff = c(0.25, 0.5, 0.75),
                          auto_th_method = "auto",
                          auto_imbalance_thresh = 0.35,
                          auto_pr_vs_roc_gate   = 0.5,
                          collector = collector) {
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  
  train_expd <- as.data.frame(train_exp)[, lassogene, drop=FALSE]
  train_expd$labels <- factor(train_labels)
  model <- tree::tree(labels ~ . - labels, data=train_expd)
  
  cutoff <- cutoff
  val_exp <- as.data.frame(test_exp[[1]])[,lassogene,drop=FALSE]
  val_lab <- labels_list[[2]]
  ids_val <- as.character(val_lab[[1]])
  comd_val <- intersect(rownames(val_exp), ids_val)
  
  th_auto <- NA; method_used <- auto_th_method
  if (length(comd_val)>=2) {
    val_x <- val_exp[comd_val,,drop=FALSE]
    prob_mat <- predict(model, newdata=val_x, type="vector")
    pos_col <- intersect(colnames(prob_mat), c("1","positive","POS","Yes","TRUE"))
    if (!length(pos_col)) pos_col <- tail(colnames(prob_mat),1)
    prob_val <- as.numeric(prob_mat[,pos_col[1]])
    y_val <- as_binary01(val_lab[match(comd_val,ids_val),2])
    
    if (identical(auto_th_method,"auto")) {
      method_used <- decide_threshold_method(prob_val, y_val,
                                             imbalance_thresh=auto_imbalance_thresh,
                                             pr_vs_roc_gate=auto_pr_vs_roc_gate)
    }
    th_auto <- choose_threshold(prob_val,y_val,method=method_used)
    if (is.finite(th_auto) && th_auto > 0.01 && th_auto < 0.99) {
      cutoff <- sort(unique(c(cutoff, th_auto)))
      message(sprintf("[lasso+DT] Auto threshold on B = %.4f (method=%s)", th_auto, method_used))
    } else {
      message("[lasso+DT] Auto threshold skipped; keep original cutoff grid.")
    }
  }else {
    message("[lasso+DT] Validation overlap < 2; keep original cutoff grid.")
  }
  
  for (cuti in cutoff) {
    prob_tr_mat <- predict(model, newdata=as.data.frame(train_exp)[,lassogene,drop=FALSE], type="vector")
    pos_col_tr <- intersect(colnames(prob_tr_mat), c("1","positive","POS","Yes","TRUE"))
    if (!length(pos_col_tr)) pos_col_tr <- tail(colnames(prob_tr_mat),1)
    prob_tr <- as.numeric(prob_tr_mat[,pos_col_tr[1]])
    
    train_result <- data.frame(
      predict_p=prob_tr,
      predict_result=factor(ifelse(prob_tr>cuti,"positive","negative")),
      real_label=factor(ifelse(train_labels==1,"positive","negative"))
    )
    
    all_result <- lapply(seq_along(test_exp), function(i){
      data_i <- as.data.frame(test_exp[[i]])[,lassogene,drop=FALSE]
      label_df <- labels_list[[i+1]]
      ids_i <- as.character(label_df[[1]])
      comd <- intersect(rownames(data_i), ids_i)
      if (!length(comd)) return(data.frame(predict_p=numeric(0),predict_result=factor(character()),real_label=factor(character())))
      expdata <- data_i[comd,,drop=FALSE]
      labelsdata <- label_df[match(comd,ids_i),2]
      prob_mat_i <- predict(model, newdata=expdata, type="vector")
      pos_col_i <- intersect(colnames(prob_mat_i), c("1","positive","POS","Yes","TRUE"))
      if (!length(pos_col_i)) pos_col_i <- tail(colnames(prob_mat_i),1)
      p <- as.numeric(prob_mat_i[,pos_col_i[1]])
      pr <- factor(ifelse(p>cuti,"positive","negative"))
      data.frame(predict_p=p,predict_result=pr,
                 real_label=factor(ifelse(as_binary01(labelsdata)==1,"positive","negative")))
    })
    all_result <- c(list(Train=train_result), all_result)
    names(all_result) <- c("DatasetA(Train)", paste0("Dataset", LETTERS[1:length(test_exp)], "(", c("Val", rep("Test", length(test_exp)-1)), ")"))
    
    result_FS <- sapply(all_result,function(x)f1_score(x$predict_result,x$real_label))
    result_acc <- sapply(all_result,function(x)mean(as.character(x$predict_result)==as.character(x$real_label)))
    result_recall <- sapply(all_result,function(x){tp<-sum(as.character(x$predict_result)=="positive"&as.character(x$real_label)=="positive"); fn<-sum(as.character(x$predict_result)=="negative"&as.character(x$real_label)=="positive"); if((tp+fn)==0)0 else tp/(tp+fn)})
    
    is_auto <- is.finite(th_auto) && abs(cuti-th_auto)<.Machine$double.eps^.5
    key <- if(is_auto) paste0("DT+Lasso-CV:",fold," fold (cutoff:auto(",formatC(th_auto,format="f",digits=4),",",method_used,"))")
    else paste0("DT+Lasso-CV:",fold," fold (cutoff:",formatC(cuti,format="f",digits=4),")")
    
    if (is.environment(collector)) .fh_collect(collector,key,result_acc,result_recall,result_FS,all_result)
    acc_local[[key]]<-result_acc; recall_local[[key]]<-result_recall; fs_local[[key]]<-result_FS; summary_local[[key]]<-all_result
  }
  invisible(list(acc=acc_local, recall=recall_local, fs=fs_local, summary=summary_local))
}

