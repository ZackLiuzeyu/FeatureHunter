#' K-Nearest Neighbors (KNN)
#'
#' @description
#'
#' Iterates through a set of k values, trains and predicts for each;
#' selects optimal threshold on validation set B (F1 / Youden / auto),
#' merges with fixed cutoffs, and evaluates F1 / ACC / REC on A/B/C/D datasets.
#' Results are stored in the \code{collector} environment.
#'
#' @param train_exp Training expression matrix (rows = samples, cols = features).  
#' @param train_labels Binary vector or factor (0/1).  
#' @param test_exp A list: first element = validation set B, others = C/D/…  
#' @param labels_list A list of label data frames (A=1, B=2, …).  
#' @param knumber Integer vector of k values (default \code{c(1,2,3,4,5)}).  
#' @param cutoff Fixed cutoff values (default \code{c(0.25,0.5,0.75)}).  
#' @param auto_th_method Threshold method: "auto" (default), "f1", or "youden".  
#' @param auto_imbalance_thresh imbalance threshold for \code{decide_threshold_method()} (default 0.35).  
#' @param auto_pr_vs_roc_gate PR-vs-ROC advantage gate (default 0.5).  
#' @param collector Environment to store results.  
#'
#' @return invisible list of results (acc, recall, fs, summary) also written into \code{collector}.  
#'
#' @examples
#' \dontrun{
#' collector <- new.env()
#' fh_knn(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]]),
#'   labels_list  = labels_list,
#'   knumber      = c(1,3,5),
#'   auto_th_method = "auto",
#'   collector    = collector
#' )
#' }
#'
#' @importFrom class knn
#' @export
fh_knn <- function(train_exp, train_labels, test_exp, labels_list,
                   knumber = c(1,2,3,4,5),
                   cutoff = c(0.25, 0.5, 0.75),
                   auto_th_method = "auto",
                   auto_imbalance_thresh = 0.35,
                   auto_pr_vs_roc_gate   = 0.5,
                   collector = collector) {
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  train_imp <- as.data.frame(train_exp)
  train_imp[is.na(train_imp)] <- colMeans(train_imp, na.rm=TRUE)
  y_fac <- factor(ifelse(train_labels==1, "1", "0"), levels=c("0","1"))
  
  for (k_i in knumber) {
    val_exp <- as.data.frame(test_exp[[1]])
    val_lab <- labels_list[[2]]
    ids_val <- as.character(val_lab[[1]])
    comd_val <- intersect(rownames(val_exp), ids_val)
    th_auto <- NA; method_used <- auto_th_method
    cutoff_now <- cutoff
    if (length(comd_val) >= 2) {
      val_x <- val_exp[comd_val,,drop=FALSE]
      val_x[is.na(val_x)] <- colMeans(val_x, na.rm=TRUE)
      pred_val <- class::knn(train=as.matrix(train_imp), test=as.matrix(val_x), cl=y_fac, k=k_i, prob=TRUE)
      p_hat <- as.numeric(attr(pred_val,"prob"))
      prob_val <- ifelse(as.character(pred_val)=="1", p_hat, 1-p_hat)
      y_val <- as_binary01(val_lab[match(comd_val, ids_val),2])
      method_used <- auto_th_method
      if (identical(auto_th_method, "auto")) {
      method_used <- decide_threshold_method(
        probs = prob_val, y_true = y_val,
        imbalance_thresh = auto_imbalance_thresh,
        pr_vs_roc_gate   = auto_pr_vs_roc_gate
      )
    }
    
    th_auto <- choose_threshold(prob_val, y_val, method = method_used)
      if (is.finite(th_auto)   && th_auto > 0.01 && th_auto < 0.99) {
        cutoff_now <- sort(unique(c(cutoff_now, th_auto)))
        message(sprintf("[KNN k=%d] Auto threshold on B = %.4f (method=%s)", k_i, th_auto, method_used))
      }else {
        message("[KNN k=%d] Auto threshold failed; keep fixed thresholds.")
      }
    }
    
    for (cuti in cutoff_now) {
      pred_tr <- class::knn(train=as.matrix(train_imp), test=as.matrix(train_imp), cl=y_fac, k=k_i, prob=TRUE)
      p_hat_tr <- as.numeric(attr(pred_tr,"prob"))
      prob_tr <- ifelse(as.character(pred_tr)=="1", p_hat_tr, 1-p_hat_tr)
      train_result <- data.frame(predict_p=prob_tr,
                                 predict_result=factor(ifelse(prob_tr>cuti,"positive","negative")),
                                 real_label=factor(ifelse(train_labels==1,"positive","negative")))
      
      all_result <- lapply(seq_along(test_exp), function(i){
        data_i <- as.data.frame(test_exp[[i]])
        label_df <- labels_list[[i+1]]
        ids_i <- as.character(label_df[[1]])
        comd  <- intersect(rownames(data_i), ids_i)
        if (!length(comd)) return(data.frame(predict_p=numeric(0),predict_result=factor(character()),real_label=factor(character())))
        expdata <- data_i[comd,,drop=FALSE]; expdata[is.na(expdata)] <- colMeans(expdata, na.rm=TRUE)
        labelsdata <- label_df[match(comd,ids_i),2]
        pred_i <- class::knn(train=as.matrix(train_imp), test=as.matrix(expdata), cl=y_fac, k=k_i, prob=TRUE)
        p_hat_i <- as.numeric(attr(pred_i,"prob"))
        prob_i <- ifelse(as.character(pred_i)=="1", p_hat_i, 1-p_hat_i)
        pr <- factor(ifelse(prob_i>cuti,"positive","negative"))
        data.frame(predict_p=prob_i, predict_result=pr,
                   real_label=factor(ifelse(as_binary01(labelsdata)==1,"positive","negative")))
      })
      all_result <- c(list(Train=train_result), all_result)
      names(all_result) <- c("DatasetA(Train)", paste0("Dataset", LETTERS[1:length(test_exp)], "(", c("Val", rep("Test",length(test_exp)-1)), ")"))
      
      result_FS <- sapply(all_result, function(x) f1_score(x$predict_result,x$real_label))
      result_acc <- sapply(all_result, function(x) mean(as.character(x$predict_result)==as.character(x$real_label)))
      result_recall <- sapply(all_result,function(x){tp<-sum(as.character(x$predict_result)=="positive"&as.character(x$real_label)=="positive"); fn<-sum(as.character(x$predict_result)=="negative"&as.character(x$real_label)=="positive"); if((tp+fn)==0)0 else tp/(tp+fn)})
      is_auto <- is.finite(th_auto) && abs(cuti-th_auto)<.Machine$double.eps^.5
      key <- if(is_auto) paste0("KNN (k=",k_i,") (cutoff:auto(",formatC(th_auto,format="f",digits=4),",",method_used,"))")
      else paste0("KNN (k=",k_i,") (cutoff:",formatC(cuti,format="f",digits=4),")")
      if(is.environment(collector)) .fh_collect(collector,key,result_acc,result_recall,result_FS,all_result)
      acc_local[[key]]<-result_acc; recall_local[[key]]<-result_recall; fs_local[[key]]<-result_FS; summary_local[[key]]<-all_result
    }
  }
}


#' @title LASSO + KNN Interface 
#'
#' @description
#' Uses only the LASSO-selected features (\code{lassogene}) to train KNN classifiers.  
#' Thresholds optimized on validation set B (F1/Youden/auto) and merged with fixed cutoffs.  
#' Metrics (F1, ACC, REC) computed on A/B/C/D… and results stored in \code{collector}.
#'
#'
#' @param train_exp Training expression matrix.  
#' @param train_labels Labels (0/1).  
#' @param test_exp A list: first=B, others=C/D…  
#' @param labels_list Label data frames (A=1, B=2…).  
#' @param lassogene Selected feature names from LASSO.  
#' @param fold Integer used in result key only.  
#' @param knumber Vector of k values (default \code{c(1,2,3,4,5)}).  
#' @param cutoff Fixed cutoff values (default \code{c(0.25,0.5,0.75)}).  
#' @param auto_th_method "auto"/"f1"/"youden".  
#' @param auto_imbalance_thresh imbalance threshold (default 0.35).  
#' @param auto_pr_vs_roc_gate PR-vs-ROC gate (default 0.5).  
#' @param collector Environment for results.  
#'
#' @return invisible list of results (acc, recall, fs, summary).  
#'
#' @examples
#' \dontrun{
#' collector <- new.env()
#' fh_lasso_knn(train_exp,train_labels,list(exp_list[[2]],exp_list[[3]]),labels_list,
#'              lassogene=lassogene,fold=10,knumber=c(1,3,5),collector=collector)
#' }
#'
#' @importFrom class knn
#' @export
fh_lasso_knn <- function(train_exp, train_labels, test_exp, labels_list,
                         lassogene, fold,
                         knumber=c(1,2,3,4,5),
                         cutoff=c(0.25,0.5,0.75),
                         auto_th_method="auto",
                         auto_imbalance_thresh=0.35,
                         auto_pr_vs_roc_gate=0.5,
                         collector=collector){
  acc_local <- list(); recall_local <- list(); fs_local <- list(); summary_local <- list()
  train_imp<-as.data.frame(train_exp)[,lassogene,drop=FALSE]; train_imp[is.na(train_imp)]<-colMeans(train_imp,na.rm=TRUE)
  y_fac<-factor(ifelse(train_labels==1,"1","0"),levels=c("0","1"))
  
  for(k_i in knumber){
    val_exp<-as.data.frame(test_exp[[1]])[,lassogene,drop=FALSE]
    val_lab<-labels_list[[2]]; ids_val<-as.character(val_lab[[1]]); comd_val<-intersect(rownames(val_exp),ids_val)
    th_auto<-NA; method_used<-auto_th_method; cutoff_now<-cutoff
    if(length(comd_val)>=2){
      val_x<-val_exp[comd_val,,drop=FALSE]; val_x[is.na(val_x)]<-colMeans(val_x,na.rm=TRUE)
      pred_val<-class::knn(train=as.matrix(train_imp),test=as.matrix(val_x),cl=y_fac,k=k_i,prob=TRUE)
      p_hat<-as.numeric(attr(pred_val,"prob")); prob_val<-ifelse(as.character(pred_val)=="1",p_hat,1-p_hat)
      y_val<-as_binary01(val_lab[match(comd_val,ids_val),2])
      method_used <- auto_th_method
      if (identical(auto_th_method, "auto")) {
        method_used <- decide_threshold_method(
          probs = prob_val, y_true = y_val,
          imbalance_thresh = auto_imbalance_thresh,
          pr_vs_roc_gate   = auto_pr_vs_roc_gate
        )
      }
      
      th_auto <- choose_threshold(prob_val, y_val, method = method_used)
      if(is.finite(th_auto)  && th_auto > 0.01 && th_auto < 0.99){
        cutoff_now<-sort(unique(c(cutoff_now,th_auto)));
        message(sprintf("[LASSO+KNN k=%d] Auto threshold on B = %.4f (method=%s)",k_i,th_auto,method_used))}
    }else {
      message("[KNN k=%d] Auto threshold failed; keep fixed thresholds.")
    }
    
    for(cuti in cutoff_now){
      pred_tr<-class::knn(train=as.matrix(train_imp),test=as.matrix(train_imp),cl=y_fac,k=k_i,prob=TRUE)
      p_hat_tr<-as.numeric(attr(pred_tr,"prob")); prob_tr<-ifelse(as.character(pred_tr)=="1",p_hat_tr,1-p_hat_tr)
      train_result<-data.frame(predict_p=prob_tr,predict_result=factor(ifelse(prob_tr>cuti,"positive","negative")),real_label=factor(ifelse(train_labels==1,"positive","negative")))
      
      all_result<-lapply(seq_along(test_exp),function(i){data_i<-as.data.frame(test_exp[[i]])[,lassogene,drop=FALSE]; label_df<-labels_list[[i+1]]; ids_i<-as.character(label_df[[1]]); comd<-intersect(rownames(data_i),ids_i); if(!length(comd))return(data.frame(predict_p=numeric(0),predict_result=factor(character()),real_label=factor(character()))); expdata<-data_i[comd,,drop=FALSE]; expdata[is.na(expdata)]<-colMeans(expdata,na.rm=TRUE); labelsdata<-label_df[match(comd,ids_i),2]; pred_i<-class::knn(train=as.matrix(train_imp),test=as.matrix(expdata),cl=y_fac,k=k_i,prob=TRUE); p_hat_i<-as.numeric(attr(pred_i,"prob")); prob_i<-ifelse(as.character(pred_i)=="1",p_hat_i,1-p_hat_i); pr<-factor(ifelse(prob_i>cuti,"positive","negative")); data.frame(predict_p=prob_i,predict_result=pr,real_label=factor(ifelse(as_binary01(labelsdata)==1,"positive","negative")))})
      all_result<-c(list(Train=train_result),all_result); names(all_result)<-c("DatasetA(Train)",paste0("Dataset",LETTERS[1:length(test_exp)],"(",c("Val",rep("Test",length(test_exp)-1)),")"))
      
      result_FS<-sapply(all_result,function(x)f1_score(x$predict_result,x$real_label))
      result_acc<-sapply(all_result,function(x)mean(as.character(x$predict_result)==as.character(x$real_label)))
      result_recall<-sapply(all_result,function(x){tp<-sum(as.character(x$predict_result)=="positive"&as.character(x$real_label)=="positive"); fn<-sum(as.character(x$predict_result)=="negative"&as.character(x$real_label)=="positive"); if((tp+fn)==0)0 else tp/(tp+fn)})
      is_auto<-is.finite(th_auto)&&abs(cuti-th_auto)<.Machine$double.eps^.5
      key<-if(is_auto)paste0("KNN+Lasso-CV:",fold," fold (k=",k_i,") (cutoff:auto(",formatC(th_auto,format="f",digits=4),",",method_used,"))") else paste0("KNN+Lasso-CV:",fold," fold (k=",k_i,") (cutoff:",formatC(cuti,format="f",digits=4),")")
      if(is.environment(collector)) .fh_collect(collector,key,result_acc,result_recall,result_FS,all_result)
      acc_local[[key]]<-result_acc; recall_local[[key]]<-result_recall; fs_local[[key]]<-result_FS; summary_local[[key]]<-all_result
    }
  }
}

