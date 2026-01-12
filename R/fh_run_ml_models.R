#' Run all 27 ML models in one shot (orchestrator wrapper)
#'
#' @description
#' This function orchestrates and sequentially calls your existing 27
#' model runners (listed below), collecting their return values with tryCatch
#' and returning a named list. It **does not** modify your existing functions or
#' helpers—just a thin wrapper for “batch run + aggregate results”.
#'
#' @section Model runners included:
#' - `fh_logistic()`, `fh_lda()`, `fh_lasso_lda()`, `fh_qda()`, `fh_lasso_qda()`,  
#'   `fh_knn()`, `fh_lasso_knn()`, `fh_tree()`, `fh_lasso_tree()`, `fh_rf()`,  
#'   `fh_lasso_rf()`, `fh_xgboost()`, `fh_lasso_xgboost_default()`,  
#'   `fh_lasso_xgboost_best()`, `fh_ridge_glmnet()`, `fh_lasso()`,  
#'   `fh_elastic_net()`, `fh_svm()`, `fh_svm_best()`, `fh_lasso_svm_best()`,  
#'   `fh_gbm_default()`, `fh_gbm_cv_best()`, `fh_gbm_lasso_default()`,  
#'   `fh_gbm_lasso_best()`, `fh_stepwise_lr()`, `fh_naive_bayes()`,  
#'   `fh_naive_bayes_lasso()`
#'
#' @param train_exp matrix/data.frame  
#' Training expression matrix/data frame (rows=samples, cols=features).
#'
#' @param train_labels vector/factor  
#' Training labels (0/1 or binary factor).
#'
#' @param test_exp list  
#' List of external datasets, e.g., `list(B, C, D)`.
#'
#' @param labels_list list  
#' List of label data.frames aligned with datasets (`[[1]]` for A/Train, `[[2]]` for B/Val).
#'
#' @param lassogene character NULL  
#' Vector of LASSO-selected feature names, passed to *lasso* variants.
#'
#' @param com_genes character NULL  
#' Shared feature set for Random Forest (if required by your rf flow).
#'
#' @param all_labels any NULL  
#' If `fh_logistic()` needs `all_labels`, it will be passed through.
#'
#' @param fold integer = 10  
#' Number of folds (default 10) for functions that need it.
#'
#' @param cutoff numeric = c(0.25,0.5,0.75)  
#' Base cutoff grid to pass into models that accept it.
#'
#' @param knumber integer = c(1,3,5)  
#' Vector of K values for KNN.
#'
#' @param alpha_all numeric = seq(0.1,0.9,0.1)  
#' Alpha grid for Elastic Net.
#'
#' @param kernel_all character = c("linear","polynomial","radial")  
#' Kernels for SVM.
#'
#' @param nround integer = 100  
#' Default number of boosting rounds for XGBoost.
#'
#' @param max_depth integer = 6  
#' Maximum depth of XGBoost trees (default 6).
#'
#' @param eta numeric = 0.5  
#' Learning rate (eta) for XGBoost (default 0.5).
#'
#' @param auto_th_method character = "youden"  
#' Auto threshold method for functions exposing it.
#'
#' @param cores integer(NULL)  
#' Parallel cores. If `NULL`, uses `max(1, detectCores()-1)`. Passed to GBM CV-based functions that need it.
#'
#' @param which_models character = "all"  
#' Which model runners to execute. `"all"` runs all 27, or pass a subset.
#'
#' @return list  
#' Named list (size ≤ 27) of each model’s return value. Errors are caught and returned as `NULL` with a message.
#'
#' @examples
#' \dontrun{
#' res_all <- fh_run_ml_models(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = test_exp,
#'   labels_list  = labels_list,
#'   lassogene    = lassogene,
#'   com_genes    = com_genes,
#'   all_labels   = all_labels,
#'   fold         = 10,
#'   cutoff       = c(0.25, 0.5, 0.75),
#'   knumber      = c(1,3,5),
#'   alpha_all    = seq(0.1,0.9,0.1),
#'   kernel_all   = c("linear","polynomial","radial"),
#'   nround       = 100,
#'   max_depth    = 6,
#'   eta          = 0.5,
#'   auto_th_method = "youden",
#'   cores        = NULL,
#'   which_models = "all"
#' )
#' }
#' @export
fh_run_ml_models <- function(
    train_exp,
    train_labels,
    test_exp,
    labels_list,
    lassogene,
    com_genes,
    all_labels,
    fold       = 10,
    cutoff     = c(0.25, 0.5, 0.75),
    knumber    = c(1,3,5),
    alpha_all  = seq(0.1, 0.9, 0.1),
    kernel_all = c("linear","polynomial","radial"),
    nround     = 100,
    max_depth  = 6,
    eta        = 0.5,
    auto_th_method = "youden",
    cores      = NULL,
    which_models = "all"
){
  # fallback cores
  if (is.null(cores)) {
    cores <- 1L
    if (requireNamespace("parallel", quietly = TRUE)) {
      cores <- max(1L, tryCatch(parallel::detectCores()-1L, error = function(e) 1L))
    }
  }
  
  # all 27 calls
  calls <- list(
    fh_logistic = list(
      train_exp = train_exp,
      train_labels = train_labels,
      test_exp = test_exp,    # B/C/D
      labels_list = labels_list,
      all_labels = all_labels,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lda = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso_lda = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = fold,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_qda = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso_qda = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = fold,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_knn = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      knumber      = c(1,2,3,4,5),
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso_knn = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = 10,
      knumber      = c(1,2,3,4,5),
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_tree = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso_tree = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = 10,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_rf = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      com_genes    = com_genes,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso_rf = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = 10,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_xgboost = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      nround       = 100,
      max_depth    = 6,
      eta          = 0.5,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso_xgboost_default = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso_xgboost_best = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      fold         = 10,
      lassogene    = lassogene,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_ridge_glmnet = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      fold         = 10,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_lasso = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      fold         = 10,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_elastic_net = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      fold         = 10,
      alpha_all    = seq(0, 1, 0.1),
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_svm = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      kernel_all   = c("linear","polynomial","radial"),
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_svm_best = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      fold         = 10,
      lassogene    = lassogene,
      collector = collector
    ),
    
    fh_lasso_svm_best = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      fold         = 10,
      lassogene    = lassogene,
      collector = collector
    ),
    
    fh_gbm_default = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_gbm_cv_best = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      fold         = 10,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_gbm_lasso_default = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = 10,
      auto_th_method = "auto",
      collector = collector
    ),
    
    fh_gbm_lasso_best = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = 10,
      collector = collector
    ),
    
    fh_stepwise_lr = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      collector = collector
    ),
    
    fh_naive_bayes = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      collector = collector
    ),
    
    fh_naive_bayes_lasso = list(
      train_exp    = train_exp,
      train_labels = train_labels,
      test_exp     = test_exp,
      labels_list  = labels_list,
      lassogene    = lassogene,
      fold         = 10,
      collector = collector
    )
  )
  
  all_names <- names(calls)
  if (identical(which_models,"all")) {
    run_names <- all_names
  } else {
    run_names <- intersect(which_models, all_names)
    if (!length(run_names)) stop("`which_models` contains no valid model names.")
  }
  
  .safe_call <- function(fname, args){
    if (!exists(fname, mode="function")) return(NULL)
    fun <- get(fname, mode="function")
    tryCatch(do.call(fun, args), error=function(e) NULL)
  }
  
  results <- setNames(vector("list", length(run_names)), run_names)
  for (nm in run_names) {
    results[[nm]] <- .safe_call(nm, calls[[nm]])
  }
  
  return(message = "CONGRATULATIONS! All or chosen machine learning models have finished running successfully.")
}












