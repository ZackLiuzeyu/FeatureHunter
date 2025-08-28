# R/zzz.R

.onLoad <- function(libname, pkgname) {
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "  ======================================================\n",
    "   Welcome to FeatureHunter v0.1.1 !!!\n",
    "   Hint: Use fh_init_tf() if you plan to use deep learning models.\n",
    "   28 Interfaces Implemented (Grouped):\n\n",
    
    "    Deep Learning\n",
    "     1) fh_mlp                              # Multilayer Perceptron (MLP)\n\n",
    
    "    Linear Models Regression\n",
    "     2) fh_logistic                         # Logistic Regression\n",
    "     3) fh_ridge_glmnet                     # Ridge Regression (cv.glmnet)\n",
    "     4) fh_lasso_glmnet                     # LASSO (cv.glmnet, alpha=1)\n",
    "     5) fh_elastic_net                      # Elastic Net\n",
    "     6) fh_stepwise_logistic                # Stepwise Logistic Regression\n\n",
    
    "    Discriminant Analysis & KNN\n",
    "     7)  fh_lda                             # Linear Discriminant Analysis (LDA)\n",
    "     8)  fh_lasso_lda                       # LASSO Preselection + LDA\n",
    "     9)  fh_qda                             # Quadratic Discriminant Analysis (QDA)\n",
    "     10) fh_lasso_qda                       # LASSO Preselection + QDA\n",
    "     11) fh_knn                             # k-Nearest Neighbors (KNN)\n",
    "     12) fh_lasso_knn                       # LASSO Preselection + KNN\n\n",
    
    "    Tree Models\n",
    "     13) fh_tree                            # Decision Tree\n",
    "     14) fh_lasso_tree                      # LASSO Preselection + Decision Tree\n",
    "     15) fh_rf                              # Random Forest\n",
    "     16) fh_lasso_rf                        # LASSO Preselection + Random Forest\n\n",
    
    "    Boosting Series\n",
    "     17) fh_xgboost                         # XGBoost\n",
    "     18) fh_lasso_xgboost_default           # LASSO Preselection + XGBoost (default)\n",
    "     19) fh_lasso_xgboost_best              # LASSO Preselection + XGBoost (best tuned)\n",
    "     20) fh_gbm                             # Gradient Boosting Machine (GBM)\n",
    "     21) fh_gbm_best                        # GBM (CV tuned best)\n",
    "     22) fh_lasso_gbm                       # LASSO Preselection + GBM\n",
    "     23) fh_lasso_gbm_best                  # LASSO Preselection + GBM (best)\n\n",
    
    "    Support Vector Machine (SVM)\n",
    "     24) fh_svm                             # Basic SVM\n",
    "     25) fh_svm_best                        # SVM (CV tuned best)\n",
    "     26) fh_lasso_svm_best                  # LASSO Preselection + SVM (best)\n\n",
    
    "    Bayesian\n",
    "     27) fh_naive_bayes                     # Naive Bayes\n",
    "     28) fh_lasso_naive_bayes               # LASSO Preselection + Naive Bayes\n\n",
    
    "   ------------------------------------------------------\n",
    "   Utilities: Feature Selection, Thresholding, Visualization\n",
    "  ======================================================\n"
  )
}

# Register global variables to silence R CMD check NOTES for <<- assignments
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "all_result_acc",
    "all_result_recall",
    "all_result_FS",
    "all_result_summary",
    "all_result_importance"
  ))
}

#' @title Initialize TensorFlow / Keras Environment
#'
#' @description
#' Checks whether TensorFlow is available in the current Python environment.  
#' If not installed, a message is displayed.  
#' If installed, Eager Execution is enabled.
#'
#' @return No return value.  
#' Side effect: TensorFlow is initialized and set to Eager Execution mode if available.
#'
#' @examples
#' \dontrun{
#' fh_init_tf()
#' }
#'
#' @export
fh_init_tf <- function() {
  if (!reticulate::py_module_available("tensorflow")) {
    stop("TensorFlow not installed. Please run keras3::install_tensorflow() first.")
  }
  tf <- tensorflow::tf
  if (!tf$executing_eagerly()) {
    tf$compat$v1$enable_eager_execution()
    message("TensorFlow eager execution enabled.")
  } else {
    message("TensorFlow already in eager execution mode.")
  }
}

# Silence R CMD check: global variable bindings used via NSE or by design
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "collector", "prob_val",
    "allfiles", "labels_files", "labels_list", "exp_files", "exp_list",
    "com_genes", "all_labels", "x_train", "y_train", "x_tests",
    "tests_list", "heatmap_mats"
  ))
}

