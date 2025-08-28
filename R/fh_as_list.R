#' Convert Result Collector Environment to Plain Lists
#'
#' @description
#' `fh_as_lists()` takes an internal environment-based collector 
#' (created via `.fh_new_collector()` and populated by model training 
#' functions such as `fh_logistic`, `fh_rf`, `fh_svm`, etc.) and 
#' converts it into plain R lists.  
#' This allows results to be easily inspected, saved, or passed to 
#' downstream visualization/analysis functions.
#'
#' @param collector environment  
#' An environment created by `.fh_new_collector()` and passed into 
#' model functions to store intermediate and final results.
#'
#' @return list  
#' A list of four elements, each converted from the environment:
#' \itemize{
#'   \item \code{all_result_summary} – Named list of data.frames of 
#'         predictions per dataset (must contain `predict_result`, 
#'         `real_label`, and optionally `predict_p`).  
#'   \item \code{all_result_acc} – Named list of Accuracy values.  
#'   \item \code{all_result_recall} – Named list of Recall values.  
#'   \item \code{all_result_FS} – Named list of F1 scores.  
#' }
#'
#' @examples
#' \dontrun{
#' # 1. Create a new collector
#' collector <- .fh_new_collector()
#'
#' # 2. Run one or more model functions with collector
#' fh_logistic(
#'   train_exp    = train_exp,
#'   train_labels = train_labels,
#'   test_exp     = list(exp_list[[2]], exp_list[[3]]),
#'   labels_list  = labels_list,
#'   all_labels   = all_labels,
#'   collector    = collector
#' )
#'
#' # 3. Convert collector to plain lists
#' results <- fh_as_lists(collector)
#'
#' # 4. Inspect results
#' str(results$all_result_acc)
#' str(results$all_result_summary[[1]])
#' }
#'
#' @export
fh_as_lists <- function(collector) {
  stopifnot(is.environment(collector))
  list(
    all_result_summary = .fh_as_list_or_empty(collector$all_result_summary),
    all_result_acc     = .fh_as_list_or_empty(collector$all_result_acc),
    all_result_recall  = .fh_as_list_or_empty(collector$all_result_recall),
    all_result_FS      = .fh_as_list_or_empty(collector$all_result_FS)
  )
}
