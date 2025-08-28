#' @title Initialize Result Containers
#'
#' @description
#' Creates and returns a list of empty lists to store results from subsequent
#' machine learning or deep learning analyses (e.g., predictions, accuracy,
#' recall, F-score, feature importance).
#'
#' @return A list with the following named components (each initialized as an empty list):  
#' \itemize{
#'   \item \code{summary}: prediction scores or class labels
#'   \item \code{acc}: accuracy results
#'   \item \code{recall}: recall results
#'   \item \code{FS}: F-score results
#'   \item \code{importance}: feature/gene importance results
#' }
#'
#' @examples
#' # Initialize empty containers
#' res <- fh_init_results()
#' str(res)
#' # $summary
#' # $acc
#' # $recall
#' # $FS
#' # $importance
#'
#' @export
fh_init_results <- function() {
  list(
    summary    = list(),
    acc        = list(),
    recall     = list(),
    FS         = list(),
    importance = list()
  )
}

