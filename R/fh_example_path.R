#' Path to example data shipped with the package
#'
#' @param ... Character components appended under inst/extdata
#' @return A file path (character scalar). Returns "" if not found.
#' @examples
#' fh_example_path("raw_data_pncs")
#' fh_example_path("raw_data_pncs", "DatasetA.txt")
#' @export
fh_example_path <- function(...) {
  system.file("extdata", ..., package = "FeatureHunter")
}