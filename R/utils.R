############### Utils ======

#' Check if a matrix looks count-like
#' @description
#' Heuristically determine whether a numeric matrix/data.frame is **count-like**:
#' values are (almost) integers, mean is reasonably large, and the maximum is high.
#' Useful to decide whether a \code{log2(x + 1)} transform might be appropriate.
#' \code{log2(x + 1)} 
#' @param mat matrix|data.frame 
#' Numeric matrix/data.frame (rows = samples, cols = features).
#' @param integer_prop_threshold numeric 
#' Minimum proportion of entries that must be (near) integers to qualify as count-like.
#' Default: \code{0.995}.
#' **** \code{0.995}
#' @param mean_threshold numeric 
#' Lower bound on the **global mean**; mean must be \emph{strictly} greater than this value.
#' Default: \code{1}.
#' @param max_threshold numeric 
#' Lower bound on the **global maximum** (at least this large). Default: \code{100}.
#' @param tol numeric 
#' Tolerance when testing "near integer" via \code{|x - round(x)| <= tol}. Default: \code{1e-8}.
#' \code{|x - round(x)| <= tol} \code{1e-8}
#' @return logical 
#' \code{TRUE} if \code{mat} is considered count-like, otherwise \code{FALSE}.
#' \code{TRUE} \code{FALSE}
#' @examples
#' x <- matrix(c(0,1,2,3, 10,50,100,200), nrow = 2)
#' is_count_like(x) # likely TRUE
#' y <- matrix(rnorm(100), 10)
#' is_count_like(y) # likely FALSE
#' @export
is_count_like <- function(mat,
                          integer_prop_threshold = 0.995,
                          mean_threshold = 1,
                          max_threshold = 100,
                          tol = 1e-8) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  if (!length(x)) return(FALSE)
  
  intish <- mean(abs(x - round(x)) <= tol)  
  mu     <- mean(x)                         
  mx     <- max(x)                          
  
  (intish >= integer_prop_threshold) &&
    (mu > mean_threshold) &&
    (mx >= max_threshold)
}

#' Check if data looks "log-scaled"
#' @description
#' log 
#' This function heuristically determines whether the input values resemble
#' log-transformed data. Log-transformed values are usually within a limited
#' numeric range (often < 20) and contain a substantial proportion of
#' non-integer values.
#' @param mat 
#' A matrix or data frame containing numeric values to be tested.
#' @param q_hi_prob (0.99) 
#' 99% 
#' Quantile probability (default 0.99). The high quantile (e.g., 99th percentile) 
#' is used to approximate the upper bound of the data.
#' @param q_hi_threshold (20) 
#' log 
#' Threshold for the high quantile (default 20). If the computed quantile 
#' does not exceed this threshold, the data are considered consistent with a 
#' log-like scale.
#' @param frac_nonint_threshold (0.2) 
#' log 
#' Fraction of non-integer values required (default 0.2). A higher fraction 
#' indicates that the data are less likely to be raw counts and more likely 
#' to be log-transformed.
#' @return TRUE FALSE"log" 
#' A logical value, TRUE if the data resemble a log-like scale, FALSE otherwise.
#' @export
is_logged_like <- function(mat,
                           q_hi_prob = 0.99,
                           q_hi_threshold = 20,
                           frac_nonint_threshold = 0.2) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  if (!length(x)) return(FALSE)
  
  q_hi <- as.numeric(stats::quantile(x, probs = q_hi_prob, names = FALSE))
  frac_nonint <- mean(abs(x - round(x)) > 0)
  (q_hi <= q_hi_threshold) && (frac_nonint >= frac_nonint_threshold)
}

#' log / Decide whether a matrix needs log transformation
#' @description
#' It heuristically decides whether the input matrix (or data frame) 
#' should be log-transformed, based on whether the data resemble raw counts 
#' or already log-transformed values.
#' @param mat 
#' A matrix or data frame containing numeric values to be tested.
#' @param force `"auto"``"always"``"never"` `"auto"` 
#' `"auto"` 
#' `"always"` log 
#' `"never"` log 
#' Mode for forcing decision, one of `"auto"`, `"always"`, `"never"` (default `"auto"`): 
#' `"auto"`: automatically decide based on heuristics 
#' `"always"`: always assume log transformation is needed 
#' `"never"`: never apply log transformation
#' @param integer_prop_threshold (0.995) 
#' Proportion threshold for integer-like values (default 0.995). 
#' If the proportion of near-integers exceeds this, data are likely raw counts.
#' @param mean_threshold (1) 
#' Mean threshold (default 1). The overall mean must exceed this to qualify as counts.
#' @param max_threshold (100) 
#' Maximum value threshold (default 100). Large maximum values support the count-like hypothesis.
#' @param tol (1e-8) 
#' Tolerance (default 1e-8). Used to decide whether a value is close enough to an integer.
#' @param q_hi_prob (0.99) 
#' "log" 
#' Quantile probability (default 0.99). Used to estimate the upper tail of the distribution for log-likeness.
#' @param q_hi_threshold (20) 
#' log 
#' Threshold for the high quantile (default 20). If below this, data are consistent with log scale.
#' @param frac_nonint_threshold (0.2) 
#' log 
#' Fraction of non-integer values required (default 0.2). 
#' Ensures data are not purely counts and thus may be log-transformed.
#' @return 
#' `need_log` log 
#' `is_count` 
#' `is_logged` log 
#' Returns a list with three logical elements: 
#' `need_log`: whether log transformation is needed 
#' `is_count`: whether data are count-like 
#' `is_logged`: whether data are already log-like
#' @export
decide_log_for_matrix <- function(mat,
                                  force = c("auto", "always", "never"),
                                  integer_prop_threshold = 0.995,
                                  mean_threshold = 1,
                                  max_threshold = 100,
                                  tol = 1e-8,
                                  q_hi_prob = 0.99,
                                  q_hi_threshold = 20,
                                  frac_nonint_threshold = 0.2) {
  force <- match.arg(force)
  if (force == "always") return(list(need_log = TRUE,
                                     is_count = NA, is_logged = NA))
  if (force == "never")  return(list(need_log = FALSE,
                                     is_count = NA, is_logged = NA))
  
  is_count  <- is_count_like(mat,
                             integer_prop_threshold,
                             mean_threshold,
                             max_threshold,
                             tol)
  is_logged <- is_logged_like(mat,
                              q_hi_prob,
                              q_hi_threshold,
                              frac_nonint_threshold)
  list(need_log = is_count && !is_logged,
       is_count = is_count,
       is_logged = is_logged)
}

#' log2(x+1) / Batch decide & apply log2(x+1)
#' @description
#' `exp_list` /
#' `log2(x + 1)` 
#' For each matrix/data frame in `exp_list`, this function heuristically decides
#' whether a log transformation is needed and, if so, applies `log2(x + 1)`.
#' It also returns a per-dataset decision/summary report.
#' @param exp_list 
#' A list of matrices or data frames. Each element is treated as an expression-like matrix.
#' @param forces `length(exp_list)`
#' `"auto"``"always"``"never"` `"auto"` 
#' Optional character vector of length `length(exp_list)`, with values
#' `"auto"`, `"always"`, or `"never"`. Defaults to `"auto"` for all datasets.
#' `"auto"` log 
#' `"auto"`: decide via heuristics 
#' `"always"` `log2(x+1)` 
#' `"always"`: unconditionally apply `log2(x+1)` 
#' `"never"` 
#' `"never"`: keep original scale unconditionally
#' @param integer_prop_threshold 0.995 
#' Proportion threshold of near-integer values for count-likeness (default 0.995).
#' @param mean_threshold 1 
#' Minimum mean used for count-likeness (default 1).
#' @param max_threshold 100 
#' Maximum value threshold supporting count-likeness (default 100).
#' @param tol 1e-8 
#' Numeric tolerance for near-integer checks (default 1e-8).
#' @param q_hi_prob "log" 0.99 
#' High-quantile probability for log-likeness (default 0.99).
#' @param q_hi_threshold "log" 20 
#' Threshold for the high quantile when assessing log-likeness (default 20).
#' @param frac_nonint_threshold "log" 0.2 
#' Fraction of non-integer values required for log-likeness (default 0.2).
#' @param verbose TRUE 
#' Whether to print progress messages for each dataset (default TRUE).
#' @details
#' / \code{decide_log_for_matrix()}
#' \code{is_count_like()} \code{is_logged_like()} 
#' "log" \code{log2(x+1)} 
#' \code{forces} `"always"` `"never"` 
#' The decision scheme delegates to \code{decide_log_for_matrix()}, which uses
#' \code{is_count_like()} and \code{is_logged_like()} heuristics. A dataset is
#' transformed if it is count-like and not already log-like. The \code{forces}
#' argument overrides automatic decisions.
#' @return 
#' A list with two elements: 
#' \item{data}{ \code{exp_list} \code{log2(x+1)} 
#' A list of transformed (or original) matrices, same length as \code{exp_list}.}
#' \item{report}{\code{data.frame} 
#' A \code{data.frame} summarizing decisions with columns: 
#' \itemize{
#' \item \code{dataset}: / dataset name 
#' \item \code{rows}, \code{cols}: / dimensions 
#' \item \code{mean}, \code{max}: / mean and max over finite entries 
#' \item \code{integer_prop}: \code{tol} / proportion of near-integers 
#' \item \code{is_count_like}: / count-like flag 
#' \item \code{is_logged_like}: log / log-like flag 
#' \item \code{need_log2p1}: \code{log2(x+1)} / whether transformation applied
#' }}
#' @seealso
#' \code{\link{decide_log_for_matrix}}, \code{\link{is_count_like}}, \code{\link{is_logged_like}}
#' @examples
#' \dontrun{
#' lst <- list(A = matrix(c(0,1,3,10), 2),
#' B = matrix(rnorm(6, 5, 1), 2))
#' res <- log_transform_list(lst)
#' res$report
#' }
#' @export
log_transform_list <- function(exp_list,
                               forces = NULL,  
                               integer_prop_threshold = 0.995,
                               mean_threshold = 1,
                               max_threshold = 100,
                               tol = 1e-8,
                               q_hi_prob = 0.99,
                               q_hi_threshold = 20,
                               frac_nonint_threshold = 0.2,
                               verbose = TRUE) {
  stopifnot(is.list(exp_list), length(exp_list) >= 1)
  
  n <- length(exp_list)
  if (is.null(forces)) forces <- rep("auto", n)
  stopifnot(length(forces) == n)
  
  nm <- names(exp_list)
  if (is.null(nm) || any(nchar(nm) == 0)) nm <- paste0("set", seq_len(n))
  
  out_list <- vector("list", n)
  report <- data.frame(
    dataset = nm,
    rows = NA_integer_,
    cols = NA_integer_,
    mean = NA_real_,
    max  = NA_real_,
    integer_prop = NA_real_,
    is_count_like = NA,
    is_logged_like = NA,
    need_log2p1 = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n)) {
    M <- as.matrix(exp_list[[i]])
    x <- as.numeric(M); x <- x[is.finite(x)]
    
    
    report$rows[i] <- nrow(M)
    report$cols[i] <- ncol(M)
    report$mean[i] <- ifelse(length(x), mean(x), NA_real_)
    report$max[i]  <- ifelse(length(x), max(x),  NA_real_)
    report$integer_prop[i] <- ifelse(length(x), mean(abs(x - round(x)) <= tol), NA_real_)
    

    dec <- decide_log_for_matrix(
      M, force = forces[i],
      integer_prop_threshold = integer_prop_threshold,
      mean_threshold = mean_threshold,
      max_threshold = max_threshold,
      tol = tol,
      q_hi_prob = q_hi_prob,
      q_hi_threshold = q_hi_threshold,
      frac_nonint_threshold = frac_nonint_threshold
    )
    
    report$is_count_like[i]  <- dec$is_count
    report$is_logged_like[i] <- dec$is_logged
    report$need_log2p1[i]    <- dec$need_log
    

    if (isTRUE(dec$need_log)) {
      out <- log2(M + 1)
      if (verbose) message(sprintf("[%s] applied log2(x+1)", nm[i]))
    } else {
      out <- M
      if (verbose) message(sprintf("[%s] kept original scale", nm[i]))
    }
    rownames(out) <- rownames(M); colnames(out) <- colnames(M)
    out_list[[i]] <- out
  }
  
  names(out_list) <- nm
  list(data = out_list, report = report)
}
#' / Standardize datasets by a reference with reporting
#' @description
#' `exp_list` /
#' \code{(x - ref_mean) / ref_sd} `"exact"` `"intersect"` 
#' For each matrix/data frame in `exp_list`, perform column-wise standardization
#' using the reference dataset's column means and standard deviations:
#' \code{(x - ref_mean) / ref_sd}. Supports two column alignment strategies
#' (`"exact"` and `"intersect"`). Returns scaled datasets and a summary report.
#' @param exp_list 
#' A list of matrices or data frames, where columns typically represent genes/features.
#' @param ref_index `exp_list` 1 
#' Integer index specifying which dataset in `exp_list` serves as the reference (default 1).
#' @param include_ref / `FALSE` 
#' Logical; whether to include the reference dataset itself in the output (it will
#' be standardized by its own mean/SD). Default is `FALSE`.
#' @param align `"exact"` `"intersect"` 
#' Alignment mode, one of `"exact"` or `"intersect"`: 
#' `"exact"` 
#' `"exact"`: requires all datasets to have identical column names and order as the reference; errors otherwise. 
#' `"intersect"` 
#' `"intersect"`: keeps only the intersection of columns across datasets, ordered as in the reference.
#' @return 
#' A list with: 
#' `data`: `include_ref` 
#' Scaled list of matrices (reference included if `include_ref=TRUE`). 
#' `report`: reference/test 
#' A data.frame report with index, name, role (reference/test), dimensions, scaling flag, and note. 
#' `ref_index` 
#' Attribute `ref_index`: index of the reference dataset. 
#' `ref_mean` 
#' Attribute `ref_mean`: column means of the reference dataset. 
#' `ref_sd` 1 
#' Attribute `ref_sd`: column standard deviations of the reference dataset (0 or non-finite replaced by 1).
#' @details
#' 
#' 0 1 
#' `na.rm = TRUE` 
#' `data` 
#' @seealso \code{\link{log_transform_list}}, \code{\link{decide_log_for_matrix}}
#' @examples
#' \dontrun{
#' set.seed(123)
#' A <- matrix(rnorm(12, 10, 3), 4, 3, dimnames = list(NULL, c("g1","g2","g3")))
#' B <- matrix(rnorm(12, 20, 5), 4, 3, dimnames = list(NULL, c("g1","g2","g3")))
#' res <- scale_by_ref_list(list(A=A, B=B), ref_index=1, include_ref=TRUE, align="exact")
#' res$report
#' }
#' @export
scale_by_ref_list <- function(exp_list, ref_index = 1,
                              include_ref = FALSE,
                              align = c("exact", "intersect")) {
  align <- match.arg(align)
  stopifnot(length(exp_list) >= 1,
            ref_index >= 1, ref_index <= length(exp_list))
  
  # coerce to numeric matrix
  to_num_mat <- function(x) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) {
      storage.mode(x) <- "double"
      return(x)
    }
    stop("Each exp_list[[i]] must be a matrix or data.frame.")
  }
  mats <- lapply(exp_list, to_num_mat)
  
  ref_cols <- colnames(mats[[ref_index]])
  if (is.null(ref_cols) || any(!nzchar(ref_cols)))
    stop("Reference set is missing column names (feature/gene names).")
  
  reorder_or_check <- function(M) {
    cols <- colnames(M)
    if (align == "exact") {
      if (!identical(cols, ref_cols))
        stop(
          "Column names/order do not match the reference set. ",
          "Please ensure all datasets have exactly the same column names in the same order; ",
          "or set align='intersect' to take the intersection."
        )
      M
    } else {
      keep <- intersect(ref_cols, cols)
      if (length(keep) == 0)
        stop("No common columns with the reference set.")
      
      M2 <- M[, keep, drop = FALSE]
      M2 <- M2[, match(keep, colnames(M2)), drop = FALSE]  # keep original order
      M2 <- M2[, match(keep, ref_cols), drop = FALSE]      # reorder to reference
      colnames(M2) <- ref_cols[ref_cols %in% keep]
      M2
    }
  }
  mats <- lapply(mats, reorder_or_check)
  
  # compute reference stats
  ref_mat  <- mats[[ref_index]]
  ref_mean <- colMeans(ref_mat, na.rm = TRUE)
  ref_sd   <- apply(ref_mat, 2, sd)
  ref_sd[!is.finite(ref_sd) | ref_sd == 0] <- 1  # avoid zeros / non-finite
  
  scale_by_ref <- function(M) {
    sweep(sweep(M, 2, ref_mean, "-"), 2, ref_sd, "/")
  }
  
  out <- mats
  idx_all  <- seq_along(mats)
  test_idx <- setdiff(idx_all, ref_index)
  
  for (i in test_idx) out[[i]] <- scale_by_ref(mats[[i]])
  if (isTRUE(include_ref)) out[[ref_index]] <- scale_by_ref(mats[[ref_index]])
  
  mk_row <- function(i, role, scaled, note = "") {
    data.frame(
      index   = i,
      name    = names(exp_list)[i] %||% paste0("set", i),
      role    = role,
      rows    = nrow(out[[i]]),
      cols    = ncol(out[[i]]),
      scaled  = scaled,
      note    = note,
      stringsAsFactors = FALSE
    )
  }
  `%||%` <- function(a,b) if (is.null(a)) b else a
  
  rep_list <- list(
    mk_row(ref_index, "reference", include_ref, "ref_mean/ref_sd from this set")
  )
  for (i in test_idx) {
    rep_list[[length(rep_list)+1]] <- mk_row(i, "test", TRUE, "")
  }
  report <- do.call(rbind, rep_list)
  attr(out, "report") <- report
  attr(out, "ref_index") <- ref_index
  attr(out, "ref_mean")  <- ref_mean
  attr(out, "ref_sd")    <- ref_sd
  out
}
#' Check positive rate of label sets
#' @description
#' For each dataset in `labels_list`, this function computes the positive rate.
#' It assumes the second column contains labels, which are coerced into binary
#' (0/1). A summary table is returned with total count, positives, negatives,
#' and positive rate.
#' @param labels_list 
#' 2 `"positive"` / `"1"` / `1` 
#' A list of data frames or matrices, each with at least two columns. The 2nd
#' column is treated as labels. Values `"positive"`, `"1"`, or numeric `1`
#' are mapped to positives; all others are treated as negatives.
#' @return 
#' A data.frame with columns: 
#' `Dataset`: Set_i / dataset name (auto-generated Set_i) 
#' `Total`: / total number of samples 
#' `Positive`: / number of positives 
#' `Negative`: / number of negatives 
#' `PosRate`: / positive rate (rounded to 4 decimals)
#' @examples
#' labs1 <- data.frame(ID=1:5, Label=c("positive","negative","1","0",1))
#' labs2 <- data.frame(ID=1:4, Label=c(0,0,1,"positive"))
#' check_positive_rate(list(labs1, labs2))
#' @export
check_positive_rate <- function(labels_list) {
  res <- lapply(seq_along(labels_list), function(i) {
    lab <- labels_list[[i]][, 2]  
    lab_bin <- ifelse(lab %in% c("1","Disease","Case","Positive","Tumor","Cancer","Yes","True", 1), 1, 0) 
    total <- length(lab_bin)
    pos <- sum(lab_bin == 1)
    neg <- sum(lab_bin == 0)
    pos_rate <- pos / total
    
    data.frame(
      Dataset = paste0("Set_", i),
      Total   = total,
      Positive = pos,
      Negative = neg,
      PosRate = round(pos_rate, 4)
    )
  })
  
  do.call(rbind, res)
}


#' Compute F1 score for binary or multiclass
#' @description
#' 
#' F1 Handles both binary and multiclass settings: returns the F1 for the
#' designated positive class in binary tasks, and the macro-averaged F1 in
#' multiclass tasks.
#' @param predicted 
#' Predicted labels (character/factor/numeric; coerced to character then factored).
#' @param expected  
#' Ground-truth labels (character/factor/numeric; coerced to character then factored).
#' @param positive.class `"Disease"` `"positive"` `"1"` 
#' Character name of the positive class for binary tasks (default `"positive"`).
#' If absent from the label set, falls back to `"1"`, else to the last level.
#' @return 
#' F1 
#' F1 
#' A numeric scalar: the F1 for the positive class (binary) or macro-averaged F1 (multiclass).
#' @details
#' 
#' /PrecisionRecall
#' The function explicitly controls factor levels, puts the positive class last,
#' guards against division-by-zero with \code{pmax}, and macro-averages across classes for multiclass tasks.
#' @examples
#' # "positive"
#' pred <- c("positive","negative","positive","negative","positive")
#' true <- c("positive","negative","negative","negative","positive")
#' f1_score(pred, true)
#' # "1"
#' pred <- c(1,0,1,0,1); true <- c(1,0,0,0,1)
#' f1_score(pred, true, positive.class = "positive")
#' # F1
#' pred <- c("cat","dog","cat","bird","dog")
#' true <- c("cat","dog","bird","bird","dog")
#' f1_score(pred, true)
#' @export
f1_score <- function(predicted, expected, positive.class = "positive") {
  predicted <- as.character(predicted)
  expected  <- as.character(expected)
  
  labs <- unique(c(expected, predicted))
  if (!positive.class %in% labs) {
    if ("1" %in% labs) {
      positive.class <- "1"
    } else {
      positive.class <- tail(labs, 1)
    }
  }
  
  other <- setdiff(labs, positive.class)
  lvl <- c(other, positive.class)
  predicted <- factor(predicted, levels = lvl)
  expected  <- factor(expected,  levels = lvl)
  
  cm <- table(expected, predicted)
  
  precision <- diag(cm) / pmax(colSums(cm), 1)
  recall    <- diag(cm) / pmax(rowSums(cm), 1)
  
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  if (nlevels(expected) == 2) {
    return(unname(f1[length(f1)]))
  } else {
    return(mean(f1, na.rm = TRUE))
  }
}
#' / Check if class imbalance handling is needed
#' @description
#' This function checks the proportion of the minority class in a label vector
#' and decides whether imbalance handling is required, based on a threshold.
#' @param y 
#' A vector of class labels (character, factor, or numeric).
#' @param thresh 0.35 
#' Numeric threshold (default 0.35). If the minority class proportion is below
#' this value, imbalance handling is considered needed.
#' @return 
#' `TRUE` 
#' `FALSE` 
#' A logical scalar: 
#' `TRUE` if imbalance handling is needed, 
#' `FALSE` otherwise.
#' @details
#' \code{minority_prop()} 
#' This function calls \code{minority_prop()} internally to compute the
#' proportion of the minority class.
#' @examples
#' y <- c(rep("A", 90), rep("B", 10))
#' need_imbalance_handling(y) # TRUE
#' need_imbalance_handling(y, 0.05) # FALSE
#' @export
need_imbalance_handling <- function(y, thresh = 0.35) { minority_prop(y) < thresh }


#' / Compute balanced class weights for binary classification
#' @description
#' Computes class weights for binary classification such that each class (0/1)
#' contributes approximately equally to the overall loss, regardless of class imbalance.
#' @param y 0 1 
#' A vector of binary class labels, with values 0 and 1 (numeric or coercible).
#' @return 
#' A list with two elements: 
#' `0`: 0 / weight for class 0 
#' `1`: 1 / weight for class 1 
#' @details
#' \deqn{w_c = N / (2 \times n_c)} 
#' \eqn{N} \eqn{n_c} \eqn{c} 
#' Thus, rarer classes receive higher weights, balancing the contribution
#' of positives and negatives in training.
#' @examples
#' y <- c(0,0,0,1,1)
#' class_weights_balanced(y)
#' # $`0` = 5 / (2*3) = 0.8333
#' # $`1` = 5 / (2*2) = 1.25
#' @export
class_weights_balanced <- function(y) {
  y <- as.numeric(y); pos <- sum(y == 1); neg <- sum(y == 0); tot <- length(y)
  list(`0` = tot / (2 * neg), `1` = tot / (2 * pos))
}
#' / Initialize bias term from labels
#' @description
#' bias
#' Computes an initial bias term for binary classification, e.g. for logistic
#' regression or neural network output layers, such that the initial predicted
#' probability matches the observed positive rate.
#' @param y 0 1 1 
#' A vector of binary labels, expected values 0 and 1 (coerced to numeric and compared to 1).
#' @return 
#' \deqn{\log \left(\frac{p}{1-p} \right)} 
#' \eqn{p} \code{[1e-6, 1 - 1e-6]} 
#' A numeric scalar giving \eqn{\log(p / (1-p))}, where \eqn{p} is the positive rate,
#' clipped to the range \code{[1e-6, 1 - 1e-6]} for stability.
#' @details
#' 0 1 
#' 
#' The clipping avoids infinities when the positive rate is near 0 or 1.
#' Using this bias initialization can accelerate convergence in imbalanced datasets.
#' @examples
#' y <- c(0,0,1,1,1)
#' init_bias_from_labels(y)
#' @export
init_bias_from_labels <- function(y) {
  p <- mean(as.numeric(y) == 1)
  p <- max(min(p, 1 - 1e-6), 1e-6)
  log(p / (1 - p))
}
#' / Compute classification metrics from predictions
#' @description
#' (ACC) (RECALL)
#' F1 
#' Computes standard binary classification metrics: accuracy (ACC), recall, and F1 score
#' from predicted and true labels.
#' @param pred 0/1/ 
#' A vector of predicted labels, ideally 0/1 (logical or integer; coerced to integer).
#' @param true 0/1 
#' A vector of ground-truth labels (0/1; coerced to numeric).
#' @return 
#' A named numeric vector with: 
#' `F1`: F1 
#' `ACC`: (accuracy) 
#' `RECALL`: (recall) 
#' @details
#' ACC = / 
#' RECALL = TP / (TP + FN) 0 0 
#' F1 = 2TP / (2TP + FP + FN) 0 0 
#' Metrics are computed as: 
#' ACC = mean of correct predictions 
#' Recall = TP / (TP+FN), returns 0 if denominator = 0 
#' F1 = 2TP / (2TP+FP+FN), returns 0 if denominator = 0
#' @examples
#' pred <- c(1,0,1,0,1)
#' true <- c(1,0,0,0,1)
#' metrics_from_preds(pred, true)
#' # F1, ACC, RECALL
#' @export
metrics_from_preds <- function(pred, true) {
  y <- as.numeric(true); p <- as.integer(pred)
  acc <- mean(p == y)
  tp <- sum(p == 1 & y == 1); fp <- sum(p == 1 & y == 0); fn <- sum(p == 0 & y == 1)
  rec <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1  <- if ((2 * tp + fp + fn) == 0) 0 else (2 * tp) / (2 * tp + fp + fn)
  c(F1 = f1, ACC = acc, RECALL = rec)
}
#' (MLP) / Build a multilayer perceptron (MLP) model
#' @description
#' Keras MLP
#' dropout sigmoid Adam
#' ROC AUC PRC AUC 
#' Builds a multilayer perceptron (MLP) for binary classification using Keras,
#' with dense layers, dropout regularization, sigmoid output, Adam optimizer,
#' binary cross-entropy loss, and evaluation metrics (ROC AUC and PRC AUC).
#' @param input_dim 
#' Input dimension (integer or vector) specifying the input shape of the model.
#' @param dropout_rate dropout (0 to 1) 
#' Dropout rate (0 to 1), applied to hidden layers to reduce overfitting.
#' @param lr Adam 
#' Learning rate for the Adam optimizer.
#' @param use_bias_init `FALSE` 
#' Logical; whether to use a custom bias initializer for the output layer (default `FALSE`).
#' @param bias_value `use_bias_init=TRUE` 
#' Numeric value for output layer bias initialization when `use_bias_init=TRUE`.
#' @return Keras \code{fit()} 
#' A compiled Keras model object, ready for training via \code{fit()}.
#' @details
#' Dense(32, ReLU), then Dropout 
#' Dense(16, ReLU), then Dropout 
#' Dense(8, ReLU), then Dropout 
#' Dense(1, Sigmoid) 
#' Adam (\code{learning_rate = lr}) 
#' 
#' AUC (ROC) AUC (PRC) 
#' This architecture is a simple MLP suitable for binary classification tasks,
#' with optional bias initialization to adjust initial output probabilities.
#' @examples
#' \dontrun{
#' model <- build_mlp(input_dim = 100, dropout_rate = 0.5, lr = 1e-3)
#' summary(model)
#' }
#' @export
build_mlp <- function(input_dim, dropout_rate, lr, use_bias_init = FALSE, bias_value = 0) {
  model <- keras3::keras_model_sequential() |>
    keras3::layer_dense(32, activation = "relu", input_shape = input_dim) |>
    keras3::layer_dropout(rate = dropout_rate) |>
    keras3::layer_dense(16, activation = "relu") |>
    keras3::layer_dropout(rate = dropout_rate) |>
    keras3::layer_dense(8,  activation = "relu") |>
    keras3::layer_dropout(rate = dropout_rate) |>
    keras3::layer_dense(
      1, activation = "sigmoid",
      bias_initializer = if (use_bias_init) keras3::initializer_constant(bias_value) else "zeros"
    )
  model %>% keras3::compile(
    optimizer = keras3::optimizer_adam(learning_rate = lr),
    loss = "binary_crossentropy",
    metrics = c(keras3::metric_auc(name = "auc_roc"),
                keras3::metric_auc(curve = "PR", name = "prc"))
  )
  model
}
#' / Evaluate model predictions with a threshold
#' @description
#' 0/1 
#' For a trained model and a list of datasets, this function generates predictions,
#' applies a threshold to convert probabilities into class labels (0/1),
#' and computes evaluation metrics. Returns both raw predictions and metrics.
#' @param model Keras `predict()` 
#' A trained Keras model (or any model supporting `predict()`).
#' @param sets 
#' A list of datasets, where each element is itself a list with: 
#' `x`: / feature matrix 
#' `y`: 0/1/ label vector (0/1) 
#' `th`: / numeric threshold for binarization
#' @return `sets` 
#' A list (same length as `sets`), each element containing: 
#' `raw`: 
#' A data.frame with predicted probabilities, binary results, and true labels 
#' `met`: \code{metrics_from_preds()} F1ACCRECALL 
#' A named vector of metrics (F1, ACC, RECALL) from \code{metrics_from_preds()}
#' @seealso \code{\link{metrics_from_preds}}
#' @examples
#' \dontrun{
#' # model dlist 
#' res <- eval_with_threshold(model, dlist)
#' str(res[[1]]$met) # 
#' head(res[[1]]$raw) # 
#' }
#' @export
eval_with_threshold <- function(model, sets) {
  lapply(sets, function(d) {
    p  <- as.numeric(predict(model, d$x))
    pr <- as.integer(p >= d$th)
    list(raw = data.frame(predict_p = p, predict_result = pr, real_label = d$y),
         met = metrics_from_preds(pr, d$y))
  })
}
#' / Extract a label vector from an object
#' @description
#' Extracts a (binary) label vector from various input types. Priority/order: 
#' 1) If `obj` is already a vector/factor, return it; 
#' 2) If a matrix, return column `col` (default 2); 
#' 3) If a data.frame: 
#' If column `V2` exists, return it; 
#' Else drop any column identical to row names (ID-like), then auto-detect
#' a binary column (two levels) and return the first match. 
#' Errors with guidance if no suitable label column can be found.
#' @param obj 
#' An object carrying labels: a vector/factor, a matrix, or a data.frame.
#' @param col `obj` 2 
#' Integer column index to use when `obj` is a matrix (default 2).
#' @return 
#' A label vector, preserving the original atomic type where possible (numeric/character/factor).
#' @details
#' For data.frames, any column identical to row names is treated as an ID column
#' and dropped before detection. A "binary column" is determined as above. Clear
#' errors are raised if the requested matrix column is absent or no binary label
#' column can be identified.
#' @examples
#' # /
#' get_labels_vec(c(0,1,1,0))
#' # 2 
#' m <- cbind(id = 1:4, y = c(0,1,1,0))
#' get_labels_vec(m, col = 2)
#' # V2
#' df1 <- data.frame(V1 = 1:4, V2 = c("positive","negative","positive","negative"))
#' get_labels_vec(df1)
#' # ID 
#' df2 <- data.frame(ID = paste0("s",1:4), lab = c(0,1,1,0), check.names = FALSE)
#' rownames(df2) <- df2$ID
#' get_labels_vec(df2)
#' @export
get_labels_vec <- function(obj, col = 2) {
  if (is.vector(obj) || is.factor(obj)) return(obj)
  if (is.matrix(obj)) {
    if (ncol(obj) >= col) return(obj[, col])
    stop("matrix has no", col, "column")
  }
  if (is.data.frame(obj)) {
    rn <- tryCatch(rownames(obj), error = function(e) NULL)
    if ("V2" %in% names(obj)) return(obj[["V2"]])
    drop_idx <- rep(FALSE, ncol(obj))
    for (j in seq_len(ncol(obj))) {
      cj <- obj[[j]]
      if (!is.null(rn) && is.character(cj) &&
          identical(unname(as.character(cj)), unname(as.character(rn)))) {
        drop_idx[j] <- TRUE
      }
    }
    cand <- obj[ , !drop_idx, drop = FALSE]
    is_binary_col <- function(v) {
      u <- unique(v)
      if (is.factor(v)) return(length(u) == 2)
      if (is.logical(v)) return(TRUE)
      if (is.numeric(v)) return(all(u %in% c(0, 1)))
      if (is.character(v)) return(length(unique(tolower(trimws(v)))) == 2)
      FALSE
    }
    bins <- which(vapply(cand, is_binary_col, TRUE))
    if (length(bins) >= 1) return(cand[[bins[1]]])
    stop("no suitable label list vector found in data.frame (no V2 column, no binary column after dropping ID-like columns)")
  }
  stop("Unrecognized tag object types:", class(obj))
}
#' Convert labels into 0/1 binary vector
#' @description
#' / `"positive"`/`"negative"``"case"`/`"control"`
#' `"tumor"`/`"normal"` 
#' Converts various label types (numeric, character, logical, factor) into a 0/1
#' encoding. Supports multiple aliases for positive/negative classes such as
#' `"positive"`/`"negative"`, `"case"`/`"control"`, `"tumor"`/`"normal"`, etc.
#' @param y 
#' Input label vector: numeric, character, logical, or factor.
#' @param pos_alias `"1","disease","case","positive","tumor","cancer","yes","true"` 
#' Character vector of aliases for the positive class (default includes `"1","disease","case","positive","tumor","cancer","yes","true"`).
#' @param neg_alias `"0","normal","control","negative","no","false"` 
#' Character vector of aliases for the negative class (default includes `"0","normal","control","negative","no","false"`).
#' @return 0 1 
#' An integer vector containing only 0 and 1.
#' @details
#' This function ensures consistent binary encoding across different input types.
#' @examples
#' # 
#' as_binary01(c(0,1,1,0))
#' as_binary01(c(2,3,2,3)) # 
#' # 
#' as_binary01(c("case","control","case"))
#' as_binary01(c("positive","negative","positive"))
#' # 
#' as_binary01(c(TRUE, FALSE, TRUE))
#' @export
as_binary01 <- function(y,
                        pos_alias = c("1","disease","case","positive","tumor","cancer","yes","true"),
                        neg_alias = c("0","normal","control","negative","no","false")) {
  
  .clean <- function(x) {
    if (is.null(x)) return(x)
    x <- gsub("[\u00A0\u2000-\u200B\uFEFF]", "", x)   
    trimws(tolower(x))
  }
  
  if (is.factor(y)) y <- as.character(y)
  if (is.logical(y)) return(as.integer(y))
  
  if (is.numeric(y)) {
    uy <- sort(unique(y))
    if (all(uy %in% c(0,1))) return(as.integer(y))
    if (length(uy) == 2)     return(as.integer(y == max(uy)))
    stop("Data labels are not binary:", paste(uy, collapse=", "))
  }
  
  if (is.character(y)) {
    z   <- .clean(y)
    pos <- z %in% pos_alias
    neg <- z %in% neg_alias
    
    both_known <- pos | neg
    if (!all(both_known, na.rm = TRUE)) {
      u <- unique(z[!both_known & !is.na(z)])
      if (length(u) == 2) {
        is_pos <- grepl("(?:^|[^a-z])(disease|case|positive|tumor|cancer|1)(?:$|[^a-z0-9])", u, ignore.case=TRUE)
        if (sum(is_pos) == 1L) {
          map <- setNames(c(0L,1L), c(u[!is_pos], u[is_pos]))
          out <- rep(NA_integer_, length(z))
          out[!is.na(z)] <- unname(map[z[!is.na(z)]])
          return(out)
        }
      }
      bad <- unique(y[!(both_known %in% TRUE)])
      stop("Character labels cannot be binarized:", paste(bad, collapse=", "))
    }
    out <- integer(length(z)); out[pos %in% TRUE] <- 1L; out[neg %in% TRUE] <- 0L
    out[is.na(z)] <- NA_integer_
    return(out)
  }
  
  stop("Unknown label type", paste(class(y), collapse = "/"))
}
#' Compute minority class proportion
#' @description
#' minority class 
#' For a binary label vector, computes the proportion of the minority class.
#' @param y01 0/1 1 
#' A binary label vector, ideally encoded as 0/1. Logical or numeric values are coerced to integer and compared to 1.
#' @return \code{[0, 0.5]} 
#' A numeric scalar in \code{[0, 0.5]}, giving the proportion of the minority class.
#' @details
#' 1 \eqn{p1} 
#' \eqn{\min(p1, 1 - p1)} 
#' 0.5 0 
#' This function is often used as a helper for imbalance detection
#' (e.g. in \code{need_imbalance_handling()}).
#' @examples
#' y <- c(0,0,0,1,1)
#' minority_prop(y) # 0.4
#' y2 <- c(1,1,1,1,0)
#' minority_prop(y2) # 0.2
#' @seealso \code{\link{need_imbalance_handling}}, \code{\link{class_weights_balanced}}
#' @export

minority_prop <- function(y01) { p1 <- mean(as.integer(y01)==1); min(p1, 1-p1) }
#' Safely coerce an object to a matrix
#' @description
#' \code{as.matrix()} 
#' Safely converts input to a matrix: 
#' If already a matrix, return as-is; 
#' If a data.frame, convert via \code{as.matrix()}; 
#' If a vector, wrap into a single-column matrix; 
#' Otherwise, throw an error.
#' @param x 
#' Input object: matrix, data.frame, or vector.
#' @return 
#' A matrix object.
#' @examples
#' safe_matrix(matrix(1:4, 2))
#' safe_matrix(data.frame(a=1:3, b=4:6))
#' safe_matrix(1:5)
#' \dontrun{
#' safe_matrix(list(a=1)) # 
#' }
#' @export
safe_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  if (is.vector(x)) return(matrix(x, ncol = 1))
  stop("Cannot convert an object to a matrix:", class(x))
}
#' Stratified train/validation split
#' @description
#' Ensures stratified splitting of indices into training and validation sets
#' for binary classification labels, preserving the class distribution.
#' @param y01 0/1 
#' Binary label vector (0/1 or coercible to integer).
#' @param val_frac (0,1) 0.2 
#' Fraction of samples assigned to validation set (between 0 and 1), default 0.2.
#' @param seed 42 
#' Random seed for reproducibility (default 42).
#' @return 
#' A list with: 
#' `train_idx`: / indices for training set 
#' `val_idx`: / indices for validation set 
#' @details
#' Stratification ensures both classes appear in validation set, with minimum of 1 per class.
#' @examples
#' y <- c(rep(0, 90), rep(1, 10))
#' split <- stratified_split_idx(y, val_frac = 0.2, seed = 1)
#' length(split$train_idx) # ~80
#' length(split$val_idx) # ~20
#' @export
stratified_split_idx <- function(y01, val_frac = 0.2, seed = 42) {
  stopifnot(val_frac > 0, val_frac < 1)
  set.seed(seed)
  y01 <- as.integer(y01)
  idx_pos <- which(y01 == 1); idx_neg <- which(y01 == 0)
  n_val_pos <- max(1, floor(length(idx_pos) * val_frac))
  n_val_neg <- max(1, floor(length(idx_neg) * val_frac))
  val_idx <- c(sample(idx_pos, n_val_pos), sample(idx_neg, n_val_neg))
  val_idx <- sort(unique(val_idx))
  train_idx <- setdiff(seq_along(y01), val_idx)
  list(train_idx = train_idx, val_idx = val_idx)
}
#' Match labels for a matrix by sample IDs
#' @description
#' Given a matrix with sample IDs as row names, this function finds the best-matching
#' label table in `labels_list` (by maximizing ID overlap), extracts the label column,
#' aligns it with the matrix rows, and returns it as a binary vector. Errors if alignment fails.
#' @param X ID 
#' A matrix with row names representing sample IDs.
#' @param labels_list 
#' A list of label data.frames or matrices. Each element should have at least two columns: 
#' column 1: sample IDs (V1) 
#' column 2: labels (V2)
#' @return 0/1 
#' A binary (0/1) label vector, aligned to the row order of `X`.
#' @details
#' `labels_list` 
#' 
#' \code{as_binary01()} 0/1 
#' The function identifies the label table with the highest overlap of sample IDs, then checks for full alignment. Labels are converted to 0/1 using \code{as_binary01()}.
#' @examples
#' X <- matrix(rnorm(9), 3, 3)
#' rownames(X) <- c("s1","s2","s3")
#' labels_list <- list(data.frame(V1=c("s1","s2","s3"), V2=c("case","control","case")))
#' get_labels_for_matrix(X, labels_list)
#' @seealso \code{\link{as_binary01}}, \code{\link{get_labels_vec}}
#' @export
get_labels_for_matrix <- function(X, labels_list) {
  samp <- rownames(X)
  if (is.null(samp)) stop("The passed expression matrix lacks a row name (sample ID) and cannot be matched in labels_list.")
  best_j <- NA; best_hit <- -1L
  for (j in seq_along(labels_list)) {
    labj <- labels_list[[j]]
    idsj <- as.character(labj[[1]])  # V1
    hit <- sum(samp %in% idsj)
    if (hit > best_hit) { best_hit <- hit; best_j <- j }
  }
  if (is.na(best_j) || best_hit <= 0) {
    stop("No label matching the sample ID of this matrix can be found in labels_list.")
  }
  lab_best <- labels_list[[best_j]]
  ids_best <- as.character(lab_best[[1]])
  labvec   <- lab_best[[2]]
  idx <- match(samp, ids_best)
  if (anyNA(idx)) {
    stop(sprintf("Labels could not be perfectly aligned to the matrix: total sample %d, match success %d.",
                 length(samp), sum(!is.na(idx))))
  }
  as_binary01(labvec[idx])
}
#' Choose optimal classification threshold
#' @description
#' Selects an optimal classification threshold given predicted probabilities
#' and true labels. Supports: 
#' `"youden"`: maximizes Youden's index from the ROC curve (requires \pkg{pROC}); 
#' `"f1"`: maximizes F1 score across candidate thresholds. 
#' Falls back to `fallback` if computation fails.
#' @param probs \code{[0,1]} 
#' Numeric vector of predicted probabilities (range \code{[0,1]}).
#' @param y_true 0/1 
#' Vector of true binary labels (0/1; coerced to numeric).
#' @param method `"youden"` `"f1"` `"youden"` 
#' Method for threshold selection, `"youden"` or `"f1"` (default `"youden"`).
#' @param fallback 0.5 
#' Numeric fallback threshold if no valid threshold can be computed (default 0.5).
#' @return 
#' A numeric scalar giving the chosen threshold.
#' @details
#' This function is useful for imbalanced classification problems where the
#' default threshold 0.5 may not be optimal.
#' @examples
#' set.seed(1)
#' probs <- runif(100)
#' labels <- rbinom(100, 1, 0.3)
#' choose_threshold(probs, labels, method = "youden")
#' choose_threshold(probs, labels, method = "f1")
#' @seealso \code{\link{f1_score}}, \code{\link{metrics_from_preds}}
#' @export
choose_threshold <- function(probs, y_true, method = c("youden", "f1"), fallback = 0.5) {
  method <- match.arg(method, choices = c("youden","f1"))
  p <- as.numeric(probs); y <- as.numeric(y_true)
  if (!length(p) || length(unique(p[is.finite(p)])) < 2) return(as.numeric(fallback))
  if (method == "youden" && requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(response = y, predictor = p, quiet = TRUE)
    ths <- pROC::coords(roc_obj, x = "all",
                        ret = c("threshold", "sensitivity", "specificity"),
                        transpose = FALSE)
    ths <- as.data.frame(ths)
    ths <- ths[is.finite(ths$threshold), , drop = FALSE]
    if (!nrow(ths)) return(as.numeric(fallback))
    youden <- ths$sensitivity + ths$specificity - 1
    return(as.numeric(ths$threshold[which.max(youden)]))
  }
  cand <- sort(unique(p))
  if (length(cand) > 200)
    cand <- quantile(cand, probs = seq(0, 1, length.out = 200), names = FALSE)
  f1_at <- function(th) {
    pred <- as.integer(p >= th)
    tp <- sum(pred == 1 & y == 1); fp <- sum(pred == 1 & y == 0); fn <- sum(pred == 0 & y == 1)
    if ((2 * tp + fp + fn) == 0) return(0)
    (2 * tp) / (2 * tp + fp + fn)
  }
  f1s <- vapply(cand, f1_at, 0.0)
  cand[which.max(f1s)]
}
#' Decide method for threshold selection
#' @description
#' Automatically decides whether to use `"youden"` (ROC Youden index) or `"f1"`
#' (maximize F1 score) as the threshold selection method, based on data balance
#' and available metrics packages.
#' @param probs 
#' Numeric vector of predicted probabilities.
#' @param y_true 0/1 
#' Vector of true binary labels (0/1; coerced to integer).
#' @param imbalance_thresh 0.35 
#' `"f1"` 
#' Numeric threshold for class imbalance (default 0.35). If the minority class
#' proportion is below this value, `"f1"` is preferred.
#' @param pr_vs_roc_gate  0.5 
#' PR AUC ROC AUC `"youden"` 
#' Numeric cutoff (default 0.5) for comparing PR AUC gain vs ROC AUC gain.
#' If PR gain is relatively low, `"youden"` is chosen.
#' @return `"youden"` `"f1"` 
#' A character string, either `"youden"` or `"f1"`, indicating the chosen method.
#' @details
#' Decision logic: 
#' Returns `"f1"` if predicted probabilities are degenerate or classes are highly imbalanced. 
#' If both \pkg{pROC} and \pkg{PRROC} are available, compares PR AUC vs ROC AUC gains. 
#' Falls back to `"youden"` otherwise.
#' @examples
#' set.seed(123)
#' probs <- runif(100)
#' labels <- rbinom(100, 1, 0.2)
#' decide_threshold_method(probs, labels)
#' @seealso \code{\link{choose_threshold}}, \code{\link{f1_score}}
#' @export
decide_threshold_method <- function(probs, y_true, imbalance_thresh = 0.35, pr_vs_roc_gate = 0.5) {
  p <- as.numeric(probs); y <- as.integer(y_true)
  ok <- is.finite(p) & !is.na(y)
  p <- p[ok]; y <- y[ok]
  if (length(p) < 2 || length(unique(p)) < 2) return("f1")
  prev1 <- mean(y == 1); minority <- min(prev1, 1 - prev1)
  if (minority < imbalance_thresh) return("f1")
  has_roc  <- requireNamespace("pROC",  quietly = TRUE)
  has_pr   <- requireNamespace("PRROC", quietly = TRUE)
  if (has_roc && has_pr) {
    roc_auc <- as.numeric(pROC::auc(pROC::roc(response = y, predictor = p, quiet = TRUE)))
    pr_obj  <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = FALSE)
    pr_auc  <- unname(pr_obj$auc.integral)
    roc_gain <- max(roc_auc - 0.5, 0); pr_gain <- max(pr_auc - prev1, 0)
    ratio <- if (roc_gain > 0) pr_gain / roc_gain else 1
    if (ratio < pr_vs_roc_gate) return("youden")
  }
  "youden"
}

#' to01
#' @description Autogenerated stub documentation. Edit as needed.
#' @param x (autodetected)
#' @return See function body for details.
#' @export


to01 <- function(x) {
  x <- tolower(as.character(x))
  ifelse(x %in% c("positive","pos","1","yes","true"), 1L,
         ifelse(x %in% c("negative","neg","0","no","false"), 0L, NA_integer_))
}
# # Global switch for auto method; you can set it to "f1" or "youden" to force.
# if (!exists("auto_th_method")) auto_th_method <- "auto"


.fh_new_collector <- function() {
  new.env(parent = emptyenv())
}

.fh_collect <- function(collector, key, acc, recall, fs, summary) {
  stopifnot(is.environment(collector), length(key) == 1, is.character(key))
  
  if (is.null(collector$all_result_acc))     collector$all_result_acc     <- list()
  if (is.null(collector$all_result_recall))  collector$all_result_recall  <- list()
  if (is.null(collector$all_result_FS))      collector$all_result_FS      <- list()
  if (is.null(collector$all_result_summary)) collector$all_result_summary <- list()
  
  collector$all_result_acc[[key]]     <- acc
  collector$all_result_recall[[key]]  <- recall
  collector$all_result_FS[[key]]      <- fs
  collector$all_result_summary[[key]] <- summary
  
  invisible(collector)
}

.fh_as_list_or_empty <- function(x) if (is.null(x)) list() else as.list(x)


.get_pos_col <- function(prob_df) {
  pos_col <- intersect(colnames(prob_df), c("1","X1","positive","POS","Yes","TRUE"))
  if (length(pos_col) == 0) pos_col <- tail(colnames(prob_df), 1)
  pos_col[1]
}
.is_degenerate_probs <- function(p, tol = 1e-12) {
  p <- as.numeric(p)
  p <- p[is.finite(p)]
  if (length(p) == 0) return(TRUE)
  rng <- max(p) - min(p)
  (rng < tol)  
}
.predict_nb_with_fallback <- function(model, newx, cutoff) {
  prob_df <- tryCatch(predict(model, newdata = newx, type = "raw"),
                      error = function(e) NULL)
  if (!is.null(prob_df)) {
    pos_col <- .get_pos_col(prob_df)
    p <- as.numeric(prob_df[, pos_col])
    if (!.is_degenerate_probs(p)) {
      pr <- factor(ifelse(p > cutoff, "positive", "negative"))
      return(list(predict_p = p, predict_result = pr))
    }
  }
  cls <- predict(model, newdata = newx, type = "class")
  pr  <- factor(ifelse(as.character(cls) %in% c("1","X1","positive","POS","Yes","TRUE"),
                       "positive", "negative"))
  return(list(predict_p = rep(NA_real_, length(pr)), predict_result = pr))
}


.align_xy <- function(X, label_df) {
  stopifnot(ncol(label_df) >= 2)
  ids <- as.character(label_df[[1]])
  comd <- intersect(rownames(X), ids)
  if (length(comd) < 2) {
    return(list(x = matrix(numeric(0), nrow=0, ncol=ncol(X)),
                y = numeric(0), n = 0, overlap = 0))
  }
  Xi <- as.matrix(X[comd, , drop = FALSE])
  yi <- as_binary01(label_df[match(comd, ids), 2])
  list(x = Xi, y = yi, n = nrow(Xi), overlap = length(comd))
}



