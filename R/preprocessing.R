#' Preprocess Expression Data (log transform + reference scaling) and Split Train/Test
#'
#' @description
#' This function wraps the preprocessing pipeline:
#' \enumerate{
#'   \item Compute positive label proportions (overall and per dataset).  
#'   \item Auto-detect and apply \code{log2(x+1)} transform for each dataset in \code{exp_list}.  
#'   \item Standardize all datasets using reference set mean/variance  
#'         (supporting column alignment by \code{"exact"} or \code{"intersect"}).  
#'   \item Split data into training and test sets according to \code{train_exp_name}.  
#' }
#'
#' @param exp_list list of matrices/data.frames. Each element is a dataset (rows = samples, cols = genes).  
#' @param labels_list list of data.frames (or a single data.frame). Each element has two columns: \code{sample_id}, \code{label}.  
#' @param exp_files character vector. File names corresponding to \code{exp_list}, used for naming and reporting.  
#' @param train_exp_name character. The file name of the training set (must exist in \code{exp_files}).  
#' @param ref_index integer. Index of the reference dataset (e.g., 1 = first dataset).  
#' @param include_ref logical. Whether to z-score the reference set itself (default = TRUE).  
#' @param align character. Column alignment strategy: \code{"exact"} (strict match) or \code{"intersect"} (intersect with reference).  
#'
#' @return list A list containing:  
#' \itemize{
#'   \item \code{pos_rate_summary} — positive label proportion summary (overall + per dataset)  
#'   \item \code{log_report} — dataset-wise log transform decisions and summaries  
#'   \item \code{scale_report} — reference scaling report  
#'   \item \code{exp_list} — processed (log + scaled) expression matrices  
#'   \item \code{train_exp} — training set expression matrix  
#'   \item \code{test_exp} — named list of test sets (named by \code{exp_files})  
#' }
#'
#' @export
fh_preprocess_data <- function(exp_list,
                               labels_list,
                               exp_files,
                               train_exp_name,
                               ref_index   = 1,
                               include_ref = TRUE,
                               align       = c("exact", "intersect")) {
  if (is.data.frame(labels_list)) labels_list <- list(labels_list)
  if (!length(exp_files)) stop("exp_files must not be empty.")
  if (!train_exp_name %in% exp_files) {
    stop("train_exp_name not found in exp_files: ", train_exp_name,
         "\navailable: ", paste(exp_files, collapse = ", "))
  }
  if (is.null(names(exp_list))) names(exp_list) <- exp_files
  
  # 1) positive rate
  if (exists("check_positive_rate")) {
    pos_rate_summary <- tryCatch(
      check_positive_rate(labels_list),
      error = function(e) .fh_pos_rate_fallback(labels_list)
    )
  } else {
    pos_rate_summary <- .fh_pos_rate_fallback(labels_list)
  }
  
  # 2) log transform
  log_report <- NULL
  if (exists("log_transform_list")) {
    res <- log_transform_list(exp_list)
    if (is.list(res) && !is.null(res$data)) {
      exp_list  <- res$data
      log_report <- res$report
    } else {
      stop("log_transform_list() must return a list with $data and $report.")
    }
  } else {
    tmp <- lapply(exp_list, .fh_log_auto_one)
    exp_list  <- lapply(tmp, `[[`, "mat")
    log_report <- do.call(rbind, lapply(tmp, `[[`, "report"))
  }
  if (!is.null(log_report) && nrow(log_report) == length(exp_list)) {
    rownames(log_report) <- names(exp_list)
  }
  
  # 3) reference scaling
  scale_report <- NULL
  if (exists("scale_by_ref_list")) {
    exp_scaled <- scale_by_ref_list(exp_list, ref_index = ref_index,
                                    include_ref = include_ref, align = align)
    scale_report <- attr(exp_scaled, "report")
    exp_list <- exp_scaled
  } else {
    out <- .fh_scale_by_ref_list_fallback(exp_list, ref_index, include_ref, align)
    exp_list <- out$data
    scale_report <- out$report
  }
  if (!is.null(scale_report) && is.data.frame(scale_report) &&
      nrow(scale_report) == length(exp_list)) {
    rownames(scale_report) <- names(exp_list)
  }
  
  # 4) split train/test
  train_idx <- which(exp_files == train_exp_name)[1]
  train_exp <- exp_list[[train_idx]]
  test_idx  <- setdiff(seq_along(exp_list), train_idx)
  test_exp  <- if (length(test_idx)) exp_list[test_idx] else list()
  if (length(test_exp)) {
    if (is.null(names(test_exp))) names(test_exp) <- exp_files[test_idx]
  }
  
  res1 <- list(
    pos_rate_summary = pos_rate_summary,
    log_report       = log_report,
    scale_report     = scale_report,
    exp_list         = exp_list,
    train_exp        = train_exp,
    test_exp         = test_exp
  )
  message("Data has been processed successfully.")
  invisible(res1)
}

#' @title LASSO/Elastic Net-based Feature Prescreening with Flexible Strictness
#'
#' @description
#' Perform feature pre-screening using LASSO or Elastic Net regression on the training set.
#' Features are selected based on non-zero coefficients at a lambda chosen from cross-validation.
#' Users can control strictness (lenient/balanced/conservative/custom) to tune the number of selected features.
#' Optionally, test matrices are subset to the selected features.
#' A cross-validation curve (CV error vs. log(lambda)) can be plotted and optionally saved.
#'
#' @details
#' - `x_train` should be a gene expression matrix (rows = samples, cols = genes).
#' - `y_train` can be factor, numeric, or character vectors, but must be binary.
#' - Missing values are imputed with the column median if `impute_median = TRUE`.
#' - `cv.glmnet` determines `lambda.min` (most predictive) and `lambda.1se` (simpler, more regularized).
#'   - **lenient** -> use `lambda.min` (more features).
#'   - **conservative** -> use `lambda.1se` (fewer features, more stable).
#'   - **balanced** -> midpoint on the log scale between `lambda.min` and `lambda.1se`.
#'   - **custom_fraction** -> interpolate between `lambda.min` and `lambda.1se` (controlled by `fraction`).
#' - Final coefficients are extracted at the chosen lambda. `min_keep` and `max_keep` bound the feature count.
#'
#' @param x_train matrix/data.frame; training expression data (rows = samples, cols = genes).
#' @param y_train factor/numeric/character; binary labels (two-level factor or convertible).
#' @param x_tests optional named list of matrices to be subset to the selected genes.
#' @param alpha numeric; Elastic Net mixing parameter (`1` = LASSO, `0 < alpha < 1` = Elastic Net).
#' @param nfolds integer; number of folds for `cv.glmnet` cross-validation (default `10`).
#' @param nlambda integer; number of candidate lambda values (default `100`).
#' @param intercept logical; whether the final glmnet fit includes an intercept (default `FALSE`).
#' @param impute_median logical; whether to impute missing numeric values with medians (default `TRUE`).
#' @param selection character; lambda selection strategy: one of
#'   `"lenient"`, `"balanced"`, `"conservative"`, `"custom_fraction"`.
#' @param fraction numeric in `、code{[0, 1]}`; used only if `selection = "custom_fraction"`.
#'   `0` = `lambda.min`, `1` = `lambda.1se`; default `0.5` = midpoint.
#' @param min_keep integer; minimum number of features to retain (default `0`).
#' @param max_keep integer; maximum number of features to retain (default `Inf`).
#' @param plot logical; if `TRUE`, draw the CV curve (default `TRUE`).
#' @param plot_file character or `NULL`; if non-`NULL`, save the plot to this path (without extension if you open a device inside).
#' @param plot_width numeric; plot width in inches when saving (default `7`).
#' @param plot_height numeric; plot height in inches when saving (default `5`).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{genes}: selected gene names.
#'   \item \code{x_train}: reduced training matrix with selected genes.
#'   \item \code{x_tests}: reduced test matrices with selected genes (if provided).
#'   \item \code{cvfit}: `cv.glmnet` object (CV curve info).
#'   \item \code{lambda_min}: optimal lambda from CV.
#'   \item \code{lambda_1se}: conservative lambda from CV.
#'   \item \code{lambda_used}: lambda actually used for selection.
#' }
#' If `plot = TRUE` and `plot_file` is non-`NULL`, a figure will also be written to disk.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(1000), nrow = 50, ncol = 20)
#' colnames(x) <- paste0("gene", 1:20)
#' y <- sample(c(0, 1), 50, replace = TRUE)
#'
#' # Lenient selection (lambda.min), with plot displayed
#' res1 <- fh_lasso_prescreen(x, y, selection = "lenient", plot = TRUE)
#'
#' # Conservative selection (lambda.1se), save plot to file
#' res2 <- fh_lasso_prescreen(
#'   x, y, selection = "conservative",
#'   plot = TRUE, plot_file = "cv_curve_lasso.pdf",
#'   plot_width = 7, plot_height = 5
#' )
#'
#' # Custom fraction (75% towards conservative)
#' res3 <- fh_lasso_prescreen(x, y, selection = "custom_fraction", fraction = 0.75)
#'
#' # Check selected genes
#' print(res1$genes)
#' }
#'
#' @seealso [glmnet::cv.glmnet], [glmnet::glmnet]
#' @export
fh_lasso_prescreen <- function(x_train,
                               y_train,
                               x_tests = NULL,
                               alpha = 1,
                               nfolds = 10,
                               nlambda = 100,
                               intercept = FALSE,
                               impute_median = TRUE,
                               selection = c("lenient","balanced","conservative","custom_fraction"),
                               fraction = 0.5,
                               min_keep = 0,
                               max_keep = Inf,
                               plot = c("none","cv","path","both"),
                               plot_file = NULL,
                               plot_width = 7,
                               plot_height = 5) {
  # -- helper: to matrix + optional median impute --
  to_matrix <- function(x) {
    if (is.data.frame(x)) {
      if (impute_median) {
        num_cols <- vapply(x, is.numeric, logical(1))
        if (any(num_cols)) {
          x[num_cols] <- lapply(x[num_cols], function(col) {
            if (anyNA(col)) col[is.na(col)] <- stats::median(col, na.rm = TRUE)
            col
          })
        }
      }
      x <- as.matrix(x)
    } else {
      if (impute_median && is.matrix(x)) {
        for (j in seq_len(ncol(x))) {
          v <- x[, j]
          if (anyNA(v)) {
            med <- stats::median(v, na.rm = TRUE)
            v[is.na(v)] <- med
            x[, j] <- v
          }
        }
      }
    }
    storage.mode(x) <- "double"
    x
  }
  
  X <- to_matrix(x_train)
  
  # y handling: keep factor if already binary; otherwise binarize then factor
  if (is.factor(y_train)) {
    if (nlevels(y_train) != 2) stop("y_train must be binary (2 levels).")
    y_glm <- y_train
  } else {
    y01 <- as_binary01(y_train)
    y_glm <- factor(y01, levels = c(0, 1))
  }
  
  # 1) CV to get lambda.min / lambda.1se
  cvfit <- glmnet::cv.glmnet(
    X, y_glm,
    family  = "binomial",
    alpha   = alpha,
    nlambda = nlambda,
    nfolds  = nfolds
  )
  lambda_min  <- cvfit$lambda.min
  lambda_1se  <- cvfit$lambda.1se
  if (!is.finite(lambda_min)) stop("cv.glmnet failed to produce lambda.min.")
  if (!is.finite(lambda_1se)) lambda_1se <- lambda_min
  
  # choose lambda according to strictness
  selection <- match.arg(selection)
  lambda_used <- switch(
    selection,
    lenient       = lambda_min,
    conservative  = lambda_1se,
    balanced      = exp(0.5 * (log(lambda_min) + log(lambda_1se))),  
    custom_fraction = {
      if (!is.finite(fraction) || fraction < 0 || fraction > 1)
        stop("fraction must be in [0,1] when selection='custom_fraction'.")
      exp((1 - fraction) * log(lambda_min) + fraction * log(lambda_1se))
    }
  )
  
  # 2) Final path fit (keep defaults except intercept = user controlled)
  fit <- glmnet::glmnet(
    X, y_glm,
    family    = "binomial",
    alpha     = alpha,
    intercept = intercept
  )
  

  plot <- match.arg(plot)
  do_plot_cv   <- plot %in% c("cv","both")
  do_plot_path <- plot %in% c("path","both")
  
  open_pdf <- function(path) {
    grDevices::pdf(path, width = plot_width, height = plot_height)
  }
  close_dev <- function() {
    try(grDevices::dev.off(), silent = TRUE)
  }
  
  if (do_plot_cv) {
    if (!is.null(plot_file)) open_pdf(paste0(plot_file, "_cv.pdf"))
    graphics::plot(cvfit)   # 只画CV曲线，不加竖线和图例
    if (!is.null(plot_file)) close_dev()
  }
  
  if (do_plot_path) {
    if (!is.null(plot_file)) open_pdf(paste0(plot_file, "_path.pdf"))
    graphics::plot(fit, xvar = "lambda")   # 只画系数路径，不加竖线和图例
    if (!is.null(plot_file)) close_dev()
  }
  # ============================
  
  # coefficients at chosen lambda
  co_use <- stats::coef(fit, s = lambda_used)
  co_use <- as.matrix(co_use)
  sel    <- rownames(co_use)[which(co_use != 0)]
  sel    <- setdiff(sel, "(Intercept)")
  
  # guard rails using top-|coef| at lambda.min if needed
  co_min <- stats::coef(fit, s = lambda_min)
  co_min <- as.matrix(co_min)
  co_min <- co_min[setdiff(rownames(co_min), "(Intercept)"), , drop = FALSE]
  ord_min <- order(abs(co_min[, 1]), decreasing = TRUE)
  top_min <- rownames(co_min)[ord_min]
  
  # ensure minimum
  if (length(sel) < min_keep) {
    need <- setdiff(top_min, sel)
    sel  <- unique(c(sel, head(need, max(0, min_keep - length(sel)))))
  }
  # cap maximum
  if (is.finite(max_keep) && length(sel) > max_keep) {
    co_use_noi <- co_use[setdiff(rownames(co_use), "(Intercept)"), , drop = FALSE]
    ord_use <- order(abs(co_use_noi[, 1]), decreasing = TRUE)
    ranked  <- rownames(co_use_noi)[ord_use]
    ranked  <- intersect(ranked, sel)
    if (length(ranked) >= max_keep) {
      sel <- ranked[seq_len(max_keep)]
    } else {
      pad <- setdiff(top_min, ranked)
      sel <- c(ranked, head(pad, max_keep - length(ranked)))
    }
  }
  
  if (!length(sel)) {
    warning("No non-zero coefficients at selected lambda; falling back to top-1 by |coef| at lambda.min.")
    sel <- head(top_min, 1)
  }
  
  keep  <- intersect(sel, colnames(X))
  X_sel <- X[, keep, drop = FALSE]
  
  tests_sel <- NULL
  if (!is.null(x_tests)) {
    tests_sel <- lapply(x_tests, function(m) {
      m <- to_matrix(m)
      cols <- intersect(keep, colnames(m))
      m[, cols, drop = FALSE]
    })
  }
  
  res <- list(
    genes       = keep,
    x_train     = X_sel,
    x_tests     = tests_sel,
    cvfit       = cvfit,
    lambda_used = lambda_used,
    lambda_min  = lambda_min,
    lambda_1se  = lambda_1se
  )
  message(
    sprintf(
      "LASSO prescreen completed: %d genes (selection=%s, lambda=%.3g).",
      length(keep), selection, lambda_used
    )
  )
  invisible(res)
}

# =========================== fallback implementations ===========================

.fh_pos_rate_fallback <- function(labels_list) {
  one <- function(df) {
    colnames(df)[1:2] <- c("sample_id","label")
    tab <- table(df$label)
    n1 <- ifelse("1" %in% names(tab), tab[["1"]], sum(df$label == 1))
    n0 <- ifelse("0" %in% names(tab), tab[["0"]], sum(df$label == 0))
    data.frame(
      n = length(df$label),
      pos = as.integer(n1),
      neg = as.integer(n0),
      pos_rate = ifelse((n1 + n0) > 0, n1 / (n1 + n0), NA_real_)
    )
  }
  per_set <- lapply(labels_list, one)
  per_set_df <- do.call(rbind, per_set)
  rownames(per_set_df) <- paste0("set", seq_along(labels_list))
  all_df <- one(do.call(rbind, labels_list))
  structure(list(per_set = per_set_df, overall = all_df), class = "fh_pos_rate")
}

.fh_log_auto_one <- function(mat) {
  m <- as.matrix(mat)
  x <- as.numeric(m)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(list(mat = m, report = data.frame(needs_log = NA, q99 = NA, frac_nonint = NA)))
  }
  q99 <- stats::quantile(x, probs = 0.99, na.rm = TRUE, names = FALSE)
  frac_nonint <- mean(abs(x - round(x)) > 0, na.rm = TRUE)
  is_count_like <- (q99 > 20) && (frac_nonint < 0.2)
  
  if (is_count_like) {
    m2 <- log2(m + 1)
    rpt <- data.frame(needs_log = TRUE, q99 = as.numeric(q99), frac_nonint = frac_nonint)
    list(mat = m2, report = rpt)
  } else {
    rpt <- data.frame(needs_log = FALSE, q99 = as.numeric(q99), frac_nonint = frac_nonint)
    list(mat = m, report = rpt)
  }
}

.fh_scale_by_ref_list_fallback <- function(exp_list, ref_index, include_ref, align = c("exact","intersect")) {
  align <- match.arg(align)
  ref <- exp_list[[ref_index]]
  if (is.null(colnames(ref))) stop("Reference dataset has no column names.")
  report <- list()
  
  proc_one <- function(x) {
    if (align == "exact") {
      if (!identical(colnames(x), colnames(ref))) {
        stop("align='exact' requires identical column names and order.")
      }
      x
    } else {
      keep <- intersect(colnames(ref), colnames(x))
      x[, keep, drop = FALSE][, colnames(ref)[colnames(ref) %in% keep], drop = FALSE]
    }
  }
  
  aligned <- lapply(exp_list, proc_one)
  
  ref_mu <- colMeans(aligned[[ref_index]], na.rm = TRUE)
  ref_sd <- apply(aligned[[ref_index]], 2, sd, na.rm = TRUE)
  ref_sd[ref_sd == 0 | is.na(ref_sd)] <- 1
  
  res <- aligned
  for (i in seq_along(res)) {
    z <- sweep(res[[i]], 2, ref_mu, FUN = "-")
    z <- sweep(z,       2, ref_sd, FUN = "/")
    if (i == ref_index && !include_ref) z <- res[[i]]
    res[[i]] <- z
    report[[i]] <- data.frame(
      idx         = i,
      nrow        = nrow(z),
      ncol        = ncol(z),
      aligned     = TRUE,
      include_ref = ifelse(i == ref_index, include_ref, NA)
    )
  }
  structure(list(data = res, report = do.call(rbind, report)),
            class = "fh_scale_by_ref")
}