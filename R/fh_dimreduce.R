#' Dimensionality reduction (PCA, ICA, PLS-DA) with gene ranking
#'
#' @description
#' Performs PCA, ICA, and supervised PLS-DA (via mixOmics), exports per-sample
#' scores and PCA variance explained, and computes a fused gene-level importance
#' by combining PCA, ICA, and PLS-DA loadings. Outputs are saved under
#' \code{file.path(getwd(), "preprocess")} with a dataset-specific prefix.
#' You can either:
#' \itemize{
#'   \item \strong{Dir mode} (preferred): provide \code{dir}, \code{train_exp_name},
#'         \code{hub_file}, and \code{positive_label}. The function will call
#'         \code{fh_load_data()} to load, align, and hub-filter data.
#'   \item \strong{Legacy mode}: provide \code{expr_path} and \code{meta_path}.
#' }
#' If the number of genes is below \code{min_genes}, the function exits early and
#' reports "no reduction needed".
#'
#' @param expr_path (Legacy) Path to expression matrix (tab-delimited).
#'   First column must be gene IDs; remaining columns are samples.
#' @param meta_path (Legacy) Path to sample metadata (tab-delimited, no header).
#'   Two columns: \code{Sample}, \code{Group}.
#' @param dir Directory containing training expression, label files, and hub file.
#'   If provided (non-NULL), the function calls \code{fh_load_data()} and ignores
#'   \code{expr_path/meta_path}.
#' @param train_exp_name Training expression file name inside \code{dir}.
#' @param hub_file Hub gene file name inside \code{dir}. Default: \code{"Common_genes.txt"}.
#' @param positive_label Positive class label for \code{fh_load_data()} to derive 0/1.
#' @param out_dir Output directory; defaults to \code{file.path(getwd(), "preprocess")}.
#'   Created if it does not exist.
#' @param max_components Maximum number of components to fit for PCA/ICA (capped
#'   internally by \code{min(ncol(X), nrow(X) - 1)}). Default: \code{30L}.
#' @param min_genes Minimum number of genes required to proceed. If the
#'   expression matrix has fewer columns than this threshold, the function skips
#'   reduction and returns early. Default: \code{200L}.
#' @param seed Random seed (L'Ecuyer-CMRG) for reproducibility. Default:
#'   \code{20240917L}.
#' @param normalize Logical; if \code{TRUE}, scale features to zero mean and unit
#'   variance before reduction. Default: \code{TRUE}.
#' @param ica_method \code{fastICA} algorithm method, e.g. \code{"C"}. Default:
#'   \code{"C"}.
#' @param pls_ncomp_default Default supervised components for PLS-DA
#'   (capped by matrix size). Default: \code{10L}.
#' @param group_levels Optional factor level order for \code{Group} (legacy mode only),
#'   e.g. \code{c("Normal","Disease")}. Default: \code{NULL}.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{out_dir}: resolved output directory
#'   \item \code{prefix}: dataset prefix (training file name in dir-mode; or basename(expr_path))
#'   \item \code{n_samples}: number of samples used
#'   \item \code{n_genes}: number of genes used
#'   \item \code{max_components}: effective component cap
#' }
#'
#' @details
#' \strong{Dir mode (via \code{fh_load_data()}):}
#' data are aligned and hub-filtered to the intersection of hub genes and
#' cross-dataset common genes. \code{Group} is taken from \code{y_train} as "0"/"1".
#'
#' \strong{Outputs written to \code{out_dir}} (prefix as above):
#' \itemize{
#'   \item \code{*_PCA_scores.csv}, \code{*_PCA_var_explained.csv}, \code{*_PCA_scree.png}
#'   \item \code{*_ICA_scores.csv}
#'   \item \code{*_PLS_scores.csv}
#'   \item \code{*_reduced_concat.csv} (with \code{Sample, Group})
#'   \item \code{*_gene_ranking_all.csv}, \code{*_gene_ranking_combo_pos.csv}
#'   \item \code{*_expr_combo_pos_scaled.tsv} (always)
#'   \item \code{*_expr_combo_pos_raw.tsv} (legacy mode only)
#'   \item \code{*_combo_pos_gene_list.txt}
#' }
#'
#' @examples
#' \dontrun{
#' # Legacy mode
#' res1 <- fh_dimreduce(
#'   expr_path = "DatasetD.txt",
#'   meta_path = "SampleD.txt",
#'   min_genes = 200,
#'   max_components = 30,
#'   pls_ncomp_default = 10,
#'   group_levels = c("Normal","Disease")
#' )
#'
#' # Dir mode (preferred)
#' res2 <- fh_dimreduce(
#'   dir = "dir_sample",
#'   train_exp_name = "DatasetD.txt",
#'   hub_file = "Common_genes.txt",
#'   positive_label = "Disease",
#'   max_components = 30,
#'   min_genes = 200
#' )
#' }
#'
#' @importFrom stats prcomp
#' @importFrom graphics plot legend par lines
#' @importFrom grDevices png dev.off
#' @importFrom utils installed.packages
#' @importFrom tools file_path_sans_ext
#' @export
fh_dimreduce <- function(
    expr_path = NULL,
    meta_path = NULL,
    dir = NULL,
    train_exp_name = NULL,
    hub_file = "Common_genes.txt",
    positive_label = NULL,
    out_dir = file.path(getwd(), "preprocess"),
    
    max_components = 30L,
    min_genes = 200L,
    seed = 20240917L,
    normalize = TRUE,
    ica_method = "C",
    pls_ncomp_default = 10L,
    group_levels = NULL
) {
  .ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  .msg_write <- function(path) message(sprintf("[%s] Wrote: %s", .ts(), path))
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # dependency checks (Imports)
  .need_pkgs <- c("data.table", "fastICA", "mixOmics")
  .miss <- .need_pkgs[! .need_pkgs %in% rownames(utils::installed.packages())]
  if (length(.miss) > 0) {
    stop(sprintf("Missing packages: %s. Please install them.", paste(.miss, collapse = ", ")))
  }
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  # -------- 0) Load data: prefer dir-mode via fh_load_data() --------
  prefix <- NULL
  expr_t <- NULL  # legacy-only raw matrix before scaling
  meta_dt <- NULL
  X <- NULL
  
  if (!is.null(dir)) {
    if (is.null(train_exp_name) || is.null(positive_label)) {
      stop("When 'dir' is provided, both 'train_exp_name' and 'positive_label' must be set.")
    }
    if (!exists("fh_load_data", mode = "function")) {
      stop("fh_load_data() not found. Please ensure it is available in your package namespace.")
    }
    dl <- fh_load_data(
      dir = dir,
      train_exp_name = train_exp_name,
      hub_file = hub_file,
      positive_label = positive_label,
      assign_to_global = FALSE
    )
    X <- dl$x_train                                   # samples x genes
    meta_dt <- data.table::data.table(
      Sample = rownames(X),
      Group  = as.character(dl$y_train)               # "0"/"1"
    )
    prefix <- tools::file_path_sans_ext(basename(train_exp_name))
  } else {
    # -------- Legacy mode: expr_path + meta_path --------
    if (is.null(expr_path) || is.null(meta_path)) {
      stop("Provide either (dir, train_exp_name, positive_label) or (expr_path, meta_path).")
    }
    if (!file.exists(expr_path)) stop("expr_path does not exist.")
    if (!file.exists(meta_path)) stop("meta_path does not exist.")
    
    expr_dt <- data.table::fread(expr_path, sep = "\t", header = TRUE)
    meta_dt <- data.table::fread(meta_path, sep = "\t", header = FALSE,
                                 col.names = c("Sample", "Group"))
    
    if (ncol(expr_dt) < 2) stop("Expression file must have gene ID column plus at least one sample column.")
    gene_ids <- expr_dt[[1]]
    expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
    rownames(expr_mat) <- gene_ids
    expr_t <- t(expr_mat)                        # samples x genes
    rownames(expr_t) <- colnames(expr_mat)
    
    common_samples <- intersect(rownames(expr_t), meta_dt$Sample)
    if (length(common_samples) == 0) stop("No overlapping samples between expression and metadata.")
    idx <- match(meta_dt$Sample, rownames(expr_t))
    expr_t <- expr_t[idx, , drop = FALSE]
    rownames(expr_t) <- meta_dt$Sample
    keep <- rowSums(is.na(expr_t)) == 0
    expr_t <- expr_t[keep, , drop = FALSE]
    meta_dt <- meta_dt[keep, ]
    
    if (!is.null(group_levels)) {
      meta_dt$Group <- factor(meta_dt$Group, levels = group_levels)
    }
    X <- expr_t
    prefix <- tools::file_path_sans_ext(basename(expr_path))
  }
  
  # ---------- 1) Guards and scaling ----------
  if (nrow(X) < 4L) stop("Too few samples (<4).")
  if (ncol(X) < min_genes) {
    message(sprintf("[%s] Number of genes = %d < min_genes = %d. Skipping reduction.",
                    .ts(), ncol(X), min_genes))
    return(invisible(list(skipped = TRUE, reason = "insufficient_genes",
                          n_samples = nrow(X), n_genes = ncol(X))))
  }
  
  expr_raw_before_scale <- expr_t  # NULL in dir-mode; matrix in legacy-mode
  X <- if (normalize) base::scale(X, center = TRUE, scale = TRUE) else X
  if (is.null(colnames(X))) stop("Expression matrix has no gene names.")
  colnames(X) <- make.unique(colnames(X))
  
  max_comps <- min(as.integer(max_components), ncol(X), nrow(X) - 1L)
  if (max_comps < 2L) stop("Not enough degrees of freedom for at least 2 components.")
  
  # y factor for PLS-DA
  y_fac <- factor(meta_dt$Group)
  
  # helpers for output filenames
  .fp <- function(name) file.path(out_dir, paste0(prefix, "_", name))
  
  # ---------- 2) PCA ----------
  pca_fit <- stats::prcomp(X, center = FALSE, scale. = FALSE, rank. = max_comps)
  pca_scores <- pca_fit$x
  colnames(pca_scores) <- paste0("PC", seq_len(ncol(pca_scores)))
  rownames(pca_scores) <- rownames(X)
  pca_scores_dt <- data.table::data.table(Sample = rownames(pca_scores), as.data.frame(pca_scores))
  p_pca_scores <- .fp("PCA_scores.csv"); data.table::fwrite(pca_scores_dt, p_pca_scores); .msg_write(p_pca_scores)
  
  pca_var <- (pca_fit$sdev ^ 2)
  pca_var_ratio <- pca_var / sum(pca_var)
  pca_var_cum <- cumsum(pca_var_ratio)
  pca_var_dt <- data.table::data.table(
    Component = paste0("PC", seq_len(length(pca_var_ratio))),
    Var_Explained = as.numeric(pca_var_ratio),
    Var_Explained_Cum = as.numeric(pca_var_cum)
  )
  p_pca_var <- .fp("PCA_var_explained.csv"); data.table::fwrite(pca_var_dt, p_pca_var); .msg_write(p_pca_var)
  
  # Scree plot (PNG)
  p_png <- .fp("PCA_scree.png")
  grDevices::png(filename = p_png, width = 1600, height = 900, res = 150)
  op <- graphics::par(mar = c(5,5,2,1))
  graphics::plot(seq_along(pca_var_ratio), pca_var_ratio, type = "b",
                 xlab = "Component", ylab = "Proportion of Variance Explained",
                 main = "PCA Scree")
  graphics::lines(seq_along(pca_var_cum), pca_var_cum, lty = 2)
  graphics::legend("bottomright", legend = c("PVE", "Cumulative PVE"), lty = c(1,2), bty = "n")
  graphics::par(op)
  grDevices::dev.off()
  .msg_write(p_png)
  
  # ---------- 3) ICA ----------
  ica_fit <- fastICA::fastICA(X, n.comp = max_comps, method = ica_method)
  ica_scores <- ica_fit$S
  colnames(ica_scores) <- paste0("IC", seq_len(ncol(ica_scores)))
  rownames(ica_scores) <- rownames(X)
  ica_scores_dt <- data.table::data.table(Sample = rownames(ica_scores), as.data.frame(ica_scores))
  p_ica_scores <- .fp("ICA_scores.csv"); data.table::fwrite(ica_scores_dt, p_ica_scores); .msg_write(p_ica_scores)
  
  # ---------- 4) PLS-DA ----------
  ncomp_cap <- min(30L, nrow(X) - 1L, ncol(X))
  ncomp_use <- max(2L, min(as.integer(pls_ncomp_default), ncomp_cap))
  fit_plsda <- mixOmics::plsda(X, y_fac, ncomp = ncomp_use)
  pls_scores <- fit_plsda$variates$X
  colnames(pls_scores) <- paste0("PLS", seq_len(ncol(pls_scores)))
  rownames(pls_scores) <- rownames(X)
  pls_scores_dt <- data.table::data.table(Sample = rownames(pls_scores), as.data.frame(pls_scores))
  p_pls_scores <- .fp("PLS_scores.csv"); data.table::fwrite(pls_scores_dt, p_pls_scores); .msg_write(p_pls_scores)
  message(sprintf("[%s] mixOmics::plsda done.", .ts()))
  
  # ---------- 5) Concatenate with meta ----------
  stopifnot(all(c("Sample","Group") %in% colnames(meta_dt)))
  merge12  <- merge(pca_scores_dt, ica_scores_dt, by = "Sample", all = TRUE)
  merge123 <- merge(merge12,       pls_scores_dt, by = "Sample", all = TRUE)
  merged   <- merge(meta_dt, merge123, by = "Sample", all.y = TRUE)
  
  feat_cols <- setdiff(colnames(merged), c("Sample","Group"))
  data.table::setcolorder(merged, c("Sample","Group", feat_cols))
  data.table::setorderv(merged, c("Group","Sample"))
  
  p_concat <- .fp("reduced_concat.csv"); data.table::fwrite(merged, p_concat); .msg_write(p_concat)
  
  # ---------- 6) Gene ranking ----------
  genes <- colnames(X)
  K_pca <- min(10L, ncol(pca_fit$rotation))
  K_ica <- min(10L, ncol(ica_scores))
  K_pls <- min(10L, ncol(pls_scores))
  
  pca_load <- pca_fit$rotation[, seq_len(K_pca), drop = FALSE]
  pca_w    <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  pca_wK   <- pca_w[seq_len(K_pca)]
  pca_imp  <- as.numeric(abs(pca_load) %*% matrix(pca_wK, ncol = 1)); names(pca_imp) <- genes
  
  if (is.null(ica_fit$A)) stop("fastICA returned NULL mixing matrix A.")
  ica_A <- t(ica_fit$A)[, seq_len(K_ica), drop = FALSE]
  ica_imp <- rowSums(abs(ica_A)); names(ica_imp) <- genes
  
  pls_load <- fit_plsda$loadings$X[, seq_len(K_pls), drop = FALSE]
  pls_imp  <- rowSums(abs(pls_load)); names(pls_imp) <- genes
  
  pca_z <- base::scale(pca_imp)[,1]
  ica_z <- base::scale(ica_imp)[,1]
  pls_z <- base::scale(pls_imp)[,1]
  combo <- (pca_z + ica_z + pls_z) / 3
  
  rank_dt <- data.table::data.table(
    Gene = genes,
    PCA_importance = as.numeric(pca_imp),
    ICA_importance = as.numeric(ica_imp),
    PLS_importance = as.numeric(pls_imp),
    PCA_z = as.numeric(pca_z),
    ICA_z = as.numeric(ica_z),
    PLS_z = as.numeric(pls_z),
    Combo_z = as.numeric(combo)
  )
  # robust in-place sort by a single column name (no NSE)
  data.table::setorderv(rank_dt, "Combo_z", order = -1L)
  
  p_rank_all <- .fp("gene_ranking_all.csv"); data.table::fwrite(rank_dt, p_rank_all); .msg_write(p_rank_all)
  
  keep_dt <- rank_dt[rank_dt[["Combo_z"]] > 0, ]
  p_rank_pos <- .fp("gene_ranking_combo_pos.csv"); data.table::fwrite(keep_dt, p_rank_pos); .msg_write(p_rank_pos)
  
  keep_genes <- keep_dt$Gene
  if (length(keep_genes) == 0L) {
    message(sprintf("[%s] No genes with Combo_z > 0.", .ts()))
  } else {
    keep_genes <- keep_genes[keep_genes %in% colnames(X)]
    
    # scaled subset (always available)
    expr_scl_keep <- X[, keep_genes, drop = FALSE]
    dt_scl <- data.table::data.table(Sample = rownames(expr_scl_keep), as.data.frame(expr_scl_keep))
    p_expr_scl <- .fp("expr_combo_pos_scaled.tsv"); data.table::fwrite(dt_scl, p_expr_scl, sep = "\t"); .msg_write(p_expr_scl)
    
    # raw subset (legacy-only)
    if (!is.null(expr_raw_before_scale)) {
      expr_raw_keep <- expr_raw_before_scale[, keep_genes, drop = FALSE]
      dt_raw <- data.table::data.table(Sample = rownames(expr_raw_keep), as.data.frame(expr_raw_keep))
      p_expr_raw <- .fp("expr_combo_pos_raw.tsv"); data.table::fwrite(dt_raw, p_expr_raw, sep = "\t"); .msg_write(p_expr_raw)
    }
    
    p_gene_txt <- .fp("combo_pos_gene_list.txt")
    data.table::fwrite(data.table::data.table(Gene = keep_genes), p_gene_txt, sep = "\t", col.names = FALSE)
    .msg_write(p_gene_txt)
    message(sprintf("[%s] Kept %d genes with Combo_z > 0.", .ts(), length(keep_genes)))
  }
  
  message(sprintf("[%s] Done. Outputs written under: %s", .ts(), out_dir))
  
  invisible(list(
    out_dir = out_dir,
    prefix = prefix,
    n_samples = nrow(X),
    n_genes = ncol(X),
    max_components = max_comps
  ))
}