#' Plot ROC and PR curves for top models from leaderboard CSV (using precrec for PR)
#'
#' @description
#' Reads a leaderboard CSV, selects models (either the first \code{n_model}
#' entries or a specific \code{pick_index}), then for each model draws ROC and
#' PR curves across available datasets using prediction probabilities stored in
#' \code{all_result_summary}. Dataset display names are derived from \code{exp_list}.
#' PR curves and AP are computed via \code{precrec} for consistency with standard practice.
#' Outputs a multi-page PDF with one model per page (left: ROC, right: PR).
#'
#' 
#' @param n_model Integer. Number of leaderboard rows to consider (default 16).
#' @param top_models_csv Character or NULL. Path to the leaderboard CSV.
#'   If NULL, a default path is constructed as
#'   file.path("heatmap", sprintf("%d_top_%s_models.csv", n_model, score_name)).
#' @param score_index Integer. Must be 1, 2, or 3. Selects which score name to
#'   use when building the default CSV filename (1 = "Accuracy",
#'   2 = "Recall", 3 = "F-score").
#' @param pick_index Optional integer. If provided, selects exactly that row
#'   (1-based) from the leaderboard instead of the first n_model.
#' @param nshow Integer. Number of leaderboard rows (defined when drawing heatmap).
#' @param out_pdf Character. Output PDF path (default "top_confusion_curves.pdf").
#'
#' @return (Invisibly) the output PDF path.
#' @export
fh_plot_PRROC <- function(n_model = 16,
                          top_models_csv = NULL,
                          score_index = 1L,
                          pick_index = NULL,
                          nshow = nshow,
                          out_pdf = "top_PR-ROC.pdf") {
  # deps
  if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("dplyr",    quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("gridExtra",quietly = TRUE)) stop("Package 'gridExtra' is required.")
  if (!requireNamespace("precrec",  quietly = TRUE)) stop("Package 'precrec' is required.")
  
  # locate required objects
  results_obj <- try(get("all_result_summary", envir = parent.frame()), silent = TRUE)
  if (inherits(results_obj, "try-error") || is.null(results_obj)) {
    results_obj <- try(get("all_result_summary", inherits = TRUE), silent = TRUE)
  }
  if (inherits(results_obj, "try-error") || is.null(results_obj)) {
    stop("Could not find 'all_result_summary' in calling or parent environments.")
  }
  
  exp_list_obj <- try(get("exp_list", envir = parent.frame()), silent = TRUE)
  if (inherits(exp_list_obj, "try-error") || is.null(exp_list_obj)) {
    exp_list_obj <- try(get("exp_list", inherits = TRUE), silent = TRUE)
  }
  if (inherits(exp_list_obj, "try-error") || is.null(exp_list_obj)) {
    stop("Could not find 'exp_list' in calling or parent environments.")
  }
  
  # derive dataset display names from exp_list length
  set_names_pref <- .make_set_names_pref(exp_list_obj)
  
  # score name mapping for default CSV path
  .score_names <- c("Accuracy","Recall","F-score")
  if (!is.null(top_models_csv)) {
    score_name <- NULL
  } else {
    if (!(score_index %in% 1:3)) stop("score_index must be 1, 2, or 3.")
    score_name <- .score_names[score_index]
    # FIX: use n_model here (was nshow)
    top_models_csv <- file.path("heatmap", sprintf("%d_top_%s_models.csv", nshow, score_name))
  }
  
  if (!file.exists(top_models_csv)) {
    stop("Leaderboard CSV not found: ", top_models_csv)
  }
  
  Top_models <- utils::read.csv(top_models_csv, stringsAsFactors = FALSE, check.names = FALSE)
  # prefer a column that contains "model", otherwise use the first column
  model_col <- names(Top_models)[grepl("model", names(Top_models), ignore.case = TRUE)]
  if (length(model_col) == 0) model_col <- names(Top_models)[1]
  if (!(model_col %in% names(Top_models))) {
    stop("Leaderboard CSV: cannot identify model column.")
  }
  if (nrow(Top_models) == 0) stop("Leaderboard CSV has no rows.")
  
  # select model keys
  model_vec <- unique(as.character(Top_models[[model_col]]))
  model_keys <- if (!is.null(pick_index)) {
    if (pick_index < 1L || pick_index > length(model_vec)) {
      stop("pick_index is out of range.")
    }
    model_vec[pick_index]
  } else {
    utils::head(model_vec, n_model)
  }
  
  # only keep keys present in results
  model_keys <- model_keys[model_keys %in% names(results_obj)]
  if (!length(model_keys)) {
    stop("No matching model keys found in 'all_result_summary' for the chosen rows.")
  }
  
  # write PDF
  grDevices::pdf(out_pdf, width = 12, height = 5)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  for (mk in model_keys) {
    sets_obj <- results_obj[[mk]]
    if (is.null(sets_obj)) next
    
    # name sets locally (do not mutate global object)
    if (is.null(names(sets_obj)) || any(is.na(names(sets_obj))) || any(nchar(names(sets_obj)) == 0)) {
      names(sets_obj) <- utils::head(set_names_pref, n = length(sets_obj))
    }
    sets_avail <- intersect(set_names_pref, names(sets_obj))
    if (!length(sets_avail)) next
    
    # collect curves
    roclist <- list(); prlist <- list(); notes <- list()
    for (sn in sets_avail) {
      cur <- .get_curves_for_set_precrec(sets_obj[[sn]])
      if (is.null(cur)) next
      if (!cur$constant && nrow(cur$roc)) {
        roclist[[sn]] <- dplyr::mutate(cur$roc, Set = sn)
      } else {
        notes[[sn]] <- paste0(sn, ": constant or missing (pos_rate=",
                              sprintf("%.3f", cur$pos_rate), ")")
      }
      if (!cur$constant && nrow(cur$pr)) {
        prlist[[sn]]  <- dplyr::mutate(cur$pr, Set = sn)
      }
    }
    
    roc_df <- if (length(roclist)) dplyr::bind_rows(roclist) else NULL
    pr_df  <- if (length(prlist))  dplyr::bind_rows(prlist)  else NULL
    
    # legends (AUC/AP per set)
    auc_leg <- NULL
    if (!is.null(roc_df) && nrow(roc_df)) {
      auc_leg <- roc_df |>
        dplyr::group_by(.data$Set) |>
        dplyr::summarise(AUC = utils::tail(stats::na.omit(.data$auc), 1)) |>
        dplyr::mutate(label = paste0(.data$Set, " | AUC=", sprintf("%.3f", .data$AUC)))
    }
    ap_leg <- NULL
    if (!is.null(pr_df) && nrow(pr_df)) {
      ap_leg <- pr_df |>
        dplyr::group_by(.data$Set) |>
        dplyr::summarise(AP = utils::tail(stats::na.omit(.data$ap), 1)) |>
        dplyr::mutate(label = paste0(.data$Set, " | AP=", sprintf("%.3f", .data$AP)))
    }
    
    # ROC plot
    pROC <- ggplot2::ggplot() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::ggtitle(paste0(mk, " - ROC")) +
      ggplot2::xlab("FPR") + ggplot2::ylab("TPR") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 3, linewidth = 0.4, color = "grey50") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 6))
    if (!is.null(roc_df) && nrow(roc_df)) {
      pROC <- pROC +
        ggplot2::geom_line(data = roc_df, ggplot2::aes(x = .data$fpr, y = .data$tpr, color = .data$Set),
                           linewidth = 0.8, alpha = 0.9) +
        ggplot2::coord_equal(xlim = c(0,1), ylim = c(0,1)) +
        ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
        ggplot2::theme(legend.position = "bottom")
    } else {
      pROC <- pROC + ggplot2::annotate("text", x = .5, y = .5, label = "Not available (constant or missing)")
    }
    if (!is.null(auc_leg) && nrow(auc_leg)) {
      pROC <- pROC + ggplot2::labs(subtitle = paste(auc_leg$label, collapse = "   "))
    }
    if (length(notes)) {
      pROC <- pROC + ggplot2::annotate("text", x = .98, y = .02, hjust = 1, vjust = 0,
                                       label = paste0(notes, collapse = " | "),
                                       size = 3.2, color = "firebrick")
    }
    
    # PR plot
    pPR <- ggplot2::ggplot() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::ggtitle(paste0(mk, " - PR")) +
      ggplot2::xlab("Recall") + ggplot2::ylab("Precision") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 6))
    if (!is.null(pr_df) && nrow(pr_df)) {
      pPR <- pPR +
        ggplot2::geom_line(data = pr_df, ggplot2::aes(x = .data$recall, y = .data$precision, color = .data$Set),
                           linewidth = 0.8, alpha = 0.9) +
        ggplot2::coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
        ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
        ggplot2::theme(legend.position = "bottom")
    } else {
      pPR <- pPR + ggplot2::annotate("text", x = .5, y = .5, label = "Not available (constant or missing)")
    }
    if (!is.null(ap_leg) && nrow(ap_leg)) {
      pPR <- pPR + ggplot2::labs(subtitle = paste(ap_leg$label, collapse = "   "))
    }
    
    # arrange
    gridExtra::grid.arrange(pROC, pPR, ncol = 2)
  }
  
  abspath <- normalizePath(out_pdf, winslash = "/", mustWork = FALSE)
  message(sprintf("[%s] Saved plot to: %s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), abspath))
  invisible(out_pdf)
}

# ---------- internal helpers (not exported) ----------

.make_set_names_pref <- function(exp_list) {
  n <- length(exp_list)
  if (n < 3) stop("exp_list must have at least 3 datasets.")
  labs <- character(n)
  labs[1] <- "DatasetA(Train)"
  labs[2] <- "DatasetB(Val)"
  if (n >= 3) {
    idx <- 3:n
    letters_idx <- LETTERS[idx]
    labs[idx] <- paste0("Dataset", letters_idx, "(Test)")
  }
  labs
}

.as_binary01_safe <- function(y) {
  if (is.factor(y)) y <- as.character(y)
  if (is.logical(y)) return(as.integer(y))
  if (is.numeric(y)) {
    uy <- sort(unique(stats::na.omit(y)))
    if (all(uy %in% c(0,1))) return(as.integer(y))
    if (length(uy) == 2)     return(as.integer(y == max(uy)))
    stop("Numeric labels are not binary.")
  }
  if (is.character(y)) {
    z <- tolower(trimws(y))
    pos_alias <- c("1","disease","case","positive","tumor","cancer","yes","true","pos")
    neg_alias <- c("0","normal","control","negative","no","false","neg")
    pos <- z %in% pos_alias; neg <- z %in% neg_alias
    if (!all(pos | neg)) stop("Character labels cannot be coerced to binary.")
    return(as.integer(pos))
  }
  stop("Unsupported label type.")
}

.roc_curve <- function(scores, y_true) {
  o <- order(scores, decreasing = TRUE)
  s <- scores[o]; y <- y_true[o]
  P <- sum(y == 1); N <- sum(y == 0)
  if (P == 0 || N == 0) return(data.frame(th=NA, fpr=NA, tpr=NA, auc=NA_real_)[0,])
  
  ths <- c(Inf, unique(s), -Inf)
  tp <- 0L; fp <- 0L; idx <- 1L
  out <- vector("list", length(ths))
  k <- 1L
  for (t in ths) {
    while (k <= length(s) && s[k] >= t) {
      if (y[k] == 1L) tp <- tp + 1L else fp <- fp + 1L
      k <- k + 1L
    }
    tpr <- tp / P
    fpr <- fp / N
    out[[idx]] <- c(th = t, fpr = fpr, tpr = tpr)
    idx <- idx + 1L
  }
  df <- as.data.frame(do.call(rbind, out))
  df2 <- df[order(df$fpr, df$tpr), ]
  auc <- sum(diff(df2$fpr) * (utils::head(df2$tpr, -1) + utils::tail(df2$tpr, -1)) / 2, na.rm = TRUE)
  df$auc <- auc
  df
}

# New: use precrec to compute PR curve and AP (replicate AP into a column)
.get_curves_for_set_precrec <- function(set_entry) {
  if (is.null(set_entry)) return(NULL)
  if (!all(c("predict_p","real_label") %in% names(set_entry))) return(NULL)
  
  scores <- as.numeric(set_entry$predict_p)
  y_true <- .as_binary01_safe(set_entry$real_label)
  
  if (!all(is.finite(scores)) || length(unique(scores)) <= 1) {
    return(list(
      roc = data.frame(th=NA, fpr=NA, tpr=NA, auc=NA_real_)[0,],
      pr  = data.frame(th=NA, precision=NA, recall=NA, ap=NA_real_)[0,],
      constant = TRUE,
      pos_rate = mean(y_true == 1)
    ))
  }
  
  # ROC via internal trapezoid (kept as-is)
  roc_df <- .roc_curve(scores, y_true)
  
  # PR via precrec
  mm <- precrec::mmdata(scores = scores, labels = y_true, modnames = "m", dsids = 1)
  ev <- precrec::evalmod(mm, mode = "pr")
  auc_tbl <- precrec::auc(ev)
  ap_val <- if ("curvetypes" %in% names(auc_tbl)) {
    auc_tbl[auc_tbl$curvetypes %in% c("PRC","PR","prc"), "aucs"][1]
  } else {
    auc_tbl$aucs[1]
  }
  
  df_all <- as.data.frame(ev)
  type_col <- intersect(names(df_all), c("curvetypes","ctype","curvetype","curve_type","type"))[1]
  if (!is.na(type_col)) {
    pr_df <- df_all[df_all[[type_col]] %in% c("PRC","PR","prc","precision-recall","precision_recall"),
                    , drop = FALSE]
  } else {
    pr_df <- df_all
  }
  
  # normalize column names to recall/precision
  if (!("recall" %in% names(pr_df))) {
    if (all(c("x","y") %in% names(pr_df))) {
      pr_df$recall    <- as.numeric(pr_df$x)
      pr_df$precision <- as.numeric(pr_df$y)
    } else {
      r_col <- intersect(names(pr_df), c("recall","tpr","sens","sensitivity","rec"))
      p_col <- intersect(names(pr_df), c("precision","ppv","prec"))
      if (length(r_col) && length(p_col)) {
        pr_df$recall    <- as.numeric(pr_df[[r_col[1]]])
        pr_df$precision <- as.numeric(pr_df[[p_col[1]]])
      } else {
        stop("Could not identify PR columns; expected recall/precision or x/y.")
      }
    }
  }
  pr_df <- pr_df[, c("recall","precision")]
  pr_df$ap <- as.numeric(rep(ap_val, nrow(pr_df)))
  
  list(
    roc = roc_df,
    pr  = pr_df,
    constant = FALSE,
    pos_rate = mean(y_true == 1)
  )
}