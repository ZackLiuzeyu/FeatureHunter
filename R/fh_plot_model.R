#' Pick a ranked model from leaderboard, match in `all_result_summary`, and plot
#'
#' @description
#' Selects the \code{pick_index}-th row from a leaderboard CSV (by rank),
#' locates the corresponding key in \code{all_result_summary}, prints
#' confusion matrices and ACC/F1/Recall, and saves PR curves, probability
#' histograms, and density plots. If the selected row appears to be an MLP
#' model (name starts with "NN-MLP ("), the function will also attempt to parse
#' hyperparameters to help match the key if exact string matching fails.
#'
#' @param all_result_summary list Nested list of per-model, per-dataset
#'   prediction data.frames. Each dataset data.frame must contain
#'   \code{predict_p}, \code{real_label}, and \code{predict_result}.
#' @param all_result_FS list Optional parallel list used only to retrieve
#'   dataset names/order when \code{ds_names} is not supplied.
#' @param nshow integer Used to build default leaderboard CSV filename.
#' @param namesS character Score names used in default filename construction,
#'   e.g., \code{c("Accuracy","Recall","F-score")}.
#' @param score_index integer Which score to use for default filename; e.g.,
#'   \code{1=Accuracy, 2=Recall, 3=F-score}.
#' @param pick_index integer 1-based index selecting the ranked row in the
#'   leaderboard CSV (no longer restricted to rows containing "MLP").
#' @param top_models_csv character Optional explicit path to leaderboard CSV;
#'   if \code{NULL}, filename is constructed from \code{nshow}, \code{namesS},
#'   and \code{score_index}.
#' @param key_pattern_prefix character Regex prefix used only when attempting
#'   to fuzzy-match MLP keys if exact matching fails. Default \code{"^NN-MLP \\("}.
#' @param learningratei,batch_sizei,epochselecti,dropoutratei,cutoffi
#'   numeric/integer/NULL Manual overrides; if \code{NULL} and the selected
#'   row is MLP, values are parsed from the leaderboard string.
#' @param ds_names character Optional dataset name order; if provided, it is
#'   enforced. Otherwise the function tries to reuse names from
#'   \code{all_result_FS[[key]]} when available; otherwise existing names are used.
#' @param cm_palette_low,cm_palette_mid,cm_palette_high character Colors for
#'   confusion-matrix heatmap.
#' @param base_size numeric Base font size for plots.
#'
#' @return list with elements:
#' \itemize{
#'   \item key matched model key (string)
#'   \item params parsed/overridden hyperparameters and file info
#'   \item metrics per-dataset confusion matrices + ACC/F1/Recall
#'   \item plots ggplot objects (confusion matrices, PR curves, histograms, density plots)
#' }
#'
#' @details
#' Exact string matching from leaderboard \code{Model} to
#' \code{names(all_result_summary)} is attempted first. If it fails and the
#' selected name appears to be an MLP, hyperparameters are parsed to narrow
#' candidates by regex. PR curves are saved to \code{pr_plots/}; histograms and
#' density plots to \code{prob_plots/}; confusion matrices to \code{confusion_matrix/}.
#'
#' @import ggplot2 precrec
#' @examples
#' \dontrun{
#' # pick the 1st ranked row (any model type), F-score leaderboard:
#' out <- fh_plot_model(
#'   all_result_summary = all_result_summary,
#'   all_result_FS      = all_result_FS,
#'   nshow        = 40,
#'   namesS       = c("Accuracy","Recall","F-score"),
#'   score_index  = 3,
#'   pick_index   = 1,
#'   cm_palette_low  = "#20ACBD",
#'   cm_palette_mid  = "#FCECCF",
#'   cm_palette_high = "#F69896",
#'   base_size = 11
#' )
#' out$key
#' out$metrics$`DatasetB(Val)`$acc
#'
#' # pick the 7th row explicitly from a given CSV:
#' out2 <- fh_plot_model(
#'   all_result_summary = all_result_summary,
#'   all_result_FS      = all_result_FS,
#'   top_models_csv = "heatmap/40_top_F-score_models.csv",
#'   pick_index   = 7
#' )
#' }
#' @export
fh_plot_model <- function(
    all_result_summary,
    all_result_FS,
    nshow,
    namesS = c("Accuracy","Recall","F-score"),
    score_index = 3,
    pick_index = 1,
    top_models_csv = NULL,
    key_pattern_prefix = "^NN-MLP \\(",
    learningratei = NULL,
    batch_sizei   = NULL,
    epochselecti  = NULL,
    dropoutratei  = NULL,
    cutoffi       = NULL,
    ds_names      = NULL,
    cm_palette_low  = "#20ACBD",
    cm_palette_mid  = "#FCECCF",
    cm_palette_high = "#F69896",
    base_size = 11
){
  stopifnot(!missing(all_result_summary))
  stopifnot(!missing(nshow))
  if (!(score_index %in% seq_along(namesS))) {
    stop("score_index out of bounds.")
  }
  if (is.null(top_models_csv)) {
    top_models_csv <- file.path("heatmap", paste0(nshow, "_top_", namesS[score_index], "_models.csv"))
  }
  if (!file.exists(top_models_csv)) {
    stop("Leaderboard CSV not found: ", top_models_csv)
  }
  
  Top_models <- utils::read.csv(top_models_csv, stringsAsFactors = FALSE)
  if (!("Model" %in% names(Top_models))) stop("Leaderboard CSV missing column 'Model'.")
  if (nrow(Top_models) < pick_index || pick_index < 1) {
    stop("pick_index out of range (nrow=", nrow(Top_models), ").")
  }
  
  model_str <- Top_models$Model[pick_index]
  
  # ---- primary: exact key match ----
  keys <- names(all_result_summary)
  if (model_str %in% keys) {
    matched_key <- model_str
  } else {
    # ---- fallback: if looks like MLP, try to parse hyperparams and fuzzy-match ----
    if (grepl("^NN-MLP \\(", model_str)) {
      extract_param <- function(pattern, text) {
        m <- regmatches(text, regexpr(pattern, text))
        if (length(m) == 0) return(NA_real_)
        as.numeric(gsub("[^0-9\\.]", "", m))
      }
      lr_p     <- extract_param("lr:[0-9\\.]+",        model_str)
      bs_p     <- extract_param("bs:[0-9]+",           model_str)
      ep_p     <- extract_param("ep:[0-9]+",           model_str)
      dr_p     <- extract_param("dropout:[0-9\\.]+",   model_str)
      co_auto1 <- regmatches(model_str, regexec("cutoff:auto\\(([0-9\\.]+),", model_str))[[1]]
      if (length(co_auto1) >= 2) {
        co_p <- as.numeric(co_auto1[2])
      } else {
        co_p <- extract_param("cutoff:(?:auto\\(([0-9\\.]+)\\)|([0-9\\.]+))", model_str)
        co_p <- stats::na.omit(co_p)[1]
      }
      # allow manual overrides to supersede
      learningratei <- if (is.null(learningratei)) lr_p else learningratei
      batch_sizei   <- if (is.null(batch_sizei))   bs_p else batch_sizei
      epochselecti  <- if (is.null(epochselecti))  ep_p else epochselecti
      dropoutratei  <- if (is.null(dropoutratei))  dr_p else dropoutratei
      cutoffi       <- if (is.null(cutoffi))       co_p else cutoffi
      
      # build candidate set by regex filters
      cand <- grep(key_pattern_prefix, keys, value = TRUE)
      fil  <- function(v, patt) v[grepl(patt, v, perl = TRUE)]
      if (!is.na(learningratei)) cand <- fil(cand, paste0("lr:", sub("\\.?0+$","", as.character(learningratei)), "(\\b|\\D)"))
      if (!is.na(batch_sizei))   cand <- fil(cand, paste0("bs:", batch_sizei, "(\\b|\\D)"))
      if (!is.na(epochselecti))  cand <- fil(cand, paste0("ep:", epochselecti, "(\\b|\\D)"))
      if (!is.na(dropoutratei))  cand <- fil(cand, paste0("dropout:", sub("\\.?0+$","", as.character(dropoutratei)), "(\\b|\\D)"))
      if (!is.na(cutoffi)) {
        fmt3  <- function(x) format(round(x, 3), nsmall = 3, trim = TRUE)
        co_s  <- fmt3(cutoffi)
        if (grepl("\\.[0-9]*[1-9]$", co_s)) co_num <- co_s else {
          co_num <- sub("0+$", "", co_s); co_num <- sub("\\.$", "", co_num)
        }
        patt_vec <- c(
          paste0("cutoff:auto(", co_s, ", youden)"),
          paste0("cutoff:auto(", co_s, ", f1)"),
          paste0("cutoff:", co_num)
        )
        keep <- Reduce(`|`, lapply(patt_vec, function(p) grepl(p, cand, fixed = TRUE)))
        cand <- cand[keep]
      }
      if (!length(cand)) stop("Could not match model key by exact or MLP fuzzy rules.")
      if (length(cand) > 1) {
        message("Multiple candidates matched; using the first:\n", paste0(" - ", cand, collapse = "\n"))
      }
      matched_key <- cand[1]
    } else {
      stop("Model string not found in all_result_summary and not an MLP name: ", model_str)
    }
  }
  
  if (!(matched_key %in% names(all_result_summary))) {
    stop("Matched key not found in all_result_summary (unexpected).")
  }
  lst_raw <- all_result_summary[[matched_key]]
  stopifnot(is.list(lst_raw), length(lst_raw) >= 2)
  
  # dataset names/order
  if (!is.null(ds_names)) {
    if (length(ds_names) != length(lst_raw)) {
      stop("ds_names length must equal the number of datasets in the matched entry.")
    }
    names(lst_raw) <- ds_names
  } else if (!missing(all_result_FS) &&
             !is.null(all_result_FS) &&
             (matched_key %in% names(all_result_FS))) {
    tmp <- names(all_result_FS[[matched_key]])
    if (!is.null(tmp) && length(tmp) == length(lst_raw)) {
      names(lst_raw) <- tmp
    } else if (is.null(names(lst_raw))) {
      stop("Dataset names missing and cannot be inferred from all_result_FS; please provide ds_names.")
    }
  } else if (is.null(names(lst_raw))) {
    stop("Dataset names missing; please provide ds_names.")
  }
  lst_named <- lst_raw
  
  to01 <- function(x) {
    x <- tolower(as.character(x))
    out <- ifelse(x %in% c("positive","pos","1","yes","true"), 1,
                  ifelse(x %in% c("negative","neg","0","no","false"), 0, NA))
    as.numeric(out)
  }
  metrics_from_df <- function(df) {
    pr <- if (is.factor(df$predict_result)) as.character(df$predict_result) else df$predict_result
    rl <- if (is.factor(df$real_label))    as.character(df$real_label)    else df$real_label
    pr <- ifelse(pr %in% c("positive","pos","1",1,TRUE), 1L, 0L)
    rl <- ifelse(rl %in% c("positive","pos","1",1,TRUE), 1L, 0L)
    TP <- sum(pr==1 & rl==1); FP <- sum(pr==1 & rl==0)
    FN <- sum(pr==0 & rl==1); TN <- sum(pr==0 & rl==0)
    ACC <- (TP+TN)/length(rl)
    PREC<- if ((TP+FP)==0) 0 else TP/(TP+FP)
    REC <- if ((TP+FN)==0) 0 else TP/(TP+FN)
    F1  <- if ((PREC+REC)==0) 0 else 2*PREC*REC/(PREC+REC)
    cm <- matrix(c(TN,FP,FN,TP), nrow=2, byrow=TRUE,
                 dimnames=list(pred=c("neg","pos"), real=c("neg","pos")))
    list(cm=cm, acc=ACC, f1=F1, rec=REC)
  }
  
  cat("===== CONFUSION MATRICES from all_result_summary =====\n")
  metrics <- list()
  for (nm in names(lst_named)) {
    out <- metrics_from_df(lst_named[[nm]])
    cat("\n-- ", nm, " --\n", sep="")
    print(out$cm)
    cat(sprintf("  ACC=%.3f  F1=%.3f  REC=%.3f\n", out$acc, out$f1, out$rec))
    metrics[[nm]] <- out
  }
  
  plot_cm <- function(cm, title = "", subtitle = NULL,
                      palette_low  = cm_palette_low,
                      palette_mid  = cm_palette_mid,
                      palette_high = cm_palette_high,
                      base_size    = 11) {
    df <- as.data.frame(as.table(cm))
    colnames(df)[1:2] <- c("pred", "actual")
    df$pred   <- factor(df$pred,   levels = rev(rownames(cm)))
    df$actual <- factor(df$actual, levels = colnames(cm))
    vmin <- min(df$Freq, na.rm = TRUE)
    vmax <- max(df$Freq, na.rm = TRUE)
    vmid <- (vmin + vmax) / 2
    
    ggplot2::ggplot(df, ggplot2::aes(x = pred, y = actual, fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = Freq),
                         color = "black", size = base_size * 0.5) +
      ggplot2::scale_fill_gradient2(low = palette_low, mid = palette_mid, high = palette_high,
                                    midpoint = vmid, limits = c(vmin, vmax), name = "Count") +
      ggplot2::coord_fixed() +
      ggplot2::labs(x = "Predicted", y = "Actual", title = title, subtitle = subtitle) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size * 0.9),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
        axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(r = 6))
      )
  }
  
  cm_dir <- "confusion_matrix"
  if (!dir.exists(cm_dir)) dir.create(cm_dir, showWarnings = FALSE, recursive = TRUE)
  
  plots_cm <- list()
  for (nm in names(lst_named)) {
    out <- metrics[[nm]]
    title_str <- sprintf("%s\nACC=%.3f  F1=%.3f  REC=%.3f\n%s",
                         nm, out$acc, out$f1, out$rec, matched_key)
    p <- plot_cm(out$cm, title = title_str, subtitle = NULL, base_size = base_size)
    plots_cm[[nm]] <- p
    fname <- paste0("Confusion_Matrix_", gsub("[^A-Za-z0-9_\\-\\.]", "_", nm), ".pdf")
    fpath <- file.path(cm_dir, fname)
    ggplot2::ggsave(filename = fpath, plot = p, width = 6, height = 6)
    abspath <- normalizePath(fpath, winslash = "/", mustWork = FALSE)
    message(sprintf("[%s] Saved plot to: %s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), abspath))
  }
  message("Finished saving Confusion Matrix plots.")
  
  if (!dir.exists("pr_plots"))   dir.create("pr_plots",   showWarnings = FALSE)
  if (!dir.exists("prob_plots")) dir.create("prob_plots", showWarnings = FALSE)
  
  to_bin <- function(x) {
    x <- tolower(as.character(x))
    out <- ifelse(x %in% c("positive","pos","1","yes","true"), 1,
                  ifelse(x %in% c("negative","neg","0","no","false"), 0, NA))
    as.numeric(out)
  }
  
  plots_pr  <- list()
  plots_hist<- list()
  plots_den <- list()
  
  for (nm in names(lst_named)) {
    df <- lst_named[[nm]]
    if (!all(c("predict_p","real_label") %in% names(df))) {
      message("[skip] ", nm, ": missing predict_p or real_label")
      next
    }
    score <- as.numeric(df$predict_p)
    label <- to_bin(df$real_label)
    keep  <- !(is.na(score) | is.na(label))
    score <- score[keep]; label <- label[keep]
    
    if (length(unique(label)) < 2) {
      message("[skip] ", nm, ": single-class labels, cannot draw PR. n=", length(label))
      next
    }
    
    mm <- precrec::mmdata(scores = score, labels = label, modnames = nm, dsids = 1)
    ev <- precrec::evalmod(mm, mode = "pr")
    auc_tbl <- precrec::auc(ev)
    ap <- if ("curvetypes" %in% names(auc_tbl)) {
      auc_tbl[auc_tbl$curvetypes %in% c("PRC","PR","prc"), "aucs"][1]
    } else {
      auc_tbl$aucs[1]
    }
    
    df_all <- as.data.frame(ev)
    type_col <- intersect(names(df_all), c("curvetypes","ctype","curvetype","curve_type","type"))[1]
    if (!is.na(type_col)) {
      df_pr <- df_all[df_all[[type_col]] %in% c("PRC","PR","prc","precision-recall","precision_recall"), , drop = FALSE]
    } else {
      df_pr <- df_all
    }
    if (!("recall" %in% names(df_pr))) {
      if (all(c("x","y") %in% names(df_pr))) {
        df_pr$recall    <- as.numeric(df_pr$x)
        df_pr$precision <- as.numeric(df_pr$y)
      } else {
        r_col <- intersect(names(df_pr), c("recall","tpr","sens","sensitivity","rec"))
        p_col <- intersect(names(df_pr), c("precision","ppv","prec"))
        if (length(r_col) && length(p_col)) {
          df_pr$recall    <- as.numeric(df_pr[[r_col[1]]])
          df_pr$precision <- as.numeric(df_pr[[p_col[1]]])
        } else {
          stop("Could not identify PR columns; expected recall/precision or x/y.")
        }
      }
    }
    
    pos_prev <- mean(label == 1L)
    pr_plot <- ggplot2::ggplot(df_pr, ggplot2::aes(x = recall, y = precision)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_hline(yintercept = pos_prev, linetype = 2) +
      ggplot2::annotate("text", x = 0.03, y = pos_prev,
                        label = sprintf("baseline = %.2f", pos_prev),
                        hjust = 0, vjust = -0.6, size = 3) +
      ggplot2::coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      ggplot2::labs(title = paste0("PR Curve - ", nm),
                    subtitle = sprintf("AUC-PR = %.3f", ap),
                    x = "Recall", y = "Precision") +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5)
      )
    plots_pr[[nm]] <- pr_plot
    pr_fname <- paste0("pr_plots/PR_", gsub("[^A-Za-z0-9_\\-\\.]", "_", nm), ".pdf")
    ggplot2::ggsave(filename = pr_fname, plot = pr_plot, width = 5.2, height = 4.4)
    pr_path <- normalizePath(pr_fname, winslash = "/", mustWork = FALSE)
    message(sprintf("[%s] Saved plot to: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), pr_path))
    
    dplot <- ggplot2::ggplot(
      data = data.frame(score = score, label = factor(label, labels = c("Negative","Positive"))),
      ggplot2::aes(x = score, fill = label)
    ) +
      ggplot2::geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
      ggplot2::coord_cartesian(xlim = c(0,1)) +
      ggplot2::labs(title = paste0("Predicted Probability Histogram - ", nm),
                    x = "Predicted probability", y = "Count", fill = "real") +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    plots_hist[[nm]] <- dplot
    hist_fname <- paste0("prob_plots/Hist_", gsub("[^A-Za-z0-9_\\-\\.]", "_", nm), ".pdf")
    ggplot2::ggsave(filename = hist_fname, plot = dplot, width = 5.2, height = 4.4)
    hist_path <- normalizePath(hist_fname, winslash = "/", mustWork = FALSE)
    message(sprintf("[%s] Saved plot to: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), hist_path))
    
    kdplot <- ggplot2::ggplot(
      data = data.frame(score = score, label = factor(label, labels = c("Negative","Positive"))),
      ggplot2::aes(x = score, color = label, fill = label)
    ) +
      ggplot2::geom_density(alpha = 0.25, adjust = 1) +
      ggplot2::coord_cartesian(xlim = c(0,1)) + 
      ggplot2::labs(title = paste0("Predicted Probability Density - ", nm),
                    x = "Predicted probability", y = "Density",
                    color = "real", fill = "real") +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    plots_den[[nm]] <- kdplot
    dens_fname <- paste0("prob_plots/Density_", gsub("[^A-Za-z0-9_\\-\\.]", "_", nm), ".pdf")
    ggplot2::ggsave(filename = dens_fname, plot = kdplot, width = 5.2, height = 4.4)
    dens_path <- normalizePath(dens_fname, winslash = "/", mustWork = FALSE)
    message(sprintf("[%s] Saved plot to: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), dens_path))
    
    message("Finished saving PR/Histogram/Density plots for ", nm)
  }
  
  invisible(list(
    key = matched_key,
    params = list(
      top_models_csv = top_models_csv,
      score_name = namesS[score_index],
      pick_index = pick_index,
      learningratei = learningratei,
      batch_sizei   = batch_sizei,
      epochselecti  = epochselecti,
      dropoutratei  = dropoutratei,
      cutoffi       = cutoffi,
      key_pattern_prefix = key_pattern_prefix,
      model_str = model_str
    ),
    metrics = metrics,
    plots = list(
      confusion_matrices = plots_cm,
      pr = plots_pr,
      hist = plots_hist,
      density = plots_den
    )
  ))
}