#' Export model names & batch-save per-model/per-dataset result CSVs and ROC PDFs
#'
#' @description
#' This function does two things:
#' (1) Exports *all* model names (from `all_result_summary`) to a CSV;
#' (2) Reads the ranking file `filtered_Accuracy_ordered_models.csv`, takes the top-N models,
#'     then for each selected model and each dataset, writes the result CSV and draws a ROC
#'     curve (via **pROC**). Score column is detected from common names (e.g., `predict_p`);
#'     if not found, it falls back to binarized `predict_result`.
#'
#' All outputs are written under a user-configurable base directory `out_dir`
#' (by default, the user's current working directory). Result CSVs go to
#' `file.path(out_dir, out_dir_result, ...)` and ROC PDFs go to
#' `file.path(out_dir, out_dir_roc, ...)`. The file with all model names is
#' written to `file.path(out_dir, out_modelnames_csv)`.
#'
#' @param all_result_summary list
#'   Global list of evaluation outputs. Each entry corresponds to a model name,
#'   with a list of per-dataset data.frames (each should contain `predict_result`
#'   and `real_label`; if `predict_p` exists it will be used for ROC).
#' @param filtered_models_csv_in character
#'   Ranking CSV (no header, first column = model names). The function selects top `top_n`.
#'   Default: `"filtered_Accuracy_ordered_models.csv"`.
#' @param top_n integer
#'   Number of top models to take. Your original used `raw[[1]][2:7]` (skip first row),
#'   which is kept by setting `skip_first=TRUE`. Default: `6`.
#' @param skip_first logical
#'   Whether to skip the first row of the ranking CSV (to match `2:7`). Default: `TRUE`.
#' @param namesroc character(NULL)
#'   Optional dataset names used for filenames. If `NULL`, try `names(resulti)`,
#'   otherwise fall back to `"Dataset1"`, `"Dataset2"`, ...
#' @param out_dir character
#'   Base output directory. All files and subfolders will be created under this path.
#'   Default: `getwd()`.
#' @param out_modelnames_csv character
#'   Output filename (relative to `out_dir`) for all model names. Default: `"modelnames.csv"`.
#' @param out_dir_roc character
#'   Subdirectory name (under `out_dir`) for ROC PDFs. Default: `"roc"`.
#' @param out_dir_result character
#'   Subdirectory name (under `out_dir`) for per-dataset result CSVs. Default: `"result"`.
#' @param heatmapcolor character(NULL)
#'   Not used in this function (kept only for signature consistency).
#'
#' @details
#' **Column conventions**
#' - Candidate probability/score columns: `predict_p`, `predict_prob`, `prob`, `pred_prob`, `score`,
#'   `pred_score`, `prob_positive`, `p1` (matched in this order).
#' - If none found, falls back to `predict_result` (mapped to 0/1).
#' - True labels should be in `real_label`, standardized to 0/1.
#'
#' **Filenames**
#' - Use `model_name_dataset_name` (whitespace/punctuation replaced by `_`).
#'   Results are saved as `file.path(out_dir, out_dir_result, "xxx_result.csv")`
#'   and `file.path(out_dir, out_dir_roc, "xxx_roc.pdf")`.
#'
#' @return invisible(list)
#'   Returns a list with `models_written` (vector of written model names) and
#'   `selected_models` (top N selected from ranking file).
#'
#' @importFrom pROC roc
#'
#' @export
fh_roc <- function(
    all_result_summary,
    filtered_models_csv_in = file.path("heatmap", "filtered_Accuracy_ordered_models.csv"),
    top_n        = 6,
    skip_first   = TRUE,
    namesroc     = NULL,
    out_dir             = getwd(),
    out_modelnames_csv  = "modelnames.csv",
    out_dir_roc         = "roc",
    out_dir_result      = "result",
    heatmapcolor        = NULL
){
  if (!is.list(all_result_summary) || !length(all_result_summary)) {
    stop("`all_result_summary` must be a non-empty list.")
  }
  
  ## Resolve absolute output paths
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  abs_modelnames_csv <- file.path(out_dir, out_modelnames_csv)
  abs_dir_roc    <- file.path(out_dir, out_dir_roc)
  abs_dir_result <- file.path(out_dir, out_dir_result)
  
  if (!dir.exists(abs_dir_roc))    dir.create(abs_dir_roc,    recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(abs_dir_result)) dir.create(abs_dir_result, recursive = TRUE, showWarnings = FALSE)
  
  # 1) export all model names
  modelnames <- names(all_result_summary)
  if (is.null(modelnames) || !length(modelnames)) {
    modelnames <- sprintf("model_%d", seq_along(all_result_summary))
  }
  utils::write.csv(modelnames, abs_modelnames_csv, quote = TRUE, row.names = FALSE)
  
  # 2) read ranking file, take Top-N
  if (!file.exists(filtered_models_csv_in)) {
    stop("Ranking file not found: ", filtered_models_csv_in)
  }
  raw <- utils::read.csv(filtered_models_csv_in, header = FALSE, stringsAsFactors = FALSE)
  if (!ncol(raw)) stop("Ranking CSV has zero columns.")
  if (nrow(raw) < (top_n + as.integer(skip_first))) {
    warning("Ranking CSV has fewer rows than requested; will use as many as available.")
  }
  if (skip_first) {
    idx <- seq.int(from = 2L, length.out = min(top_n, max(0, nrow(raw) - 1L)))
  } else {
    idx <- seq_len(min(top_n, nrow(raw)))
  }
  selected_models <- raw[[1]][idx]
  
  # helpers
  to01 <- function(x) {
    x <- tolower(as.character(x))
    out <- ifelse(x %in% c("positive","pos","1","yes","true"), 1,
                  ifelse(x %in% c("negative","neg","0","no","false"), 0, NA))
    as.numeric(out)
  }
  pick_score <- function(df) {
    cand <- c("predict_p","predict_prob","prob","pred_prob","score",
              "pred_score","prob_positive","p1")
    hit <- cand[cand %in% names(df)]
    if (length(hit) > 0) {
      return(as.numeric(df[[hit[1]]]))
    } else if ("predict_result" %in% names(df)) {
      return(as.numeric(to01(df$predict_result)))
    } else {
      stop("No probability/score column found, and no predict_result. Columns: ", paste(names(df), collapse = ", "))
    }
  }
  
  # iterate selected models
  for (roci in selected_models) {
    resulti <- all_result_summary[[roci]]
    if (is.null(resulti)) {
      message("[skip] Model not found: ", roci)
      next
    }
    
    # dataset names
    local_namesroc <- namesroc
    if (is.null(local_namesroc)) {
      local_namesroc <- names(resulti)
      if (is.null(local_namesroc) || !length(local_namesroc)) {
        local_namesroc <- paste0("Dataset", seq_along(resulti))
      }
    }
    
    for (dataseti in seq_along(resulti)) {
      dataset <- as.data.frame(resulti[[dataseti]])
      
      # fill default column names if missing (3 columns case)
      if (length(names(dataset)) == 0) {
        if (ncol(dataset) == 3) {
          names(dataset) <- c("predict_p","predict_result","real_label")
        } else {
          stop(sprintf("Dataset %d has no column names and %d columns, cannot auto-detect.", dataseti, ncol(dataset)))
        }
      }
      
      # safe filename
      ds_name <- if (dataseti <= length(local_namesroc)) local_namesroc[dataseti] else paste0("Dataset", dataseti)
      knamess <- gsub("[[:punct:][:space:]]", "_", paste0(roci, "_", ds_name))
      
      # write result CSV
      utils::write.csv(dataset, file.path(abs_dir_result, paste0(knamess, "_result.csv")), row.names = FALSE)
      
      # labels & score
      if (!("real_label" %in% names(dataset))) {
        stop("No real_label column found. Columns: ", paste(names(dataset), collapse = ", "))
      }
      y <- to01(dataset$real_label)
      score <- pick_score(dataset)
      
      keep <- !(is.na(y) | is.na(score))
      y <- y[keep]; score <- score[keep]
      
      if (length(unique(y)) < 2) {
        message(sprintf("[skip] %s: only one class in labels, cannot draw ROC. n=%d", knamess, length(y)))
        next
      }
      
      # ROC
      roc_obj <- pROC::roc(response = y, predictor = score, levels = c(0,1), direction = "<", quiet = TRUE)
      grDevices::pdf(file.path(abs_dir_roc, paste0(knamess, "_roc.pdf")), width = 5.5, height = 5.5)
      plot(roc_obj, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE,
           col = "black", auc.polygon.col = "#B84D64")
      grDevices::dev.off()
    }
  }
  
  invisible(list(
    models_written  = modelnames,
    selected_models = selected_models
  ))
}