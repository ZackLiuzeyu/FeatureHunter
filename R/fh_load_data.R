#' Load Expression and Label Data with Hub-Gene Filtering
#'
#' @description
#' Load expression/label files, filter by a hub-gene list, and split into train
#' (first expression file = `train_exp_name`) and external test sets.
#' Returns a list. Optionally assigns key objects to a target environment
#' (e.g., `.GlobalEnv`) via `assign_to_global = TRUE`.
#'
#' @param dir character. Directory that contains raw files.
#' @param train_exp_name character. Filename of the training expression dataset (must exist in `dir`).
#' @param hub_file character. Filename of hub-gene list (TSV, first column = gene ids). Default: "hub_gene.txt".
#' @param positive_label character. Value in label files representing the positive class.
#' @param assign_to_global logical. If TRUE, assign key objects to `target_env`. Default: FALSE.
#' @param target_env environment. Target environment for assignment when `assign_to_global = TRUE`.
#'                   Default: .GlobalEnv.
#'
#' @return list with:
#' \itemize{
#'   \item x_train: training expression matrix
#'   \item y_train: training labels (factor 0/1)
#'   \item x_tests: named list of external test matrices
#'   \item genes: used hub-gene set (intersection)
#'   \item all_labels: combined label table (standardized to 0/1 in `label`)
#'   \item files: list(allfiles, labels_files, exp_files)
#'   \item exp_list: list of processed expression matrices (after filtering/transposing)
#'   \item com_genes: intersection genes used
#' }
#' @export
fh_load_data <- function(dir,
                         train_exp_name,
                         hub_file = "Common_genes.txt",
                         positive_label,
                         assign_to_global = FALSE,
                         target_env = .GlobalEnv) {
  
  .read_table <- function(path) {
    utils::read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  }
  .read_labels <- function(path) {
    # try to read a few lines as UTF-8 to check encoding
    test_enc <- tryCatch(readLines(path, n = 5, encoding = "UTF-8"), error = function(e) NULL)
    if (is.null(test_enc)) {
      warning("This txt file may have been exported from another format and is not UTF-8 encoded. ",
              "Please copy all text into a newly created txt file and re-import.")
    }
    utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  }
  
  if (!dir.exists(dir)) stop("Directory not found: ", dir)
  allfiles <- list.files(dir, full.names = FALSE)
  if (!length(allfiles)) stop("No files found under: ", dir)
  
  # identify label files and expression files
  labels_files <- allfiles[grepl("Sample", allfiles, ignore.case = TRUE)]
  if (!length(labels_files)) stop("No label files matched pattern 'Sample' under: ", dir)
  
  if (!file.exists(file.path(dir, hub_file))) {
    stop("Hub-gene file not found: ", file.path(dir, hub_file))
  }
  
  # expression files = all - (labels + hub)
  exp_files <- setdiff(allfiles, c(labels_files, hub_file))
  if (!length(exp_files)) stop("No expression files found after excluding labels & hub file.")
  
  # put training file first
  if (!(train_exp_name %in% exp_files))
    stop("`train_exp_name` not found among expression files: ", train_exp_name)
  exp_files <- c(train_exp_name, setdiff(exp_files, train_exp_name))
  
  # read expression matrices (genes in rows, samples in columns), then transpose
  exp_list <- lapply(file.path(dir, exp_files), .read_table)
  
  # hub genes = intersection with common genes across datasets
  common_across <- Reduce(intersect, lapply(exp_list, rownames))
  hub_path <- file.path(dir, hub_file)
  first_lines <- tryCatch(readLines(hub_path, n = 5, encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(first_lines)) {
    warning("This txt file may have been exported from another format and is not UTF-8 encoded. ",
            "Please copy all text into a newly created txt file and re-import.")
  }
  hub_vec <- utils::read.table(hub_path, sep = "\t", stringsAsFactors = FALSE)[, 1]
  com_genes <- intersect(unique(common_across), unique(hub_vec))
  if (!length(com_genes)) stop("No common hub genes found across datasets.")
  
  # filter and transpose to (samples x genes)
  exp_list <- lapply(exp_list, function(x) t(x[com_genes, , drop = FALSE]))
  
  # read labels and standardize to 0/1
  # read labels and standardize to 0/1 (robust to factor/character/whitespace)
  labels_list <- lapply(file.path(dir, labels_files), .read_labels)
  all_labels <- do.call(rbind, labels_list)
  if (ncol(all_labels) < 2) stop("Label files must have at least two columns: sample_id, label_raw")
  colnames(all_labels)[1:2] <- c("sample_id", "label_raw")
  rownames(all_labels) <- all_labels$sample_id
  
  # 统一成字符并去空格，正类标签也统一成字符比较
  lab_raw_chr <- trimws(as.character(all_labels$label_raw))
  pos_chr    <- trimws(as.character(positive_label))
  
  # 真·0/1，整数型，不走因子
  all_labels$label <- as.integer(lab_raw_chr == pos_chr)
  
  # train set
  x_train <- exp_list[[1]]
  common_ids <- intersect(rownames(x_train), rownames(all_labels))
  if (!length(common_ids)) stop("No overlapping samples between train expression and labels.")
  x_train <- x_train[common_ids, , drop = FALSE]
  
  # y：直接用整数 0/1；若下游某些绘图需要因子，也额外给一个 factor 版本
  y_train_num    <- all_labels[common_ids, "label"]
  y_train_factor <- factor(y_train_num, levels = c(0, 1), labels = c("0","1"))
  y_train <- y_train_factor  # 下游建模用这个
  
  # external tests
  x_tests <- if (length(exp_list) > 1) exp_list[-1] else list()
  if (length(x_tests)) names(x_tests) <- exp_files[-1]
  
  # bundle results
  res <- list(
    x_train = x_train,
    x_tests = x_tests,
    genes = com_genes,
    all_labels = all_labels,
    files = list(allfiles = allfiles, labels_files = labels_files, exp_files = exp_files),
    exp_list = exp_list,
    com_genes = com_genes,
    y_train = y_train_factor,
    pos_rate_summary = {
      # 汇总每个数据集的阳性率，方便 sanity check
      sets <- c("Set_1", if (length(exp_list) > 1) paste0("Set_", seq_len(length(exp_list))[-1]) else character(0))
      lbl <- all_labels[colnames(t(exp_list[[1]])), , drop=FALSE]  # 先用训练集名作为占位
      data.frame(
        Dataset  = c("Set_1"),
        Total    = length(y_train_num),
        Positive = sum(y_train_num == 1L),
        Negative = sum(y_train_num == 0L),
        PosRate  = round(mean(y_train_num == 1L), 4)
      )
    }
  )
  
  # optional: assign into a target environment (e.g., .GlobalEnv)
  if (isTRUE(assign_to_global)) {
    assign("allfiles",     allfiles,     envir = target_env)
    assign("labels_files", labels_files, envir = target_env)
    assign("labels_list",  labels_list,  envir = target_env)
    assign("exp_files",    exp_files,    envir = target_env)
    assign("exp_list",     exp_list,     envir = target_env)
    assign("com_genes",    com_genes,    envir = target_env)
    assign("all_labels",   all_labels,   envir = target_env)
    assign("x_train",      x_train,      envir = target_env)
    assign("y_train",      y_train,      envir = target_env)
    assign("x_tests",      x_tests,      envir = target_env)
  }
  
  invisible(res)
}