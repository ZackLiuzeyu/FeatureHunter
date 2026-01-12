# 看每个 label 表的“样子”：维度、前几行、唯一取值、NA 数
lapply(seq_along(labels_list), function(i){
  df <- labels_list[[i]]
  list(
    idx = i,
    dim = dim(df),
    head = head(df, 3),
    uniques = unique(df[[2]]),
    n_na = sum(is.na(df[[2]]))
  )
})

check_align <- function(i){
  Xi <- as.data.frame(exp_list[[i]])      # i=1 是 A, i=2 是 B ...
  Li <- labels_list[[i]]
  ids <- as.character(Li[[1]])
  c(
    rows_X = nrow(Xi),
    rows_L = nrow(Li),
    intersect_n = length(intersect(rownames(Xi), ids)),
    missing_in_label = sum(!(rownames(Xi) %in% ids))
  )
}
sapply(seq_along(exp_list), check_align)

.clean_str <- function(x){
  x <- gsub("[\u00A0\u2000-\u200B\uFEFF]", "", x)  # 不可见空白
  trimws(tolower(x))
}
normalize_labels <- function(df, positive_alias = c("1","disease","tumor","case","positive","pos","yes","true")){
  stopifnot(ncol(df) >= 2)
  df[[1]] <- .clean_str(as.character(df[[1]]))
  lab_raw <- .clean_str(as.character(df[[2]]))
  lab01 <- ifelse(lab_raw %in% positive_alias, 1L,
                  ifelse(lab_raw %in% c("0","control","normal","negative","neg","no","false"), 0L, NA_integer_))
  data.frame(df[[1]], lab01, stringsAsFactors = FALSE,
             row.names = NULL, check.names = FALSE)
}
labels_list_fixed <- lapply(labels_list, normalize_labels)

pos_rate_summary2 <- do.call(rbind, lapply(seq_along(labels_list_fixed), function(i){
  y <- labels_list_fixed[[i]][[2]]
  cbind(
    Dataset = paste0("Set_", i),
    Total   = length(y),
    Positive= sum(y == 1, na.rm=TRUE),
    Negative= sum(y == 0, na.rm=TRUE),
    PosRate = mean(y == 1, na.rm=TRUE)
  )
}))
pos_rate_summary2
