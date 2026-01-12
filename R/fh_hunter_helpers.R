#' Internal helpers for FeatureHunter
#' @keywords internal

.ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

.zcol <- function(x) {
    sdx <- stats::sd(x, na.rm = TRUE)
    if (is.na(sdx) || sdx == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / sdx
}

.zrob <- function(x) {
    m <- stats::median(x, na.rm = TRUE)
    s <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- 1
    (x - m) / s
}

.num1 <- function(pat, txt) {
    m <- regmatches(txt, regexpr(pat, txt, perl = TRUE))
    if (length(m) == 0) {
        return(NA_real_)
    }
    as.numeric(gsub("[^0-9\\.]", "", m))
}

.safe_set_names <- function(v, tgt_names) {
    v <- as.numeric(v)
    if (length(v) != length(tgt_names)) {
        # 截断到共同长度，避免长度不一致报错
        k <- min(length(v), length(tgt_names))
        v <- v[seq_len(k)]
        tgt_names <- tgt_names[seq_len(k)]
    }
    names(v) <- tgt_names
    v
}

.cutoff_from_modelstr <- function(txt) {
    mm <- regexec("cutoff:auto\\(([0-9\\.]+),", txt)
    hit <- regmatches(txt, mm)[[1]]
    if (length(hit) >= 2) as.numeric(hit[2]) else .num1("cutoff\\s*:\\s*[0-9\\.]+", txt)
}

.boot_ci <- function(v, R = 1000, alpha = 0.05, seed = 424, progress = FALSE, pb_style = 3) {
    v <- as.numeric(v)
    v <- v[is.finite(v)]
    n <- length(v)

    if (n <= 1L) {
        m <- if (n == 0L) NA_real_ else mean(v, na.rm = TRUE)
        return(c(mean = m, lo = m, hi = m))
    }

    set.seed(seed)

    bstat <- numeric(R)
    pb <- NULL
    if (isTRUE(progress)) {
        pb <- utils::txtProgressBar(min = 0, max = R, style = pb_style)
    }

    for (i in seq_len(R)) {
        idx <- sample.int(n, n, replace = TRUE)
        bstat[i] <- mean(v[idx], na.rm = TRUE)
        if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    if (!is.null(pb)) close(pb)

    m <- mean(v, na.rm = TRUE)
    lo <- as.numeric(stats::quantile(bstat, alpha / 2, na.rm = TRUE, names = FALSE))
    hi <- as.numeric(stats::quantile(bstat, 1 - alpha / 2, na.rm = TRUE, names = FALSE))

    c(mean = m, lo = lo, hi = hi)
}

.find_knee <- function(x) {
    n <- length(x)
    if (n < 3) {
        return(n)
    }
    i <- 2:(n - 1)
    d1 <- x[i - 1] - x[i]
    d2 <- d1[-1] - d1[-length(d1)]
    j <- which.max(pmax(d2, 0))
    if (!length(j)) {
        return(n)
    }
    i[j + 1]
}

.detect_model_type <- function(s) {
    s_low <- tolower(s)
    first_pos <- function(pat) {
        m <- regexpr(pat, s_low, perl = TRUE, ignore.case = TRUE)
        if (m[1] == -1) Inf else as.integer(m[1])
    }
    primary <- list(
        svm = "\\bsvm\\b",
        rf = "(\\brf\\b|random\\s*forest|\\branger\\b)",
        xgb = "(\\bxgboost\\b|\\bxgb\\b)",
        gbm = "\\bgbm\\b",
        mlp = "(\\bnn-mlp\\b|\\bmlp\\b|keras|tensorflow)",
        lda = "\\blda\\b",
        qda = "\\bqda\\b",
        nb = "(naive\\s*bayes|\\bnaivebayes\\b|\\bnb\\b)",
        knn = "\\bknn\\b",
        dt = "(\\bdt\\b|decision\\s*tree|\\btree\\b)",
        lr = "(^|[^a-z])lr([^a-z]|$)|\\blogistic\\b",
        stepwise_lr = "(stepwise.*lr|step\\s*wise.*lr|aic\\+lr)"
    )
    secondary_glmnet <- "(\\benr\\b|elastic\\s*net|\\brr\\b|\\bridge\\b|lasso\\.?r?\\b|\\bglmnet\\b)"
    pos_primary <- vapply(names(primary), function(k) first_pos(primary[[k]]), numeric(1))
    if (!all(is.infinite(pos_primary))) {
        typ <- names(pos_primary)[which.min(pos_primary)]
        if (typ == "svm") {
            if (grepl("kernel\\s*:\\s*linear|linear\\s*kernel", s_low, perl = TRUE)) {
                return("svm_linear")
            } else {
                return("svm_rbf")
            }
        }
        return(typ)
    }
    if (regexpr(secondary_glmnet, s_low, perl = TRUE)[1] != -1) {
        return("glmnet")
    }
    "unknown"
}
