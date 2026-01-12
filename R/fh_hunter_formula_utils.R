#' Internal formula utilities for FeatureHunter
#' @keywords internal
#' @importFrom stats predict coef median mad quantile reorder sd glm cov var
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv

.build_linear_formula <- function(intercept, coefs, names_vec, prefix = "logit(p) = ") {
    fs <- paste0(prefix, round(intercept, 6))
    for (i in seq_along(coefs)) {
        sgn <- ifelse(coefs[i] >= 0, " + ", " - ")
        fs <- paste0(fs, sgn, abs(round(coefs[i], 6)), " * ", names_vec[i])
    }
    fs
}

.build_quadratic_formula <- function(intercept, lin_coefs, quad_mat, names_vec,
                                     prefix = "logit(p) = ", tol = 1e-10) {
    fs <- paste0(prefix, round(intercept, 6))
    if (length(lin_coefs)) {
        for (i in seq_along(lin_coefs)) {
            if (!is.finite(lin_coefs[i]) || abs(lin_coefs[i]) < tol) next
            sgn <- ifelse(lin_coefs[i] >= 0, " + ", " - ")
            fs <- paste0(fs, sgn, abs(round(lin_coefs[i], 6)), " * ", names_vec[i])
        }
    }
    if (!is.null(quad_mat) && all(dim(quad_mat) == c(length(names_vec), length(names_vec)))) {
        p <- length(names_vec)
        for (i in seq_len(p)) {
            qii <- quad_mat[i, i]
            if (is.finite(qii) && abs(qii) >= tol) {
                sgn <- ifelse(qii >= 0, " + ", " - ")
                fs <- paste0(fs, sgn, abs(round(qii, 6)), " * ", names_vec[i], "^2")
            }
        }
        for (i in seq_len(p)) {
            for (j in seq((i + 1), p)) {
                qij <- quad_mat[i, j] + quad_mat[j, i]
                if (is.finite(qij) && abs(qij) >= tol) {
                    sgn <- ifelse(qij >= 0, " + ", " - ")
                    fs <- paste0(fs, sgn, abs(round(qij, 6)), " * ", names_vec[i], " * ", names_vec[j])
                }
            }
        }
    }
    fs
}

.fit_logit_or_ridge <- function(x, y,
                                alpha_ridge = 0, # ridge = 0
                                firth_ci = FALSE,
                                firth_maxit = 1000,
                                firth_pl_maxit = 1000,
                                # FIX: 可调阈值
                                coef_guard = 20, # |beta| 超过就触发回退
                                pred_guard_pct = 0.99, # 极端预测比例阈值
                                pred_eps = 1e-6) {
    x <- as.matrix(x)
    if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))

    # FIX: 规范 y 为 {0,1}
    if (is.factor(y)) {
        lev <- levels(y)
        if (all(lev %in% c("0", "1"))) {
            y <- as.integer(as.character(y)) # "0"/"1" -> 0/1
        } else {
            y <- as.integer(y == lev[2L]) # 以第二个水平为阳性 -> 0/1
        }
    } else {
        y <- as.integer(y)
        if (!all(y %in% c(0L, 1L))) {
            y <- as.integer(y == max(y, na.rm = TRUE))
        }
    }

    # 清理非有限值行
    ok_row <- rowSums(!is.finite(x)) == 0 & is.finite(y)
    if (!all(ok_row)) {
        x <- x[ok_row, , drop = FALSE]
        y <- y[ok_row]
    }

    # FIX: 去除近零方差列（z 后仍可能有 ~0 方差）
    sds <- apply(x, 2, function(col) {
        s <- stats::sd(col, na.rm = TRUE)
        if (!is.finite(s)) 0 else s
    })
    nzv <- which(apply(x, 2, function(col) is.finite(sd(col)) && sd(col) < 1e-6))
    if (length(nzv) > 0) x <- x[, -nzv, drop = FALSE]

    logistic_data <- data.frame(y = y, x)
    sep_flag <- FALSE
    nonconv_flag <- FALSE

    glm_fit <- withCallingHandlers(
        stats::glm(
            y ~ .,
            data = logistic_data,
            family = stats::binomial(link = "logit"),
            control = stats::glm.control(maxit = 200, epsilon = 1e-8)
        ),
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("fitted probabilities numerically 0 or 1 occurred", msg)) {
                sep_flag <<- TRUE
                invokeRestart("muffleWarning")
            } else if (grepl("algorithm did not converge", msg)) {
                nonconv_flag <<- TRUE
                invokeRestart("muffleWarning")
            }
        }
    )

    coef_vec <- tryCatch(stats::coef(glm_fit), error = function(e) rep(NA_real_, ncol(x) + 1))
    bad_coef <- any(!is.finite(coef_vec))

    # ---- 定义一个统一的回退函数 ----
    .fallback_ridge <- function() {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
            return(NULL)
        }
        cvfit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = alpha_ridge, standardize = TRUE)
        beta <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
        intercept <- as.numeric(beta[1, 1])
        coefs <- as.numeric(beta[-1, 1])
        names(coefs) <- rownames(beta)[-1]
        coefs <- coefs[colnames(x)]
        names(coefs) <- colnames(x)
        coef_summary <- cbind(
            Estimate     = c(intercept, coefs),
            `Std. Error` = NA_real_,
            `z value`    = NA_real_,
            `Pr(>|z|)`   = NA_real_
        )
        rownames(coef_summary) <- c("(Intercept)", names(coefs))
        list(
            intercept = intercept,
            coefs = coefs,
            coef_summary = coef_summary,
            source = "ridge_logit_fallback"
        )
    }

    .fallback_firth <- function() {
        if (!requireNamespace("logistf", quietly = TRUE)) {
            return(NULL)
        }
        lf <- tryCatch(
            logistf::logistf(
                y ~ .,
                data = logistic_data,
                control = logistf::logistf.control(maxit = firth_maxit),
                pl = isTRUE(firth_ci),
                plcontrol = if (isTRUE(firth_ci)) logistf::logistpl.control(maxit = firth_pl_maxit) else NULL
            ),
            error = function(e) NULL
        )
        if (is.null(lf)) {
            return(NULL)
        }
        cf <- tryCatch(stats::coef(lf), error = function(e) NULL)
        if (is.null(cf)) cf <- tryCatch(lf$coefficients, error = function(e) NULL)
        if (is.null(cf)) {
            return(NULL)
        }
        varmat <- tryCatch(lf$var, error = function(e) NULL)
        se <- if (!is.null(varmat) && is.matrix(varmat) &&
            nrow(varmat) == length(cf) && ncol(varmat) == length(cf)) {
            sqrt(pmax(diag(varmat), 0))
        } else {
            rep(NA_real_, length(cf))
        }
        zval <- ifelse(is.finite(se) & se > 0, cf / se, NA_real_)
        pval <- ifelse(is.finite(zval), 2 * stats::pnorm(-abs(zval)), NA_real_)
        coef_summary <- cbind(
            Estimate     = as.numeric(cf),
            `Std. Error` = as.numeric(se),
            `z value`    = as.numeric(zval),
            `Pr(>|z|)`   = as.numeric(pval)
        )
        rownames(coef_summary) <- names(cf)

        intercept <- unname(cf["(Intercept)"])
        if (!is.finite(intercept)) intercept <- 0
        coef_vec_only <- cf[setdiff(names(cf), "(Intercept)")]
        coef_vec_only <- coef_vec_only[colnames(x)]
        names(coef_vec_only) <- colnames(x)

        list(
            intercept = as.numeric(intercept),
            coefs = as.numeric(coef_vec_only),
            coef_summary = coef_summary,
            source = "firth_logistic"
        )
    }

    # 先看有没有显式报错/未收敛
    if (sep_flag || nonconv_flag || bad_coef) {
        firth_res <- .fallback_firth()
        if (!is.null(firth_res)) {
            return(firth_res)
        }
        ridge_res <- .fallback_ridge()
        if (!is.null(ridge_res)) {
            return(ridge_res)
        }
        # 实在不行才兜底返回 glm 结果
        intercept <- coef_vec[1]
        coefs <- coef_vec[-1]
        coefs[!is.finite(coefs)] <- 0
        names(coefs) <- colnames(x)
        coef_summary <- tryCatch(summary(glm_fit)$coefficients, error = function(e) {
            cm <- cbind(
                Estimate     = c(intercept, coefs),
                `Std. Error` = NA_real_,
                `z value`    = NA_real_,
                `Pr(>|z|)`   = NA_real_
            )
            rownames(cm) <- c("(Intercept)", colnames(x))
            cm
        })
        return(list(
            intercept = intercept,
            coefs = coefs,
            coef_summary = coef_summary,
            source = "glm_fallback_nonfinite"
        ))
    }

    # FIX: 数值卫兵（即使没有 warning 也检测异常）
    intercept <- coef_vec[1]
    coefs <- coef_vec[-1]
    names(coefs) <- colnames(x)
    p_hat <- tryCatch(stats::predict(glm_fit, type = "response"), error = function(e) NULL)
    extreme_pred <- !is.null(p_hat) && mean(p_hat < pred_eps | p_hat > 1 - pred_eps) > pred_guard_pct
    too_large <- (is.finite(intercept) && abs(intercept) > coef_guard) ||
        any(is.finite(coefs) & abs(coefs) > coef_guard)

    if (extreme_pred || too_large) {
        firth_res <- .fallback_firth()
        if (!is.null(firth_res)) {
            return(firth_res)
        }
        ridge_res <- .fallback_ridge()
        if (!is.null(ridge_res)) {
            return(ridge_res)
        }
        # 兜底仍然返回 glm，但标记来源
        coef_summary <- summary(glm_fit)$coefficients
        return(list(
            intercept = intercept,
            coefs = coefs,
            coef_summary = coef_summary,
            source = "native_glm_guard_triggered_no_fallback"
        ))
    }

    # 正常返回 glm 结果
    coef_summary <- summary(glm_fit)$coefficients
    list(
        intercept = intercept,
        coefs = coefs,
        coef_summary = coef_summary,
        source = "native_glm"
    )
}

.qda_to_logit <- function(X, y) {
    X <- as.matrix(X)
    y <- .y_to_fac(y)
    cls <- levels(y)
    if (length(cls) != 2) stop("QDA requires binary y.")
    x0 <- X[y == cls[1], , drop = FALSE]
    x1 <- X[y == cls[2], , drop = FALSE]
    mu0 <- colMeans(x0)
    mu1 <- colMeans(x1)
    S0 <- stats::cov(x0)
    S1 <- stats::cov(x1)
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    pi0 <- n0 / (n0 + n1)
    pi1 <- 1 - pi0

    # regularize if needed
    reg <- function(S) {
        if (!all(is.finite(S)) || det(S) == 0) S <- S + diag(1e-6, ncol(S))
        S
    }
    S0 <- reg(S0)
    S1 <- reg(S1)

    iS0 <- tryCatch(solve(S0), error = function(e) MASS::ginv(S0))
    iS1 <- tryCatch(solve(S1), error = function(e) MASS::ginv(S1))

    # logit = c + b^T x + x^T A x
    A <- 0.5 * (iS0 - iS1)
    b <- as.numeric(iS1 %*% mu1 - iS0 %*% mu0)
    c <- -0.5 * (as.numeric(crossprod(mu1, iS1 %*% mu1)) -
        as.numeric(crossprod(mu0, iS0 %*% mu0))) +
        0.5 * as.numeric(determinant(S0, logarithm = TRUE)$modulus -
            determinant(S1, logarithm = TRUE)$modulus) +
        log(pi1 / pi0)

    list(intercept = c, linear = b, quad = A)
}

.nb_gaussian_to_logit <- function(X, y) {
    X <- as.matrix(X)
    y <- .y_to_fac(y)
    cls <- levels(y)
    if (length(cls) != 2) stop("NB requires binary y.")
    x0 <- X[y == cls[1], , drop = FALSE]
    x1 <- X[y == cls[2], , drop = FALSE]
    mu0 <- colMeans(x0)
    mu1 <- colMeans(x1)
    v0 <- apply(x0, 2, stats::var)
    v1 <- apply(x1, 2, stats::var)
    v0[!is.finite(v0) | v0 <= 0] <- 1e-6
    v1[!is.finite(v1) | v1 <= 0] <- 1e-6
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    pi0 <- n0 / (n0 + n1)
    pi1 <- 1 - pi0

    # logit = c + sum(b_i x_i) + sum(a_i x_i^2), a_i only on diagonal (no cross terms)
    a <- 0.5 * (1 / v0 - 1 / v1)
    b <- (mu1 / v1) - (mu0 / v0)
    c <- -0.5 * (sum((mu1^2) / v1) - sum((mu0^2) / v0)) -
        0.5 * sum(log(v1 / v0)) + log(pi1 / pi0)

    list(intercept = c, linear = b, quad_diag = a)
}

.generate_formula_report <- function(model_type, x_all, y_all, parsed, svm_cost = 1) {
    formula_source <- NULL
    formula_str <- NULL
    coef_summary <- NULL

    if (identical(model_type, "glmnet")) {
        a_use <- parsed$alpha
        if (is.null(a_use) || is.na(a_use)) a_use <- 1
        y_all_num <- .y_to_num(y_all)
        cvfit <- glmnet::cv.glmnet(x_all, y_all_num,
            family = "binomial",
            alpha = a_use, standardize = TRUE
        )
        beta <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
        intercept <- as.numeric(beta[1, 1])

        # 先用beta自带的行名（去掉截距），再与x_all列名对齐，避免长度不一致
        raw_names <- rownames(beta)[-1]
        raw_coefs <- as.numeric(beta[-1, 1])
        names(raw_coefs) <- raw_names

        # 按 x_all 的列顺序对齐（有些列可能被惩罚到0，但仍保留在coef里）
        aligned <- setNames(rep(0, ncol(x_all)), colnames(x_all))
        common <- intersect(names(raw_coefs), colnames(x_all))
        aligned[common] <- raw_coefs[common]
        coefs <- .safe_set_names(aligned, colnames(x_all))

        coef_summary <- cbind(
            Estimate     = c(intercept, as.numeric(coefs)),
            `Std. Error` = NA_real_,
            `z value`    = NA_real_,
            `Pr(>|z|)`   = NA_real_
        )
        rownames(coef_summary) <- c("(Intercept)", colnames(x_all))

        formula_str <- .build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
        cat(formula_str, "\n")
        formula_source <- "native_glmnet"
    } else if (identical(model_type, "lr") || identical(model_type, "stepwise_lr")) {
        # native logistic on selected features (with robust fallback)
        res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
        intercept <- res$intercept
        coefs <- res$coefs
        coef_summary <- res$coef_summary
        feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
        names(coefs) <- feat_names_out
        rownames(coef_summary) <- c("(Intercept)", names(coefs))
        formula_str <- .build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
        cat(formula_str, "\n")
        formula_source <- paste0(model_type, "_native_logit")
    } else if (identical(model_type, "lda")) {
        if (length(unique(y_all)) != 2) {
            # fallback: 用逻辑回归/岭回归
            res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
            intercept <- res$intercept
            coefs <- .safe_set_names(res$coefs, colnames(x_all))
            coef_summary <- cbind(
                Estimate     = c(intercept, as.numeric(coefs)),
                `Std. Error` = NA_real_,
                `z value`    = NA_real_,
                `Pr(>|z|)`   = NA_real_
            )
            rownames(coef_summary) <- c("(Intercept)", colnames(x_all))
            formula_str <- .build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
            cat(formula_str, "\n")
            formula_source <- paste0("lda_", res$source)
        } else {
            # 原生 LDA → 等价的logit系数
            cls <- sort(unique(y_all))
            x0 <- x_all[y_all == cls[1], , drop = FALSE]
            x1 <- x_all[y_all == cls[2], , drop = FALSE]
            mu0 <- colMeans(x0)
            mu1 <- colMeans(x1)
            S0 <- stats::cov(x0)
            S1 <- stats::cov(x1)
            n0 <- nrow(x0)
            n1 <- nrow(x1)
            Sp <- ((n0 - 1) * S0 + (n1 - 1) * S1) / (n0 + n1 - 2)
            if (any(!is.finite(Sp)) || det(Sp) == 0) {
                Sp <- Sp + diag(1e-6, ncol(Sp))
            }
            beta_vec <- tryCatch(solve(Sp, (mu1 - mu0)),
                error = function(e) {
                    MASS::ginv(Sp) %*% (mu1 - mu0)
                }
            )
            beta_vec <- as.numeric(beta_vec)
            coefs <- .safe_set_names(beta_vec, colnames(x_all))

            pi1 <- n1 / (n0 + n1)
            pi0 <- 1 - pi1
            intercept <- as.numeric(-0.5 * crossprod((mu1 + mu0), as.numeric(coefs)) + log(pi1 / pi0))

            coef_summary <- cbind(
                Estimate     = c(intercept, as.numeric(coefs)),
                `Std. Error` = NA_real_,
                `z value`    = NA_real_,
                `Pr(>|z|)`   = NA_real_
            )
            rownames(coef_summary) <- c("(Intercept)", colnames(x_all))

            formula_str <- .build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
            cat(formula_str, "\n")
            formula_source <- "native_lda_logit"
        }
    } else if (identical(model_type, "qda")) {
        q <- .qda_to_logit(x_all, y_all)
        intercept <- q$intercept
        lin <- q$linear
        names(lin) <- colnames(x_all)
        quad <- q$quad
        rownames(quad) <- colnames(x_all)
        colnames(quad) <- colnames(x_all)
        coef_summary <- cbind(
            Estimate = c(intercept, lin),
            `Std. Error` = NA_real_,
            `z value` = NA_real_,
            `Pr(>|z|)` = NA_real_
        )
        rownames(coef_summary) <- c("(Intercept)", names(lin))
        formula_str <- .build_quadratic_formula(intercept, lin, quad, colnames(x_all), prefix = "logit(p) = ")
        cat(formula_str, "\n")
        formula_source <- "native_qda_quadratic"
    } else if (identical(model_type, "nb")) {
        nbp <- .nb_gaussian_to_logit(x_all, y_all)
        intercept <- nbp$intercept
        lin <- nbp$linear
        names(lin) <- colnames(x_all)
        quad_diag <- nbp$quad_diag
        Q <- diag(as.numeric(quad_diag), nrow = ncol(x_all))
        rownames(Q) <- colnames(x_all)
        colnames(Q) <- colnames(x_all)
        coef_summary <- cbind(
            Estimate = c(intercept, lin),
            `Std. Error` = NA_real_,
            `z value` = NA_real_,
            `Pr(>|z|)` = NA_real_
        )
        rownames(coef_summary) <- c("(Intercept)", names(lin))
        formula_str <- .build_quadratic_formula(intercept, lin, Q, colnames(x_all), prefix = "logit(p) = ")
        cat(formula_str, "\n")
        formula_source <- "native_nb_gaussian_quadratic"
    } else if (identical(model_type, "svm_linear")) {
        svm_fit <- e1071::svm(
            x = x_all, y = factor(y_all),
            kernel = "linear", probability = TRUE,
            cost = svm_cost, scale = TRUE
        )
        if (is.null(svm_fit$coefs) || is.null(svm_fit$SV)) {
            res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
            intercept <- res$intercept
            coefs <- res$coefs
            coef_summary <- res$coef_summary
            feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
            names(coefs) <- feat_names_out
            formula_str <- .build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
            cat(formula_str, "\n")
            formula_source <- "svm_linear_fallback_logit"
        } else {
            wv <- as.numeric(crossprod(svm_fit$coefs, svm_fit$SV))
            b <- -as.numeric(svm_fit$rho)
            names(wv) <- colnames(x_all)
            margin_str <- .build_linear_formula(b, wv, names(wv), prefix = "margin = ")
            cat(margin_str, "\n")
            coef_summary <- cbind(
                Estimate = c(b, wv),
                `Std. Error` = NA_real_,
                `z value` = NA_real_,
                `Pr(>|z|)` = NA_real_
            )
            rownames(coef_summary) <- c("(Intercept)", names(wv))
            formula_str <- margin_str
            formula_source <- "native_svm_linear_margin"
        }
    } else if (model_type %in% c("mlp", "rf", "xgb", "gbm", "svm_rbf", "dt", "tree", "knn")) {
        # black-box models: interpretable linear surrogate
        res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
        intercept <- res$intercept
        coefs <- res$coefs
        coef_summary <- res$coef_summary
        feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
        names(coefs) <- feat_names_out
        formula_str <- .build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
        cat(formula_str, "\n")
        formula_source <- paste0(model_type, "_surrogate_", res$source)
    } else {
        # unknown: fallback to interpretable logistic
        res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
        intercept <- res$intercept
        coefs <- res$coefs
        coef_summary <- res$coef_summary
        feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
        names(coefs) <- feat_names_out
        formula_str <- .build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
        cat(formula_str, "\n")
        formula_source <- "unknown_surrogate_logit"
    }

    list(
        formula_source = formula_source,
        formula_str = formula_str,
        coef_summary = coef_summary
    )
}
