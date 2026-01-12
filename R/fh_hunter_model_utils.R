#' Internal model utilities for FeatureHunter
#' @keywords internal
#' @importFrom stats predict coef median mad quantile reorder sd glm
#' @importFrom utils txtProgressBar setTxtProgressBar

.uses_shap_for_type <- function(model_type) {
    # parametric, intrinsically interpretable models: skip SHAP
    skip_set <- c("glmnet", "lr", "stepwise_lr", "svm_linear", "lda", "qda", "nb")
    !(model_type %in% skip_set)
}

.pred_proba <- function(model, X, model_type) {
    X <- as.matrix(X)

    .pos_ix <- function(pr, prefer = c("1", "pos", "positive", "TRUE", "Case")) {
        ix <- integer(0)
        if (!is.null(colnames(pr))) ix <- which(colnames(pr) %in% prefer)
        if (!length(ix)) ix <- ncol(pr)
        ix[1]
    }
    .clamp01 <- function(v) {
        v <- as.numeric(v)
        v[!is.finite(v)] <- NA_real_
        pmin(pmax(v, 0), 1)
    }

    switch(model_type,

        ## -------------------- Deep / Linear families --------------------
        mlp = {
            .clamp01(stats::predict(model, X, verbose = 0))
        },
        glmnet = {
            if (inherits(model, "cv.glmnet")) {
                .clamp01(as.numeric(predict(model, newx = X, s = "lambda.min", type = "response")))
            } else {
                .clamp01(as.numeric(predict(model, newx = X, type = "response")))
            }
        },
        lr = {
            .clamp01(stats::predict(model, newdata = data.frame(X), type = "response"))
        },
        stepwise_lr = {
            .clamp01(stats::predict(model, newdata = data.frame(X), type = "response"))
        },
        lda = {
            pr <- predict(model, data.frame(X))$posterior
            .clamp01(pr[, .pos_ix(pr)])
        },
        qda = {
            pr <- predict(model, data.frame(X))$posterior
            .clamp01(pr[, .pos_ix(pr)])
        },
        nb = {
            pr <- predict(model, data.frame(X), type = "raw")
            .clamp01(pr[, .pos_ix(pr)])
        },

        ## -------------------- Tree / Ensemble families --------------------
        rf = {
            pr <- predict(model, data = data.frame(X), type = "response")$predictions
            if (is.matrix(pr)) .clamp01(pr[, .pos_ix(pr)]) else .clamp01(pr)
        },
        xgb = {
            # predict.xgb.Booster
            .clamp01(stats::predict(model, newdata = X))
        },
        gbm = {
            ntrees <- if (!is.null(model$n.trees)) {
                model$n.trees
            } else if (!is.null(model$trees)) {
                length(model$trees)
            } else {
                100L
            }
            .clamp01(predict(model, newdata = data.frame(X), n.trees = ntrees, type = "response"))
        },
        dt = {
            if (inherits(model, "rpart")) {
                pr <- predict(model, newdata = data.frame(X), type = "prob")
                if (is.null(dim(pr))) pr <- cbind(`0` = 1 - pr, `1` = pr)
                .clamp01(pr[, .pos_ix(pr)])
            } else if (inherits(model, "party") || inherits(model, "constparty")) {
                # partykit
                if (requireNamespace("partykit", quietly = TRUE)) {
                    pr <- stats::predict(model, newdata = data.frame(X), type = "prob")
                    if (is.list(pr)) {
                        cls <- names(pr[[1]])
                        M <- matrix(NA_real_, nrow = length(pr), ncol = length(cls))
                        colnames(M) <- cls
                        for (i in seq_along(pr)) M[i, ] <- pr[[i]]
                        pr <- M
                    }
                    .clamp01(pr[, .pos_ix(pr)])
                } else {
                    stop("[DT] partykit not installed for this tree model.")
                }
            } else {
                pr <- tryCatch(predict(model, newdata = data.frame(X), type = "prob"),
                    error = function(e) NULL
                )
                if (is.null(pr)) stop("[DT] unsupported tree object for probability prediction.")
                if (is.null(dim(pr))) pr <- cbind(`0` = 1 - pr, `1` = pr)
                .clamp01(pr[, .pos_ix(pr)])
            }
        },
        tree = { # alias of dt
            if (inherits(model, "rpart")) {
                pr <- predict(model, newdata = data.frame(X), type = "prob")
                if (is.null(dim(pr))) pr <- cbind(`0` = 1 - pr, `1` = pr)
                .clamp01(pr[, .pos_ix(pr)])
            } else {
                pr <- tryCatch(predict(model, newdata = data.frame(X), type = "prob"),
                    error = function(e) NULL
                )
                if (is.null(pr)) stop("[tree] unsupported tree object for probability prediction.")
                if (is.null(dim(pr))) pr <- cbind(`0` = 1 - pr, `1` = pr)
                .clamp01(pr[, .pos_ix(pr)])
            }
        },

        ## -------------------- SVM families --------------------
        svm_linear = {
            pr <- attr(predict(model, data.frame(X), probability = TRUE), "probabilities")
            if (is.null(pr)) stop("[SVM-linear] Train with probability=TRUE")
            if (!is.null(model$levels) && all(model$levels %in% colnames(pr))) {
                pos_lab <- intersect(model$levels, c("1", "pos", "positive", "TRUE", "Case"))
                if (!length(pos_lab)) pos_lab <- tail(model$levels, 1)
                .clamp01(pr[, pos_lab[1]])
            } else {
                .clamp01(pr[, .pos_ix(pr)])
            }
        },
        svm_rbf = {
            pr <- attr(predict(model, data.frame(X), probability = TRUE), "probabilities")
            if (is.null(pr)) stop("[SVM-rbf] Train with probability=TRUE")
            if (!is.null(model$levels) && all(model$levels %in% colnames(pr))) {
                pos_lab <- intersect(model$levels, c("1", "pos", "positive", "TRUE", "Case"))
                if (!length(pos_lab)) pos_lab <- tail(model$levels, 1)
                .clamp01(pr[, pos_lab[1]])
            } else {
                .clamp01(pr[, .pos_ix(pr)])
            }
        },

        ## -------------------- KNN (自定义“模型对象”) --------------------
        knn = {
            # 你的“模型”是 list(kind="knn", X_train, y, k, center, scale)
            if (!is.list(model) || is.null(model$kind) || model$kind != "knn") {
                stop("[knn] Expect a list(kind='knn', ...).")
            }
            Xs <- scale(X, center = model$center, scale = model$scale)
            Xtr <- model$X_train
            k <- as.integer(model$k)
            # 优先用 FNN::get.knnx（快）；否则用简易欧氏距离
            if (requireNamespace("FNN", quietly = TRUE)) {
                nn <- FNN::get.knnx(data = Xtr, query = Xs, k = k)
                idx <- nn$nn.index
            } else {
                # 简单暴力：逐点计算距离（n*m，可能慢，但通用）
                idx <- matrix(NA_integer_, nrow = nrow(Xs), ncol = k)
                for (i in seq_len(nrow(Xs))) {
                    d <- rowSums((t(t(Xtr) - Xs[i, ]))^2)
                    ord <- order(d, decreasing = FALSE)[seq_len(min(k, length(d)))]
                    if (length(ord) < k) ord <- c(ord, rep(ord[length(ord)], k - length(ord)))
                    idx[i, ] <- ord
                }
            }
            yy <- as.integer(model$y)
            # 计算邻居中“正类”比例作为概率
            p_hat <- rowMeans(matrix(yy[idx], nrow = nrow(idx)) == 1, na.rm = TRUE)
            .clamp01(p_hat)
        },

        ## -------------------- default --------------------
        {
            stop("Unsupported model_type in .pred_proba: ", model_type)
        }
    )
}

.pred_proba_ensemble <- function(models, X, model_type) {
    if (length(models) == 1) {
        return(.pred_proba(models[[1]], X, model_type))
    }
    P <- vapply(models, function(m) .pred_proba(m, X, model_type), numeric(nrow(as.matrix(X))))
    rowMeans(cbind(P), na.rm = TRUE)
}

.perm_importance_fun <- function(
    predict_fun, X, y,
    metric = c("f1", "prauc"),
    nrep = 3L,
    stratified = TRUE,
    seed = 424,
    progress = FALSE,
    pb_style = 3) {
    metric <- match.arg(metric)
    set.seed(seed)

    base_p <- predict_fun(X)
    base_s <- switch(metric,
        f1 = {
            pr <- as.integer(base_p >= 0.5)
            tp <- sum(pr == 1 & y == 1)
            fp <- sum(pr == 1 & y == 0)
            fn <- sum(pr == 0 & y == 1)
            if (tp + fp + fn == 0) 0 else (2 * tp) / (2 * tp + fp + fn)
        },
        prauc = PRROC::pr.curve(scores.class0 = base_p[y == 1], scores.class1 = base_p[y == 0])$auc.integral
    )

    p <- ncol(X)
    drop <- numeric(p)
    idx0 <- which(y == 0)
    idx1 <- which(y == 1)

    pb <- NULL
    if (isTRUE(progress)) {
        pb <- utils::txtProgressBar(min = 0, max = p, style = pb_style)
    }

    for (j in seq_len(p)) {
        s <- 0
        for (r in seq_len(nrep)) {
            Xp <- X
            if (isTRUE(stratified)) {
                Xp[idx0, j] <- sample(X[idx0, j])
                Xp[idx1, j] <- sample(X[idx1, j])
            } else {
                Xp[, j] <- sample(X[, j])
            }
            prp <- predict_fun(Xp)
            s <- s + switch(metric,
                f1 = {
                    pr <- as.integer(prp >= 0.5)
                    tp <- sum(pr == 1 & y == 1)
                    fp <- sum(pr == 1 & y == 0)
                    fn <- sum(pr == 0 & y == 1)
                    if (tp + fp + fn == 0) 0 else (2 * tp) / (2 * tp + fp + fn)
                },
                prauc = PRROC::pr.curve(scores.class0 = prp[y == 1], scores.class1 = prp[y == 0])$auc.integral
            )
        }
        drop[j] <- base_s - s / nrep
        if (!is.null(pb)) utils::setTxtProgressBar(pb, j)
    }
    if (!is.null(pb)) close(pb)

    pmax(drop, 0)
}

.inner_importance <- function(model, X, model_type, rbf_grad_subsample = 200L) {
    X <- as.matrix(X)
    p <- ncol(X)
    cn <- colnames(X)

    safe_vec <- function(v, names_ref) {
        v[!is.finite(v)] <- 0
        names(v) <- names_ref
        as.numeric(v)
    }

    out <- switch(model_type,

        # 1) MLP: use the first linear layer "slayer" absolute weights
        mlp = {
            if (!requireNamespace("keras3", quietly = TRUE)) stop("Package 'keras3' needed for MLP importance.")
            w <- tryCatch(keras3::get_weights(keras3::get_layer(model, "slayer"))[[1]], error = function(e) NULL)
            if (is.null(w)) rep(0, p) else safe_vec(abs(as.numeric(w)), cn)
        },

        # 2) Random Forest (ranger): variable.importance
        rf = {
            vi <- model$variable.importance
            if (is.null(vi)) {
                out0 <- rep(0, p)
                names(out0) <- cn
                as.numeric(out0)
            } else {
                v <- rep(0, p)
                names(v) <- cn
                v[names(vi)] <- as.numeric(vi)
                safe_vec(v[cn], cn)
            }
        },

        # 3) GLMNET (cv.glmnet or glmnet): |beta| * sd(X)
        glmnet = {
            # glmnet is loaded
            bvec <- tryCatch(as.numeric(stats::coef(model, s = "lambda.min")), error = function(e) as.numeric(stats::coef(model)))
            if (length(bvec) == p + 1) {
                b <- bvec[-1]
                sx <- apply(X, 2, stats::sd)
                sx[!is.finite(sx) | sx == 0] <- 1
                safe_vec(abs(b) * sx, cn)
            } else {
                rep(0, p)
            }
        },

        # 4) XGBoost: gain-based importance mapped to columns
        xgb = {
            if (!requireNamespace("xgboost", quietly = TRUE)) stop("Package 'xgboost' needed for XGB importance.")
            im <- tryCatch(xgboost::xgb.importance(feature_names = cn, model = model), error = function(e) NULL)
            if (is.null(im) || nrow(im) == 0) {
                rep(0, p)
            } else {
                v <- setNames(im$Gain, im$Feature)
                vv <- v[cn]
                vv[is.na(vv)] <- 0
                safe_vec(as.numeric(vv), cn)
            }
        },

        # 5) GBM (gbm package): relative influence
        gbm = {
            if (!requireNamespace("gbm", quietly = TRUE)) stop("Package 'gbm' needed for GBM importance.")
            ri <- tryCatch(
                {
                    nt <- if (!is.null(model$n.trees)) model$n.trees else NULL
                    if (is.null(nt)) gbm::relative.influence(model, normalized = TRUE) else gbm::relative.influence(model, n.trees = nt, normalized = TRUE)
                },
                error = function(e) NULL
            )
            if (is.null(ri)) {
                rep(0, p)
            } else {
                v <- rep(0, p)
                names(v) <- cn
                # names(ri) are feature names
                v[names(ri)] <- as.numeric(ri)
                safe_vec(v[cn], cn)
            }
        },

        # 6) SVM (linear): |w|
        svm_linear = {
            if (is.null(model$coefs) || is.null(model$SV)) {
                rep(0, p)
            } else {
                w <- as.numeric(crossprod(model$coefs, model$SV))
                safe_vec(abs(w), cn)
            }
        },

        # 7) SVM (RBF): gradient-based proxy averaged over a subsample
        svm_rbf = {
            if (is.null(model$coefs) || is.null(model$SV)) {
                rep(0, p)
            } else {
                gamma <- if (!is.null(model$gamma)) model$gamma else 1 / ncol(X)
                Xs <- X
                if (isTRUE(model$scaled) && !is.null(model$x.scale)) {
                    cen <- model$x.scale$`scaled:center`
                    scl <- model$x.scale$`scaled:scale`
                    scl[is.na(scl) | scl == 0] <- 1
                    Xs <- sweep(sweep(X, 2, cen, "-"), 2, scl, "/")
                }
                n <- nrow(Xs)
                idx <- if (n > rbf_grad_subsample) sample(n, rbf_grad_subsample) else seq_len(n)
                Xsub <- Xs[idx, , drop = FALSE]
                SV <- model$SV
                al <- as.numeric(model$coefs)
                imp_sum <- numeric(p)
                for (k in seq_len(nrow(Xsub))) {
                    xk <- Xsub[k, ]
                    diff <- sweep(SV, 2, xk, "-")
                    d2 <- rowSums(diff * diff)
                    Kvec <- exp(-gamma * d2)
                    wsum <- colSums((al * Kvec) * (SV - matrix(xk, nrow(SV), p, byrow = TRUE)))
                    gk <- (2 * gamma) * wsum
                    imp_sum <- imp_sum + abs(gk)
                }
                safe_vec(imp_sum / max(1L, length(idx)), cn)
            }
        },

        # 8) LDA: |scaling|
        lda = {
            sc <- model$scaling
            if (is.null(sc)) {
                rep(0, p)
            } else {
                v <- rep(0, p)
                names(v) <- cn
                v[rownames(sc)] <- abs(as.numeric(sc[, 1]))
                safe_vec(v[cn], cn)
            }
        },

        # 9) QDA: |mean1 - mean0| (uses model$means)
        qda = {
            means <- tryCatch(model$means, error = function(e) NULL)
            if (is.null(means) || ncol(means) != p) {
                rep(0, p)
            } else {
                if (nrow(means) < 2) {
                    rep(0, p)
                } else {
                    v <- abs(as.numeric(means[2, ] - means[1, ]))
                    names(v) <- colnames(means)
                    vv <- rep(0, p)
                    names(vv) <- cn
                    vv[names(v)] <- v
                    safe_vec(vv[cn], cn)
                }
            }
        },

        # 10) Naive Bayes (e1071): for numeric features use |mu1 - mu0| / pooled_sd if available
        nb = {
            tabs <- tryCatch(model$tables, error = function(e) NULL)
            if (is.null(tabs)) {
                rep(0, p)
            } else {
                v <- rep(0, p)
                names(v) <- cn
                for (nm in intersect(names(tabs), cn)) {
                    tb <- tabs[[nm]]
                    # numeric features have a 2x2 (or 2x?) matrix with mean/sd per class
                    if (is.matrix(tb) && all(c("mean", "sd") %in% tolower(colnames(tb)))) {
                        cnms <- tolower(colnames(tb))
                        mu0 <- as.numeric(tb[1, which(cnms == "mean")[1]])
                        mu1 <- if (nrow(tb) >= 2) as.numeric(tb[2, which(cnms == "mean")[1]]) else NA
                        sd0 <- as.numeric(tb[1, which(cnms == "sd")[1]])
                        sd1 <- if (nrow(tb) >= 2) as.numeric(tb[2, which(cnms == "sd")[1]]) else NA
                        if (is.finite(mu0) && is.finite(mu1) && is.finite(sd0) && is.finite(sd1)) {
                            sp <- sqrt((sd0^2 + sd1^2) / 2)
                            if (!is.finite(sp) || sp == 0) sp <- 1
                            v[nm] <- abs(mu1 - mu0) / sp
                        }
                    } else if (is.matrix(tb) && nrow(tb) >= 2 && ncol(tb) >= 2) {
                        # fallback: absolute mean diff without sd
                        v[nm] <- abs(as.numeric(tb[2, 1] - tb[1, 1]))
                    } else {
                        v[nm] <- 0
                    }
                }
                safe_vec(v[cn], cn)
            }
        },

        # 11) Logistic regression (stats::glm binomial): |beta| * sd(X)
        lr = {
            co <- tryCatch(stats::coef(model), error = function(e) NULL)
            if (is.null(co)) {
                rep(0, p)
            } else {
                b <- co[names(co) %in% paste0(cn)] # match by exact names if possible
                # fallback: drop intercept and assume order
                if (length(b) == 0L && length(co) >= 2L) b <- co[-1]
                b <- as.numeric(b)
                if (length(b) != p) {
                    # align by colnames when possible
                    tmp <- rep(0, p)
                    names(tmp) <- cn
                    nb <- names(co)
                    nb <- nb[nb != "(Intercept)"]
                    tmp[nb] <- as.numeric(co[nb])
                    b <- as.numeric(tmp)
                }
                sx <- apply(X, 2, stats::sd)
                sx[!is.finite(sx) | sx == 0] <- 1
                safe_vec(abs(b) * sx, cn)
            }
        },

        # 12) Stepwise logistic: same as lr
        stepwise_lr = {
            co <- tryCatch(stats::coef(model), error = function(e) NULL)
            if (is.null(co)) {
                rep(0, p)
            } else {
                b <- co[names(co) %in% paste0(cn)]
                if (length(b) == 0L && length(co) >= 2L) b <- co[-1]
                b <- as.numeric(b)
                if (length(b) != p) {
                    tmp <- rep(0, p)
                    names(tmp) <- cn
                    nb <- names(co)
                    nb <- nb[nb != "(Intercept)"]
                    tmp[nb] <- as.numeric(co[nb])
                    b <- as.numeric(tmp)
                }
                sx <- apply(X, 2, stats::sd)
                sx[!is.finite(sx) | sx == 0] <- 1
                safe_vec(abs(b) * sx, cn)
            }
        },

        # 13) Decision tree (rpart/partykit): variable importance
        dt = {
            # rpart
            vi <- tryCatch(model$variable.importance, error = function(e) NULL)
            if (!is.null(vi)) {
                v <- rep(0, p)
                names(v) <- cn
                v[names(vi)] <- as.numeric(vi)
                safe_vec(v[cn], cn)
            } else {
                # partykit varimp if available
                v2 <- tryCatch(
                    {
                        if (requireNamespace("partykit", quietly = TRUE)) {
                            as.numeric(partykit::varimp(model))
                        } else {
                            NULL
                        }
                    },
                    error = function(e) NULL
                )
                if (is.null(v2)) {
                    rep(0, p)
                } else {
                    nv <- partykit::varimp(model)
                    v <- rep(0, p)
                    names(v) <- cn
                    v[names(nv)] <- as.numeric(nv)
                    safe_vec(v[cn], cn)
                }
            }
        },
        tree = { # alias
            vi <- tryCatch(model$variable.importance, error = function(e) NULL)
            if (!is.null(vi)) {
                v <- rep(0, p)
                names(v) <- cn
                v[names(vi)] <- as.numeric(vi)
                safe_vec(v[cn], cn)
            } else {
                rep(0, p)
            }
        },

        # 14) KNN: no intrinsic model-based importance; return zeros (perm/SHAP will cover)
        knn = {
            rep(0, p)
        },

        # default
        {
            stop("Unsupported model_type in .inner_importance: ", model_type)
        }
    )

    names(out) <- cn
    out[!is.finite(out)] <- 0
    as.numeric(out)
}

.shap_importance_fun <- function(
    predict_fun,
    X,
    shap_nsim = 10L,
    shap_subsample = 100L,
    seed = 424,
    progress = TRUE,
    pb_style = 3,
    update_every = 5) {
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    set.seed(seed)
    idx <- if (n > shap_subsample) sample(n, shap_subsample) else seq_len(n)
    Xsub <- X[idx, , drop = FALSE]

    acc <- rep(0, p)
    names(acc) <- colnames(X)

    pb <- NULL
    if (isTRUE(progress)) {
        pb <- utils::txtProgressBar(min = 0, max = shap_nsim, style = pb_style)
    }

    for (i in seq_len(shap_nsim)) {
        sv <- fastshap::explain(
            object = list(),
            X = as.data.frame(Xsub),
            pred_wrapper = function(object, newdata) predict_fun(as.matrix(newdata)),
            nsim = 1L,
            adjust = FALSE
        )
        svm <- colMeans(abs(as.matrix(sv)), na.rm = TRUE)
        acc[names(svm)] <- acc[names(svm)] + svm

        if (!is.null(pb) && (i %% update_every == 0 || i == shap_nsim)) {
            utils::setTxtProgressBar(pb, i)
        }
    }
    if (!is.null(pb)) close(pb)

    imp <- acc / shap_nsim
    imp[match(colnames(X), names(imp))]
}
