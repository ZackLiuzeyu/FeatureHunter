#' Universal Feature Hunter across Multiple ML Models
#'
#' @description
#' Automatically detect the selected model type (MLP / RF / GLMNET / XGBoost / SVM / LDA / QDA / Naive Bayes)
#' from a leaderboard, train the corresponding model, and compute three types of feature importance:
#' \itemize{
#'   \item Model-based importance (weights, coefficients, gain, etc.)
#'   \item Permutation importance (performance drop by shuffling a feature)
#'   \item SHAP importance (model-agnostic explanation)
#' }
#' The three importance scores are fused via robust z-scores (median/MAD) into a composite score to rank features.
#' To reduce variance, multiple seeds are trained, predictions are ensembled first, then permutation/SHAP are computed
#' on the ensemble predictor; error bars use bootstrap CI over per-run composite scores.
#'
#' @param train_exp Numeric matrix of predictors (samples x features).
#' @param train_labels Response vector (0/1, factor, or convertible).
#' @param nshow Number of models shown in leaderboard.
#' @param namesS Character vector of metric names shown in leaderboard
#'   (default: \code{c("Accuracy","Recall","F-score")}).
#' @param score_index Integer index of metric to rank by (default: 3 = F-score).
#' @param pick_index Row index of leaderboard to pick (default: 1). Works for any model type.
#' @param top_models_csv Path to leaderboard CSV (auto-inferred if \code{NULL}).
#' @param num_runs Number of repeat runs (seeds) (default: 10).
#' @param num_coregene Number of core genes to keep in stability/importance plots.
#' @param n_likes Integer or "auto".Number of top genes to use in downstream analysis (UMAP, logistic). When "auto", UMAP/logistic will use the count of genes whose 95% CI does not cross 0 within the composite bar plot candidate set.
#' @param n_interest Number of genes to show in the stability heatmap.
#' @param seed Random seed (default: 424).
#' @param out_dir Output directory for plots (default: working dir).
#'
#' @param learningratei,batch_sizei,epochselecti,dropoutratei,cutoffi
#' Optional overrides for MLP hyperparameters (if \code{NULL}, use internal defaults
#' or values parsed from the leaderboard string where applicable).
#'
#' @param perm_metric Performance metric for permutation importance (\code{"f1"} or \code{"prauc"}).
#' @param perm_nrep Number of repetitions for permutation (default: 3).
#' @param shap_nsim Number of Monte Carlo simulations for SHAP (default: 10).
#' @param shap_subsample Number of samples for SHAP subsampling (default: 100).
#'
#' @param svm_cost Cost parameter for SVM (default: 1).
#' @param svm_gamma Gamma for SVM-RBF (default: \code{1/p} if \code{NULL}).
#'
#' @param xgb_nrounds,xgb_eta,xgb_max_depth Hyperparameters for XGBoost
#'   (defaults: \code{200}, \code{0.05}, \code{4}).
#'
#' @param rf_num_trees,rf_mtry Hyperparameters for Random Forest
#'   (defaults: \code{500} trees, \code{sqrt(p)} mtry if \code{NULL}).
#'
#' @param shuffle Logical; whether to shuffle data each epoch in MLP training (default: \code{TRUE}).
#'
#' @param standardize logical; if TRUE, standardize using train-set mean/sd and apply to val/test.
#'   Default: NULL (auto).
#' @param hidden_units integer vector; number of units per hidden layer,
#'   e.g. c(32L,16L,8L) or c(64L,32L,16L,8L). Default: NULL (auto).
#' @param activation character; activation function for hidden layers,
#'   e.g. "relu". Default: NULL (auto).
#' @param use_batchnorm logical; whether to insert BatchNorm after each hidden Dense.
#'   Default: NULL (auto).
#' @param l2 numeric; L2 regularization strength (e.g. 1e-4). Default: NULL (auto).
#' @param gaussian_noise_sd numeric; stddev for input GaussianNoise layer.
#'   Default: NULL (disabled).
#' @param min_lr numeric; minimum learning rate for ReduceLROnPlateau.
#'   Default: NULL (auto).
#' @param plateau_factor numeric; factor for ReduceLROnPlateau (e.g. 0.5).
#'   Default: NULL (auto).
#' @param plateau_patience integer; patience (epochs) for ReduceLROnPlateau.
#'   Default: NULL (auto).
#' @param early_patience integer; patience (epochs) for EarlyStopping.
#'   Default: NULL (auto).
#' @param imbalance_thresh numeric in \eqn{[0,1]}; threshold for enabling imbalance handling.
#'   Default: NULL (auto).
#' @param auto_th_method character; "youden", "f1", or "auto" for thresholding.
#'   Default: NULL (auto).
#'
#' @param method_weights length-3 named numeric vector for composite weighting,
#'   names must be c("model","perm","shap"). Default: c(model=1, perm=1, shap=1).
#' @param ci_level numeric; CI level for composite bootstrap over runs (default 0.95).
#' @param strict_min_freq numeric; minimum cross-run hit frequency in per-run Top-K to keep a gene
#'   (e.g. 0.6 to 0.7). NULL disables.
#' @param strict_min_effect numeric; minimum composite mean to keep a gene on the composite scale.
#'   NULL disables.
#' @param strict_ci_gate logical; if TRUE, require composite CI lower bound > 0.
#'   Default: FALSE.
#' @param strict_knee logical; if TRUE, apply a knee cutoff on sorted composite means.
#'   Default: FALSE.
#' @param strict_nmax integer; optional hard cap after knee (e.g. 10). NULL disables.
#' @param apply_selection_to_plots logical; if TRUE, bar/UMAP use selected genes;
#'   default FALSE keeps original plotting behavior.
#'
#' @return A list with elements:
#' \item{params}{List of parsed parameters and metadata}
#' \item{importance_df}{Data frame with per-gene importance scores (3 methods + composite)}
#' \item{top_list}{List of top genes per run}
#' \item{final_top}{Vector of selected top genes after selection controls}
#' \item{composite_mat}{Matrix of per-run composite scores}
#' \item{glm_summary}{Summary of logistic regression on top genes (coefficients + formula)}
#'
#' @details
#' Output files in \code{out_dir}:
#' \itemize{
#'   \item FI_Boxplot_Top20.pdf: Top-20 feature importance distributions
#'   \item FI_Bar_TopComposite.pdf: Top-N composite feature importance bar plot (mean with 95% CI)
#'   \item FI_Density_AllMethods.pdf: Distribution comparison of all importance methods
#'   \item UMAP_TopSignatureGenes.pdf: UMAP projection of top genes
#'   \item Stability_TopGenes_Heatmap.pdf: Stability heatmap of top genes across runs
#' }
#'
#' @examples
#' \dontrun{
#' res <- fh_hunter(
#'   train_exp = X, 
#'   train_labels = Y,
#'   nshow = 40, score_index = 3, pick_index = 1,
#'   num_runs = 5, num_coregene = 50, n_likes = 30, n_interest = 30,# n_likes = "auto"
#'   strict_min_freq = 0.65, strict_ci_gate = TRUE, strict_knee = TRUE, strict_nmax = 10
#' )
#' head(res$importance_df)
#' }
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict coef cov median mad quantile reorder sd glm
#' @importFrom graphics plot
#' @importFrom FNN get.knnx
#' @importFrom rpart rpart
#' @importFrom partykit varimp
#' @importFrom logistf logistf
#' @importFrom magrittr %>%
#' @export
fh_hunter <- function(
    train_exp,
    train_labels,
    nshow,
    namesS = c("Accuracy","Recall","F-score"),
    score_index = 3,
    pick_index = 1,
    top_models_csv = NULL,
    num_runs = 10,
    num_coregene,
    n_likes,
    n_interest,
    seed = 424,
    out_dir = getwd(),
    # manual overrides for MLP (simple)
    learningratei = NULL,
    batch_sizei   = NULL,
    epochselecti  = NULL,
    dropoutratei  = NULL,
    cutoffi       = NULL,
    # permutation / SHAP
    perm_metric   = "f1",
    perm_nrep     = 3L,
    shap_nsim     = 10L,
    shap_subsample = 100L,
    # SVM / XGB / RF defaults
    svm_cost      = 1,
    svm_gamma     = NULL,
    xgb_nrounds   = 200,
    xgb_eta       = 0.05,
    xgb_max_depth = 4,
    rf_num_trees  = 500,
    rf_mtry       = NULL,
    # MLP trainer toggles (advanced; all NULL by default)
    shuffle            = TRUE,
    standardize        = NULL,
    hidden_units       = NULL,
    activation         = NULL,
    use_batchnorm      = NULL,
    l2                 = NULL,
    gaussian_noise_sd  = NULL,
    min_lr             = NULL,
    plateau_factor     = NULL,
    plateau_patience   = NULL,
    early_patience     = NULL,
    imbalance_thresh   = NULL,
    auto_th_method     = NULL,
    # >>> NEW: selection controls and weighting (defaults preserve old behavior) <<<
    method_weights = c(model = 1, perm = 1, shap = 1),
    ci_level = 0.95,
    strict_min_freq = NULL,
    strict_min_effect = NULL,
    strict_ci_gate = FALSE,
    strict_knee = FALSE,
    strict_nmax = NULL,
    apply_selection_to_plots = FALSE
){
  `%||%` <- function(a,b) if (is.null(a)) b else a
  .ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  ## >>> NEW: holders for MLP training curves (saved later) <<<
  p_train_loss <- NULL
  p_train_acc  <- NULL
  uses_shap_for_type <- function(model_type, method_weights) {
    parametric <- c("glmnet","lr","stepwise_lr","svm_linear","lda","qda","nb")
    if (model_type %in% parametric && (isTRUE(method_weights["shap"] == 0))) return(FALSE)
    if (model_type %in% parametric) return(FALSE)
    TRUE
  }
  .zcol <- function(x) { sdx <- stats::sd(x, na.rm = TRUE); if (is.na(sdx) || sdx == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / sdx }
  .zrob <- function(x) { m <- stats::median(x, na.rm = TRUE); s <- stats::mad(x, constant = 1.4826, na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1; (x - m) / s }
  .num1 <- function(pat, txt){
    m <- regmatches(txt, regexpr(pat, txt, perl=TRUE))
    if (length(m) == 0) return(NA_real_)
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
  .cutoff_from_modelstr <- function(txt){
    mm <- regexec("cutoff:auto\\(([0-9\\.]+),", txt)
    hit <- regmatches(txt, mm)[[1]]
    if (length(hit) >= 2) as.numeric(hit[2]) else .num1("cutoff\\s*:\\s*[0-9\\.]+", txt)
  }
  .pred_proba <- function(model, X, model_type){
    X <- as.matrix(X)
    
    .pos_ix <- function(pr, prefer = c("1","pos","positive","TRUE","Case")) {
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
                              error = function(e) NULL)
               if (is.null(pr)) stop("[DT] unsupported tree object for probability prediction.")
               if (is.null(dim(pr))) pr <- cbind(`0` = 1 - pr, `1` = pr)
               .clamp01(pr[, .pos_ix(pr)])
             }
           },
           
           tree = {  # alias of dt
             if (inherits(model, "rpart")) {
               pr <- predict(model, newdata = data.frame(X), type = "prob")
               if (is.null(dim(pr))) pr <- cbind(`0` = 1 - pr, `1` = pr)
               .clamp01(pr[, .pos_ix(pr)])
             } else {
               pr <- tryCatch(predict(model, newdata = data.frame(X), type = "prob"),
                              error = function(e) NULL)
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
               pos_lab <- intersect(model$levels, c("1","pos","positive","TRUE","Case"))
               if (!length(pos_lab)) pos_lab <- tail(model$levels, 1)
               .clamp01(pr[, pos_lab[1]])
             } else .clamp01(pr[, .pos_ix(pr)])
           },
           
           svm_rbf = {
             pr <- attr(predict(model, data.frame(X), probability = TRUE), "probabilities")
             if (is.null(pr)) stop("[SVM-rbf] Train with probability=TRUE")
             if (!is.null(model$levels) && all(model$levels %in% colnames(pr))) {
               pos_lab <- intersect(model$levels, c("1","pos","positive","TRUE","Case"))
               if (!length(pos_lab)) pos_lab <- tail(model$levels, 1)
               .clamp01(pr[, pos_lab[1]])
             } else .clamp01(pr[, .pos_ix(pr)])
           },
           
           ## -------------------- KNN (自定义“模型对象”) --------------------
           knn = {
             # 你的“模型”是 list(kind="knn", X_train, y, k, center, scale)
             if (!is.list(model) || is.null(model$kind) || model$kind != "knn")
               stop("[knn] Expect a list(kind='knn', ...).")
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
    if (length(models) == 1) return(.pred_proba(models[[1]], X, model_type))
    P <- vapply(models, function(m) .pred_proba(m, X, model_type), numeric(nrow(as.matrix(X))))
    rowMeans(cbind(P), na.rm = TRUE)
  }
  .perm_importance_fun <- function(
    predict_fun, X, y,
    metric = c("f1","prauc"),
    nrep = 3L,
    stratified = TRUE,
    seed = 424,
    progress = FALSE,
    pb_style = 3
  ) {
    metric <- match.arg(metric)
    set.seed(seed)
    
    base_p <- predict_fun(X)
    base_s <- switch(
      metric,
      f1 = {
        pr <- as.integer(base_p >= 0.5)
        tp <- sum(pr == 1 & y == 1); fp <- sum(pr == 1 & y == 0); fn <- sum(pr == 0 & y == 1)
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
        s <- s + switch(
          metric,
          f1 = {
            pr <- as.integer(prp >= 0.5)
            tp <- sum(pr == 1 & y == 1); fp <- sum(pr == 1 & y == 0); fn <- sum(pr == 0 & y == 1)
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
    X <- as.matrix(X); p <- ncol(X)
    cn <- colnames(X)
    
    safe_vec <- function(v, names_ref) {
      v[!is.finite(v)] <- 0
      names(v) <- names_ref
      as.numeric(v)
    }
    
    out <- switch(
      model_type,
      
      # 1) MLP: use the first linear layer "slayer" absolute weights
      mlp = {
        w <- tryCatch(keras3::get_weights(keras3::get_layer(model, "slayer"))[[1]], error = function(e) NULL)
        if (is.null(w)) rep(0, p) else safe_vec(abs(as.numeric(w)), cn)
      },
      
      # 2) Random Forest (ranger): variable.importance
      rf = {
        vi <- model$variable.importance
        if (is.null(vi)) {
          out0 <- rep(0, p); names(out0) <- cn; as.numeric(out0)
        } else {
          v <- rep(0, p); names(v) <- cn
          v[names(vi)] <- as.numeric(vi)
          safe_vec(v[cn], cn)
        }
      },
      
      # 3) GLMNET (cv.glmnet or glmnet): |beta| * sd(X)
      glmnet = {
        bvec <- tryCatch(as.numeric(stats::coef(model, s = "lambda.min")), error = function(e) as.numeric(stats::coef(model)))
        if (length(bvec) == p + 1) {
          b <- bvec[-1]
          sx <- apply(X, 2, stats::sd); sx[!is.finite(sx) | sx == 0] <- 1
          safe_vec(abs(b) * sx, cn)
        } else {
          rep(0, p)
        }
      },
      
      # 4) XGBoost: gain-based importance mapped to columns
      xgb = {
        im <- tryCatch(xgboost::xgb.importance(feature_names = cn, model = model), error = function(e) NULL)
        if (is.null(im) || nrow(im) == 0) {
          rep(0, p)
        } else {
          v <- setNames(im$Gain, im$Feature)
          vv <- v[cn]; vv[is.na(vv)] <- 0
          safe_vec(as.numeric(vv), cn)
        }
      },
      
      # 5) GBM (gbm package): relative influence
      gbm = {
        ri <- tryCatch({
          nt <- if (!is.null(model$n.trees)) model$n.trees else NULL
          if (is.null(nt)) gbm::relative.influence(model, normalized = TRUE) else gbm::relative.influence(model, n.trees = nt, normalized = TRUE)
        }, error = function(e) NULL)
        if (is.null(ri)) {
          rep(0, p)
        } else {
          v <- rep(0, p); names(v) <- cn
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
            scl <- model$x.scale$`scaled:scale`; scl[is.na(scl) | scl == 0] <- 1
            Xs <- sweep(sweep(X, 2, cen, "-"), 2, scl, "/")
          }
          n <- nrow(Xs)
          idx <- if (n > rbf_grad_subsample) sample(n, rbf_grad_subsample) else seq_len(n)
          Xsub <- Xs[idx, , drop = FALSE]
          SV  <- model$SV
          al  <- as.numeric(model$coefs)
          imp_sum <- numeric(p)
          for (k in seq_len(nrow(Xsub))) {
            xk <- Xsub[k, ]
            diff <- sweep(SV, 2, xk, "-")
            d2   <- rowSums(diff * diff)
            Kvec <- exp(-gamma * d2)
            wsum <- colSums( (al * Kvec) * (SV - matrix(xk, nrow(SV), p, byrow = TRUE)) )
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
          v <- rep(0, p); names(v) <- cn
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
            vv <- rep(0, p); names(vv) <- cn
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
          v <- rep(0, p); names(v) <- cn
          for (nm in intersect(names(tabs), cn)) {
            tb <- tabs[[nm]]
            # numeric features have a 2x2 (or 2x?) matrix with mean/sd per class
            if (is.matrix(tb) && all(c("mean","sd") %in% tolower(colnames(tb)))) {
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
          b <- co[names(co) %in% paste0(cn)]  # match by exact names if possible
          # fallback: drop intercept and assume order
          if (length(b) == 0L && length(co) >= 2L) b <- co[-1]
          b <- as.numeric(b)
          if (length(b) != p) {
            # align by colnames when possible
            tmp <- rep(0, p); names(tmp) <- cn
            nb <- names(co); nb <- nb[nb != "(Intercept)"]
            tmp[nb] <- as.numeric(co[nb])
            b <- as.numeric(tmp)
          }
          sx <- apply(X, 2, stats::sd); sx[!is.finite(sx) | sx == 0] <- 1
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
            tmp <- rep(0, p); names(tmp) <- cn
            nb <- names(co); nb <- nb[nb != "(Intercept)"]
            tmp[nb] <- as.numeric(co[nb])
            b <- as.numeric(tmp)
          }
          sx <- apply(X, 2, stats::sd); sx[!is.finite(sx) | sx == 0] <- 1
          safe_vec(abs(b) * sx, cn)
        }
      },
      
      # 13) Decision tree (rpart/partykit): variable importance
      dt = {
        # rpart
        vi <- tryCatch(model$variable.importance, error = function(e) NULL)
        if (!is.null(vi)) {
          v <- rep(0, p); names(v) <- cn
          v[names(vi)] <- as.numeric(vi)
          safe_vec(v[cn], cn)
        } else {
          # partykit varimp if available
          v2 <- tryCatch({
            if (requireNamespace("partykit", quietly = TRUE)) {
              as.numeric(partykit::varimp(model))
            } else NULL
          }, error = function(e) NULL)
          if (is.null(v2)) rep(0, p) else {
            nv <- partykit::varimp(model)
            v <- rep(0, p); names(v) <- cn
            v[names(nv)] <- as.numeric(nv)
            safe_vec(v[cn], cn)
          }
        }
      },
      
      tree = {  # alias
        vi <- tryCatch(model$variable.importance, error = function(e) NULL)
        if (!is.null(vi)) {
          v <- rep(0, p); names(v) <- cn
          v[names(vi)] <- as.numeric(vi)
          safe_vec(v[cn], cn)
        } else rep(0, p)
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
    update_every = 5
  ) {
    X <- as.matrix(X)
    n <- nrow(X); p <- ncol(X)
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
  
  # decide whether to compute SHAP for a given model_type
  .use_shap <- function(model_type) {
    # parametric, intrinsically interpretable models: skip SHAP
    skip_set <- c("glmnet","lr","stepwise_lr","svm_linear","lda","qda","nb")
    !(model_type %in% skip_set)
  }
  
  # >>> NEW: knee finder helper for selection <<<
  find_knee <- function(x) {
    n <- length(x); if (n < 3) return(n)
    i <- 2:(n - 1)
    d1 <- x[i - 1] - x[i]
    d2 <- d1[-1] - d1[-length(d1)]
    j  <- which.max(pmax(d2, 0))
    if (!length(j)) return(n)
    i[j + 1]
  }
  
  stopifnot(!missing(train_exp), !missing(train_labels))
  stopifnot(is.numeric(num_runs) && num_runs >= 1)
  stopifnot(!missing(num_coregene), !missing(n_likes), !missing(n_interest))
  
  perm_metric <- perm_metric[1]
  
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  plots_dir <- file.path(out_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  set.seed(seed)
  train_exp <- as.matrix(train_exp)
  n <- nrow(train_exp); p <- ncol(train_exp)
  train_y <- as.numeric(train_labels)
  if (any(!train_y %in% c(0,1))) {
    if (is.factor(train_labels) || is.character(train_labels)) {
      ty <- tolower(as.character(train_labels))
      train_y <- ifelse(ty %in% c("1","pos","positive","yes","true","case"), 1,
                        ifelse(ty %in% c("0","neg","negative","no","false","control"), 0, NA))
      if (any(is.na(train_y))) stop("train_labels cannot be converted to 0/1.")
    } else stop("train_labels must be 0/1 (or convertible).")
  }
  
  if (is.null(top_models_csv)) {
    top_models_csv <- file.path("heatmap", paste0(nshow, "_top_", namesS[score_index], "_models.csv"))
  }
  if (!file.exists(top_models_csv)) stop("Leaderboard CSV not found: ", top_models_csv)
  Top_models <- utils::read.csv(top_models_csv, stringsAsFactors = FALSE)
  if (!("Model" %in% names(Top_models))) stop("Leaderboard CSV is missing column 'Model'.")
  
  pick_row <- pick_index
  model_str <- Top_models$Model[pick_row]
  if (is.na(model_str)) stop("Invalid pick index: ", pick_row)
  
  detect_type <- function(s) {
    s_low <- tolower(s)
    first_pos <- function(pat) {
      m <- regexpr(pat, s_low, perl = TRUE, ignore.case = TRUE)
      if (m[1] == -1) Inf else as.integer(m[1])
    }
    primary <- list(
      svm   = "\\bsvm\\b",
      rf    = "(\\brf\\b|random\\s*forest|\\branger\\b)",
      xgb   = "(\\bxgboost\\b|\\bxgb\\b)",
      gbm   = "\\bgbm\\b",
      mlp   = "(\\bnn-mlp\\b|\\bmlp\\b|keras|tensorflow)",
      lda   = "\\blda\\b",
      qda   = "\\bqda\\b",
      nb    = "(naive\\s*bayes|\\bnaivebayes\\b|\\bnb\\b)",
      knn   = "\\bknn\\b",
      dt    = "(\\bdt\\b|decision\\s*tree|\\btree\\b)",
      lr    = "(^|[^a-z])lr([^a-z]|$)|\\blogistic\\b",
      stepwise_lr = "(stepwise.*lr|step\\s*wise.*lr|aic\\+lr)"
    )
    secondary_glmnet <- "(\\benr\\b|elastic\\s*net|\\brr\\b|\\bridge\\b|lasso\\.?r?\\b|\\bglmnet\\b)"
    pos_primary <- vapply(names(primary), function(k) first_pos(primary[[k]]), numeric(1))
    if (!all(is.infinite(pos_primary))) {
      typ <- names(pos_primary)[which.min(pos_primary)]
      if (typ == "svm") {
        if (grepl("kernel\\s*:\\s*linear|linear\\s*kernel", s_low, perl = TRUE)) return("svm_linear") else return("svm_rbf")
      }
      return(typ)
    }
    if (regexpr(secondary_glmnet, s_low, perl = TRUE)[1] != -1) return("glmnet")
    "unknown"
  }
  model_type <- detect_type(model_str)
  if (model_type == "unknown") stop("Unknown model type: ", model_str)
  
  cutoff_p <- .cutoff_from_modelstr(model_str)
  parsed <- list(model_type = model_type, cutoff = cutoff_p)
  
  if (model_type == "mlp") {
    lr_p <- .num1("lr\\s*:\\s*[0-9\\.]+", model_str)
    bs_p <- .num1("bs\\s*:\\s*[0-9]+",    model_str)
    ep_p <- .num1("ep\\s*:\\s*[0-9]+",    model_str)
    dr_p <- .num1("dropout\\s*:\\s*[0-9\\.]+", model_str)
    learningratei <- learningratei %||% lr_p
    batch_sizei   <- batch_sizei   %||% bs_p
    epochselecti  <- epochselecti  %||% ep_p
    dropoutratei  <- dropoutratei  %||% dr_p
    cutoffi       <- cutoffi       %||% cutoff_p
    if (any(is.na(c(learningratei, batch_sizei, epochselecti, dropoutratei, cutoffi)))) {
      stop("[MLP] failed to parse params; please provide learningratei/batch_sizei/epochselecti/dropoutratei/cutoffi")
    }
    parsed <- c(parsed, list(lr=learningratei, bs=batch_sizei, ep=epochselecti, dropout=dropoutratei))
  } else if (model_type == "rf") {
    rf_mtry <- rf_mtry %||% .num1("mtry\\s*=\\s*[0-9]+", model_str) %||% floor(sqrt(p))
    rf_num_trees <- rf_num_trees %||% .num1("(ntree|num\\.trees)\\s*=\\s*[0-9]+", model_str) %||% 500
    parsed <- c(parsed, list(mtry=rf_mtry, num.trees=rf_num_trees))
  } else if (model_type == "glmnet") {
    s_low <- tolower(model_str)
    alpha_p <- .num1("alpha\\s*[:=]\\s*[0-9\\.]+", model_str)
    if (is.na(alpha_p)) {
      if (grepl("lasso", s_low)) alpha_p <- 1 else if (grepl("ridge|\\brr\\b", s_low)) alpha_p <- 0 else if (grepl("elastic\\s*net|\\benr\\b", s_low)) alpha_p <- 0.5
    }
    parsed <- c(parsed, list(alpha = alpha_p))
  } else if (model_type == "xgb") {
    md_p <- .num1("max_depth\\s*=\\s*[0-9]+", model_str); if (is.finite(md_p)) xgb_max_depth <- md_p
    et_p <- .num1("(eta|lr)\\s*[:=]\\s*[0-9\\.]+", model_str); if (is.finite(et_p)) xgb_eta <- et_p
    nr_p <- .num1("(nrounds|trees)\\s*=\\s*[0-9]+", model_str); if (is.finite(nr_p)) xgb_nrounds <- nr_p
    parsed <- c(parsed, list(max_depth=xgb_max_depth, eta=xgb_eta, nrounds=xgb_nrounds))
  } else if (startsWith(model_type,"svm")) {
    c_p <- .num1("C\\s*=\\s*[0-9\\.]+", model_str); if (is.finite(c_p)) svm_cost <- c_p
    g_p <- .num1("gamma\\s*=\\s*[0-9\\.]+", model_str); if (is.finite(g_p)) svm_gamma <- g_p
    if (is.null(svm_gamma)) svm_gamma <- 1 / ncol(train_exp)
    parsed <- c(parsed, list(cost=svm_cost, gamma=svm_gamma))
  }
  
  message("[hunter] pick=", pick_index, " | type=", model_type, " | parsed=", paste(names(parsed), parsed, collapse=", "))
  
  dimn <- list(colnames(train_exp), NULL)
  A_model <- matrix(NA_real_, nrow = p, ncol = num_runs, dimnames = dimn)
  B_perm  <- matrix(NA_real_, nrow = p, ncol = num_runs, dimnames = dimn)
  C_shap  <- matrix(0, nrow = p, ncol = num_runs, dimnames = dimn)
  fits    <- vector("list", length = num_runs)
  top_list <- vector("list", length = num_runs)
  message("Step 1/3: Computing model-based importance (multi-run with seeds)")
  pb <- txtProgressBar(min = 0, max = num_runs, style = 3)
  
  for (run in seq_len(num_runs)) {
    set.seed(seed + run)
    
    fit <- switch(
      model_type,
      
      mlp = {
        get_or <- function(sym_name, default) {
          val <- get0(sym_name, inherits = TRUE, ifnotfound = NULL)
          if (is.null(val)) return(default)
          if (is.logical(val) && length(val) == 1 && is.na(val)) return(default)
          if (is.numeric(val) && length(val) == 1 && is.na(val)) return(default)
          val
        }
        
        seed_i <- as.integer(get_or("seed", 0)) + run
        set.seed(seed_i)
        if (reticulate::py_module_available("tensorflow")) {
          tensorflow::set_random_seed(seed_i)
          tensorflow::tf$random$set_seed(seed_i)
        }
        
        bn_on       <- isTRUE(get_or("use_batchnorm", TRUE))
        l2_reg      <- get_or("l2", 1e-4)
        gnoise_sd   <- get_or("gaussian_noise_sd", 0.0)
        act_hidden  <- get_or("activation", "relu")
        min_lr_i    <- get_or("min_lr", 1e-5)
        plat_fac    <- get_or("plateau_factor", 0.5)
        plat_pat    <- as.integer(get_or("plateau_patience", 6))
        es_pat      <- as.integer(get_or("early_patience", 9))
        shuffle_i   <- isTRUE(get_or("shuffle", TRUE))
        
        units_vec <- get_or("hidden_units", c(32L, 16L, 8L))
        units_vec <- as.integer(units_vec[is.finite(units_vec) & units_vec > 0])
        
        y01 <- as.integer(train_y)
        pos_rate <- mean(y01 == 1); neg_rate <- 1 - pos_rate
        use_imbalance <- (min(pos_rate, neg_rate) < (get_or("imbalance_thresh", 0.35)))
        class_weighti <- if (use_imbalance) list(`0` = 0.5 / max(neg_rate, 1e-12), `1` = 0.5 / max(pos_rate, 1e-12)) else NULL
        out_bias <- if (use_imbalance) keras3::initializer_constant(log(pos_rate / max(neg_rate, 1e-12))) else "zeros"
        
        use_l2 <- is.numeric(l2_reg) && length(l2_reg) == 1 && is.finite(l2_reg) && l2_reg > 0
        reg <- if (use_l2) keras3::regularizer_l2(l2_reg) else NULL
        
        inputs <- keras3::layer_input(shape = c(p))
        x <- inputs
        if (is.numeric(gnoise_sd) && length(gnoise_sd) == 1 && is.finite(gnoise_sd) && gnoise_sd > 0) {
          x <- x |> keras3::layer_gaussian_noise(stddev = gnoise_sd)
        }
        x <- x |> keras3::layer_dense(units = 1, activation = "linear", name = "slayer", kernel_regularizer = reg)
        
        if (length(units_vec)) {
          for (u in units_vec) {
            x <- x |> keras3::layer_dense(units = as.integer(u), activation = act_hidden, kernel_regularizer = reg)
            if (bn_on) x <- x |> keras3::layer_batch_normalization()
            x <- x |> keras3::layer_dropout(rate = dropoutratei)
          }
        }
        
        outputs <- x |> keras3::layer_dense(units = 1, activation = "sigmoid", bias_initializer = out_bias)
        model <- keras3::keras_model(inputs = inputs, outputs = outputs)
        
        opt <- keras3::optimizer_adam(learning_rate = learningratei, clipnorm = 1.0)
        keras3::compile(
          model,
          loss = "binary_crossentropy",
          optimizer = opt,
          metrics = c(
            keras3::metric_binary_accuracy(name = "accuracy"),
            keras3::metric_auc(name = "auc_roc"),
            keras3::metric_auc(curve = "PR", name = "prc")
          )
        )
        
        cbs <- list(
          keras3::callback_reduce_lr_on_plateau(monitor = "val_prc", mode = "max", factor = plat_fac, patience = plat_pat, min_lr = min_lr_i, verbose = 0),
          keras3::callback_early_stopping(monitor = "val_prc", mode = "max", patience = es_pat, restore_best_weights = TRUE)
        )
        
        history <- keras3::fit(model,
                               x = train_exp,
                               y = train_y,
                               epochs = epochselecti,
                               batch_size = batch_sizei,
                               shuffle = shuffle_i,
                               validation_split = 0.2,
                               class_weight = class_weighti,
                               callbacks = cbs,
                               verbose = 0)
        
        ## >>> CHANGED: cache ggplot objects for curves (only once at run==1); do NOT save yet <<<
        if (run == 1) {
          loss_vec <- as.numeric(history$metrics$loss)
          acc_vec  <- if (!is.null(history$metrics$accuracy)) {
            as.numeric(history$metrics$accuracy)
          } else if (!is.null(history$metrics$binary_accuracy)) {
            as.numeric(history$metrics$binary_accuracy)
          } else {
            NULL
          }
          n_ep <- length(loss_vec)
          if (!is.null(acc_vec)) n_ep <- min(n_ep, length(acc_vec))
          if (n_ep > 0) {
            train_metrics <- data.frame(
              Epoch    = seq_len(n_ep),
              Loss     = loss_vec[seq_len(n_ep)],
              Accuracy = if (is.null(acc_vec)) NA_real_ else acc_vec[seq_len(n_ep)]
            )
            
            ## build but don't save
            p_train_loss <- ggplot2::ggplot(train_metrics, ggplot2::aes(x = Epoch, y = Loss)) +
              ggplot2::geom_line(color = "firebrick") +
              ggplot2::labs(title = "Training Loss Curve", y = "Loss") +
              ggplot2::theme_minimal()
            
            if (!all(is.na(train_metrics$Accuracy))) {
              p_train_acc <- ggplot2::ggplot(train_metrics, ggplot2::aes(x = Epoch, y = Accuracy)) +
                ggplot2::geom_line(color = "dodgerblue") +
                ggplot2::labs(title = "Training Accuracy Curve", y = "Accuracy") +
                ggplot2::theme_minimal()
            }
          }
        }
        
        model
      },
      
      rf = {
        ranger::ranger(
          dependent.variable.name = "y",
          data  = data.frame(y = factor(train_y), train_exp),
          probability = TRUE,
          importance  = "impurity",
          num.trees   = as.integer(parsed$num.trees),
          mtry        = as.integer(parsed$mtry),
          seed        = seed + run
        )
      },
      
      glmnet = {
        glmnet::cv.glmnet(
          x = train_exp,
          y = train_y,
          family = "binomial",
          alpha  = parsed$alpha %||% 1,
          type.measure = "class",
          standardize = TRUE
        )
      },
      
      xgb = {
        dtrain <- xgboost::xgb.DMatrix(data = train_exp, label = train_y)
        xgboost::xgboost(
          data = dtrain,
          objective = "binary:logistic",
          nrounds = as.integer(xgb_nrounds),
          max_depth = as.integer(xgb_max_depth),
          eta = xgb_eta,
          subsample = 0.8, colsample_bytree = 0.8,
          verbose = 0
        )
      },
      
      svm_linear = {
        e1071::svm(x = train_exp, y = factor(train_y), kernel = "linear", probability = TRUE, cost = svm_cost, scale = TRUE)
      },
      
      svm_rbf = {
        gamma_final <- if (is.null(svm_gamma)) 1 / ncol(train_exp) else as.numeric(svm_gamma)
        parsed$gamma <- gamma_final
        e1071::svm(x = train_exp, y = factor(train_y), kernel = "radial", probability = TRUE, cost = svm_cost, gamma = gamma_final, scale = TRUE)
      },
      
      lda = MASS::lda(x = train_exp, grouping = factor(train_y)),
      qda = MASS::qda(x = train_exp, grouping = factor(train_y)),
      nb  = e1071::naiveBayes(x = train_exp, y = factor(train_y)),
      gbm = {
        gbm::gbm(
          formula = y ~ .,
          data = data.frame(y = train_y, train_exp),
          distribution = "bernoulli",
          n.trees = 200L,              # 可暴露为参数
          interaction.depth = 3L,      # 可暴露为参数
          shrinkage = 0.05,            # 可暴露为参数
          n.minobsinnode = 10L,
          bag.fraction = 0.8,
          train.fraction = 1.0,
          keep.data = FALSE,
          verbose = FALSE
        )
      },
      
      lr = {
        stats::glm(y ~ ., data = data.frame(y = train_y, train_exp),
                   family = stats::binomial(link = "logit"))
      },
      
      stepwise_lr = {
        base_fit <- stats::glm(y ~ ., data = data.frame(y = train_y, train_exp),
                               family = stats::binomial(link = "logit"))
        if (requireNamespace("MASS", quietly = TRUE)) {
          suppressWarnings(MASS::stepAIC(base_fit, trace = FALSE))
        } else {
          base_fit
        }
      },
      
      dt = {
        if (requireNamespace("rpart", quietly = TRUE)) {
          rpart::rpart(y ~ ., data = data.frame(y = factor(train_y), train_exp),
                       method = "class", parms = list(split = "gini"))
        } else {
          stop("rpart not installed but model_type=dt/tree was selected.")
        }
      },
      
      tree = {
        if (requireNamespace("rpart", quietly = TRUE)) {
          rpart::rpart(y ~ ., data = data.frame(y = factor(train_y), train_exp),
                       method = "class", parms = list(split = "gini"))
        } else {
          stop("rpart not installed but model_type=dt/tree was selected.")
        }
      },
      
      knn = {
        # “训练”就是保存标准化后的训练集与 y，以及 k
        k_knn <- max(3L, min(15L, as.integer(round(sqrt(nrow(train_exp))))))
        cen <- colMeans(train_exp)
        scl <- apply(train_exp, 2, stats::sd); scl[!is.finite(scl) | scl == 0] <- 1
        Xs  <- scale(train_exp, center = cen, scale = scl)
        list(kind = "knn", X_train = Xs, y = factor(train_y), k = k_knn,
             center = cen, scale = scl)
      }
    )
    
    fits[[run]] <- fit
    W_model <- .inner_importance(fit, train_exp, model_type)
    A_model[, run] <- W_model
    
    composite_run_tmp <- .zrob(A_model[, run])
    top_idx <- order(composite_run_tmp, decreasing = TRUE)[seq_len(min(num_coregene, length(composite_run_tmp)))]
    top_list[[run]] <- colnames(train_exp)[top_idx]
    
    setTxtProgressBar(pb, run)
  }
  close(pb)
  
  
  message("Step 2/3: Computing permutation importance (feature shuffling)")
  # Ensemble-first evaluation for permutation and SHAP
  predict_ens <- function(X) .pred_proba_ensemble(fits, X, model_type)
  
  .pb <- function(n) utils::txtProgressBar(min = 0, max = n, style = 3)
  .bump <- function(pb, i) utils::setTxtProgressBar(pb, i)
  
  steps <- c("perm_ensemble", "shap_ensemble", "broadcast", "composite", "aggregate", "bootstrap")
  pb_blk <- .pb(length(steps)); st <- 0
  
  
  ## 1) perm per-run
  st <- st + 1; .bump(pb_blk, st)
  pb_perm <- utils::txtProgressBar(min = 0, max = num_runs, style = 3)
  for (r in seq_len(num_runs)) {
    # 针对第 r 个模型的预测函数
    predict_fun_r <- function(X) .pred_proba(fits[[r]], X, model_type)
    # 逐 run 计算 permutation importance（不同 seed，避免完全一致）
    B_perm[, r] <- .perm_importance_fun(
      predict_fun = predict_fun_r,
      X = train_exp, y = train_y,
      metric = perm_metric, nrep = perm_nrep,
      stratified = TRUE, seed = seed + 1000 + r,
      progress = FALSE
    )
    utils::setTxtProgressBar(pb_perm, r)
  }
  close(pb_perm)
  
  message("Step 3/3: Computing SHAP importance (skipped for parametric models) (policy=", ifelse(.use_shap(model_type) || identical(model_type,"xgb"), "enabled", "skipped"), ")")
  ## 2) shap per-run
  st <- st + 1; .bump(pb_blk, st)
  pb_shap <- utils::txtProgressBar(min = 0, max = num_runs, style = 3)
  
  for (r in seq_len(num_runs)) {
    if (identical(model_type, "xgb")) {
      # XGBoost: use native TreeSHAP contributions
      contrib <- stats::predict(fits[[r]], newdata = train_exp, predcontrib = TRUE)
      if (ncol(contrib) == ncol(train_exp) + 1) {
        contrib <- contrib[, -ncol(contrib), drop = FALSE]  # drop bias term
      }
      C_shap[, r] <- as.numeric(colMeans(abs(contrib), na.rm = TRUE))
      names(C_shap[, r]) <- colnames(train_exp)
      
    } else if (.use_shap(model_type)) {
      # nonparametric/black-box models: approximate SHAP via fastshap
      predict_fun_r <- function(X) .pred_proba(fits[[r]], X, model_type)
      C_shap[, r] <- .shap_importance_fun(
        predict_fun = predict_fun_r,
        X = train_exp,
        shap_nsim = shap_nsim,
        shap_subsample = shap_subsample,
        seed = seed + 2000 + r,
        progress = FALSE
      )
      
    } else {
      # parametric models: skip SHAP to save time (fill zeros to keep shapes)
      C_shap[, r] <- rep(0, ncol(train_exp))
      names(C_shap[, r]) <- colnames(train_exp)
    }
    
    utils::setTxtProgressBar(pb_shap, r)
  }
  close(pb_shap)
  close(pb_blk)
  
  ## --- Robust z-score per run (by column) ---
  zA_mat <- apply(A_model, 2, .zrob)  # p x num_runs
  zB_mat <- apply(B_perm,  2, .zrob)  # p x num_runs
  zC_mat <- apply(C_shap,  2, .zrob)  # p x num_runs
  
  if (is.vector(zA_mat)) zA_mat <- matrix(zA_mat, nrow = p)
  if (is.vector(zB_mat)) zB_mat <- matrix(zB_mat, nrow = p)
  if (is.vector(zC_mat)) zC_mat <- matrix(zC_mat, nrow = p)
  
  dimnames(zA_mat) <- dimnames(A_model)
  dimnames(zB_mat) <- dimnames(B_perm)
  dimnames(zC_mat) <- dimnames(C_shap)
  
  ## --- Helper: renormalize weights over enabled components ---
  .renorm_weights <- function(w, enabled = c("model","perm","shap")) {
    w2 <- w
    disable <- setdiff(names(w), enabled)
    if (length(disable)) w2[disable] <- 0
    s <- sum(w2)
    if (s > 0) w2 <- w2 / s
    w2
  }
  
  ## --- Check and normalize weights ---
  w <- method_weights
  if (is.null(names(w)) || !all(sort(names(w)) == c("model","perm","shap"))) {
    stop("method_weights must be a named numeric vector with names c('model','perm','shap').")
  }
  
  ## --- Detect whether SHAP is effectively available this run ---
  ## Policy-free detection: if C_shap has any finite non-zero entry, treat as effective.
  shap_has_signal <- isTRUE(any(is.finite(C_shap) & (abs(C_shap) > 0)))
  shap_effective  <- shap_has_signal
  
  ## --- Finalize weights: drop SHAP if ineffective, then renormalize ---
  if (!shap_effective) {
    w <- .renorm_weights(w, enabled = c("model","perm"))
  } else {
    w <- .renorm_weights(w, enabled = c("model","perm","shap"))
  }
  
  ## --- Fuse per run (align columns r for zA/zB/zC) ---
  composite_mat <- matrix(NA_real_, nrow = p, ncol = num_runs,
                          dimnames = list(colnames(train_exp), NULL))
  
  for (r in seq_len(num_runs)) {
    zA <- zA_mat[, r]
    zB <- zB_mat[, r]
    zC <- if (shap_effective) zC_mat[, r] else rep(0, p)  # avoid 0 * NA
    
    ## Gentle clipping to avoid single-run outliers
    cap <- 5
    zA <- pmax(pmin(zA, cap), -cap)
    zB <- pmax(pmin(zB, cap), -cap)
    zC <- pmax(pmin(zC, cap), -cap)
    
    composite_mat[, r] <- w["model"] * zA + w["perm"] * zB + w["shap"] * zC
  }
  
  ## --- Per-method summaries ---
  savg <- rowMeans(A_model, na.rm = TRUE)
  ssd  <- apply(A_model, 1, stats::sd, na.rm = TRUE)
  perm_avg <- rowMeans(B_perm, na.rm = TRUE)
  perm_sd  <- apply(B_perm,  1, stats::sd, na.rm = TRUE)
  
  ## SHAP summary: if ineffective, report zeros to keep the columns while signaling "not used".
  if (shap_effective) {
    shap_avg <- rowMeans(C_shap, na.rm = TRUE)
    shap_sd  <- apply(C_shap,  1, stats::sd, na.rm = TRUE)
  } else {
    shap_avg <- rep(0, p)
    shap_sd  <- rep(0, p)
  }
  
  ## --- Composite mean (final ranking basis) ---
  composite_avg <- rowMeans(composite_mat, na.rm = TRUE)
  
  ## --- Bootstrap CI over runs for the composite ---
  alpha <- 1 - ci_level
  comp_ci <- t(vapply(
    seq_len(p),
    function(j) .boot_ci(composite_mat[j, ], R = 1000, alpha = alpha),
    numeric(3L)
  ))
  colnames(comp_ci) <- c("Composite_mean", "Composite_lo", "Composite_hi")
  
  ## --- Assemble importance table ---
  importance_df <- data.frame(
    Gene = colnames(train_exp),
    Model_mean = savg, Model_sd = ssd,
    Perm_mean  = perm_avg, Perm_sd  = perm_sd,
    SHAP_mean  = shap_avg, SHAP_sd  = shap_sd,
    Composite  = composite_avg,
    Composite_lo = comp_ci[, "Composite_lo"],
    Composite_hi = comp_ci[, "Composite_hi"],
    stringsAsFactors = FALSE
  )
  
  f_imp <- file.path(out_dir, "Importance_Table.csv")
  utils::write.csv(importance_df, f_imp, row.names = FALSE)
  message(sprintf("[%s] Saved importance table to: %s", .ts(), f_imp))
  
  # --- Selection filters: apply strict_* to ALL plots ---
  selected_genes <- NULL
  rank_df <- importance_df[order(importance_df$Composite, decreasing = TRUE), ]
  rank_df$FreqProp <- 0
  # build a temporary frequency map from top_list (per-run Top-K)
  tmp_tab <- as.data.frame(table(unlist(lapply(top_list, function(genelist) head(genelist, n_interest)))))
  if (nrow(tmp_tab)) {
    colnames(tmp_tab) <- c("Gene","Frequency")
    freq_map <- setNames(tmp_tab$Frequency / num_runs, as.character(tmp_tab$Gene))
    rank_df$FreqProp <- ifelse(rank_df$Gene %in% names(freq_map), freq_map[rank_df$Gene], 0)
  }
  keep <- rep(TRUE, nrow(rank_df))
  if (is.numeric(strict_min_freq) && is.finite(strict_min_freq)) keep <- keep & (rank_df$FreqProp >= strict_min_freq)
  if (is.numeric(strict_min_effect) && is.finite(strict_min_effect)) keep <- keep & (rank_df$Composite >= strict_min_effect)
  if (isTRUE(strict_ci_gate)) keep <- keep & (rank_df$Composite_lo > 0)
  cand_df <- rank_df[keep, , drop = FALSE]
  if (isTRUE(strict_knee) && nrow(cand_df) >= 3) {
    k <- find_knee(cand_df$Composite)
    cand_df <- cand_df[seq_len(k), , drop = FALSE]
  }
  if (is.numeric(strict_nmax) && is.finite(strict_nmax) && strict_nmax > 0) {
    cand_df <- cand_df[seq_len(min(nrow(cand_df), as.integer(strict_nmax))), , drop = FALSE]
  }
  if (!nrow(cand_df)) {
    cand_df <- rank_df[seq_len(min(nrow(rank_df), n_likes)), , drop = FALSE]
  }
  selected_genes <- cand_df$Gene
  plot_set <- if (isTRUE(apply_selection_to_plots)) selected_genes else colnames(train_exp)
  
  # Determine effective n_likes for downstream (UMAP/logistic)
  # When n_likes == "auto", count genes within plot_set whose 95% CI does not cross 0.
  if (is.character(n_likes) && tolower(n_likes) == "auto") {
    n_likes_eff <- sum(importance_df$Gene %in% plot_set & importance_df$Composite_lo > 0, na.rm = TRUE)
    n_likes_eff <- max(2L, as.integer(n_likes_eff))  # keep at least 2 for UMAP
    message(sprintf(
      "[fh_hunter] n_likes='auto' -> using %d genes with %.0f%% CI > 0 within plot_set",
      n_likes_eff, ci_level * 100
    ))
  } else {
    n_likes_eff <- as.integer(n_likes)
  }
  
  # Top lists summary (respect plot_set)
  top_df <- as.data.frame(table(unlist(lapply(top_list, function(genelist) {
    head(intersect(genelist, plot_set), n_interest)
  }))))
  colnames(top_df) <- c("Gene", "Frequency")
  if (nrow(top_df)) {
    top_df$Gene <- factor(top_df$Gene, levels = top_df$Gene[order(top_df$Frequency, decreasing = TRUE)])
  }
  
  f_top <- file.path(out_dir, "FI_TopGeneFrequencies.csv")
  
  
  ## >>> NEW: save MLP training curves at the end (only if MLP and plots exist) <<<
  if (identical(model_type, "mlp") && (!is.null(p_train_loss) || !is.null(p_train_acc))) {
    f_loss <- file.path(plots_dir, "Training_Loss.pdf")
    f_acc  <- file.path(plots_dir, "Training_Accuracy.pdf")
    
    if (!is.null(p_train_loss)) {
      ggplot2::ggsave(f_loss, p_train_loss, width = 6, height = 4)
      message(sprintf("[%s] Saved plot to: %s", .ts(), f_loss))
    }
    if (!is.null(p_train_acc)) {
      ggplot2::ggsave(f_acc, p_train_acc, width = 6, height = 4)
      message(sprintf("[%s] Saved plot to: %s", .ts(), f_acc))
    }
  }
  ## Decide whether SHAP was computed for this model_type
  include_shap <- isTRUE(uses_shap_for_type(model_type, method_weights))
  
  ## Boxplot (Top20) across runs
  topN <- min(20, nrow(importance_df))
  topN_genes <- importance_df |>
    dplyr::arrange(dplyr::desc(Composite)) |>
    dplyr::pull(Gene)
  topN_genes <- intersect(topN_genes, plot_set)
  topN_genes <- head(topN_genes, topN)
  
  if (length(topN_genes) > 0) {
    df_list_box <- list(
      tibble::tibble(
        Gene = rep(topN_genes, each = num_runs),
        Method = "Model",
        Value = as.vector(A_model[topN_genes, , drop = FALSE])
      ),
      tibble::tibble(
        Gene = rep(topN_genes, each = num_runs),
        Method = "Permutation",
        Value = as.vector(B_perm[topN_genes, , drop = FALSE])
      )
    )
    if (uses_shap_for_type(model_type, method_weights)) {
      df_list_box[[length(df_list_box) + 1]] <- tibble::tibble(
        Gene = rep(topN_genes, each = num_runs),
        Method = "SHAP",
        Value = as.vector(C_shap[topN_genes, , drop = FALSE])
      )
    }
    long_df <- dplyr::bind_rows(df_list_box)
    
    p_box <- ggplot2::ggplot(long_df, ggplot2::aes(x = Gene, y = Value, fill = Method)) +
      ggplot2::geom_boxplot() +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Top ", topN, " genes importance distributions"),
                    y = "Importance (per-run)", x = "Gene") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    f_box <- file.path(plots_dir, "FI_Boxplot_Top20.pdf")
    ggplot2::ggsave(f_box, p_box, width = 9, height = 7)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_box))
  }
  
  # Bar (Top Composite) with 95% CI
  top_genes_plot <- importance_df |>
    dplyr::arrange(dplyr::desc(Composite))
  top_genes_plot <- dplyr::filter(top_genes_plot, Gene %in% plot_set)
  top_genes_plot <- dplyr::slice_head(top_genes_plot, n = num_coregene)
  if (nrow(top_genes_plot)) {
    p_bar <- ggplot2::ggplot(top_genes_plot,
                             ggplot2::aes(x = stats::reorder(Gene, Composite), y = Composite)) +
      ggplot2::geom_col(fill = "#0072B2", alpha = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = Composite_lo, ymax = Composite_hi),
        width = 0.4, color = "#D55E00"
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0(num_coregene, " Genes by Composite Importance"),
                    x = "Gene", y = "Composite (mean with 95% CI)") +
      ggplot2::theme_minimal()
    f_bar <- file.path(plots_dir, "FI_Bar_TopComposite.pdf")
    ggplot2::ggsave(f_bar, p_bar, width = 8, height = 6)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_bar))
  }
  
  # Density (all)
  # Build combined importance long df, dropping SHAP if not used
  include_shap <- isTRUE(uses_shap_for_type(model_type, method_weights))
  A_use <- A_model[plot_set, , drop = FALSE]
  B_use <- B_perm[plot_set, , drop = FALSE]
  
  df_list <- list(
    data.frame(Method = "Model", Value = as.vector(A_use)),
    data.frame(Method = "Permutation", Value = as.vector(B_use))
  )
  if (include_shap) {
    C_use <- C_shap[plot_set, , drop = FALSE]
    df_list[[length(df_list) + 1]] <- data.frame(Method = "SHAP", Value = as.vector(C_use))
  }
  combined_imp <- dplyr::bind_rows(df_list)
  
  p_den <- ggplot2::ggplot(combined_imp, ggplot2::aes(x = Value, fill = Method)) +
    ggplot2::geom_density(alpha = 0.6, color = NA) +
    ggplot2::geom_vline(
      data = combined_imp %>%
        dplyr::group_by(Method) %>%
        dplyr::summarise(Median = stats::median(Value, na.rm = TRUE)),
      ggplot2::aes(xintercept = Median, color = Method), linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~Method, scales = "free") +
    ggplot2::labs(title = "Feature Importance Distribution Comparison",
                  x = "Importance Score", y = "Density") +
    ggplot2::theme_bw()
  f_den <- file.path(plots_dir, "FI_Density_AllMethods.pdf")
  ggplot2::ggsave(f_den, p_den, width = 9, height = 6)
  message(sprintf("[%s] Saved plot to: %s", .ts(), f_den))
  
  # UMAP (consistent with plot_set)
  final_importance <- importance_df |>
    dplyr::arrange(dplyr::desc(Composite))
  final_top <- intersect(final_importance$Gene, plot_set)
  final_top <- head(final_top, n_likes_eff)
  set.seed(seed)
  if (length(final_top) >= 2) {
    um <- umap::umap(train_exp[, final_top, drop = FALSE])
    umap_df <- data.frame(
      UMAP1 = um$layout[, 1],
      UMAP2 = um$layout[, 2],
      Label = factor(train_y, labels = c("Control", "Case"))
    )
    p_umap <- ggplot2::ggplot(umap_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = Label)) +
      ggplot2::geom_point(alpha = 0.7, size = 2.5) +
      ggplot2::scale_color_manual(values = c("#377eb8", "#e41a1c")) +
      ggplot2::labs(title = "UMAP Projection using top Signature Genes", color = "Group") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    f_umap <- file.path(plots_dir, "UMAP_TopSignatureGenes.pdf")
    ggplot2::ggsave(f_umap, p_umap, width = 6, height = 5)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_umap))
  }
  
  # Stability heatmap (already aligned via plot_set-aware top_df)
  if (nrow(top_df)) {
    p_stab <- ggplot2::ggplot(top_df, ggplot2::aes(x = Gene, y = "Runs", fill = Frequency)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = Frequency), color = "black", size = 4) +
      ggplot2::scale_fill_gradient(low = "#E6F5FF", high = "#1F77B4", name = "Frequency") +
      ggplot2::labs(title = paste0("Top Gene Stability Across ", num_runs, " Runs"),
                    x = "Gene", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    f_stab <- file.path(plots_dir, "Stability_TopGenes_Heatmap.pdf")
    ggplot2::ggsave(f_stab, p_stab, width = 9, height = 6)
    message(sprintf("[%s] Saved plot to: %s", .ts(), f_stab))
  }
  
  
  ### formula construction for top n_likes_eff genes
  feat_names <- final_top[1:min(n_likes_eff, length(final_top))]
  x_all <- as.matrix(train_exp[, feat_names, drop = FALSE])
  y_all <- as.numeric(train_y)
  
  build_linear_formula <- function(intercept, coefs, names_vec, prefix = "logit(p) = ") {
    fs <- paste0(prefix, round(intercept, 6))
    for (i in seq_along(coefs)) {
      sgn <- ifelse(coefs[i] >= 0, " + ", " - ")
      fs <- paste0(fs, sgn, abs(round(coefs[i], 6)), " * ", names_vec[i])
    }
    fs
  }
  
  build_quadratic_formula <- function(intercept, lin_coefs, quad_mat, names_vec,
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
      for (i in seq_len(p)) for (j in seq((i + 1), p)) {
        qij <- quad_mat[i, j] + quad_mat[j, i]
        if (is.finite(qij) && abs(qij) >= tol) {
          sgn <- ifelse(qij >= 0, " + ", " - ")
          fs <- paste0(fs, sgn, abs(round(qij, 6)), " * ", names_vec[i], " * ", names_vec[j])
        }
      }
    }
    fs
  }
  
  .fit_logit_or_ridge <- function(x, y,
                                  alpha_ridge = 0,        # ridge = 0
                                  firth_ci = FALSE,
                                  firth_maxit = 1000,
                                  firth_pl_maxit = 1000,
                                  # FIX: 可调阈值
                                  coef_guard = 20,        # |beta| 超过就触发回退
                                  pred_guard_pct = 0.99,  # 极端预测比例阈值
                                  pred_eps = 1e-6) {
    
    x <- as.matrix(x)
    if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
    
    # FIX: 规范 y 为 {0,1}
    if (is.factor(y)) {
      lev <- levels(y)
      if (all(lev %in% c("0","1"))) {
        y <- as.integer(as.character(y))         # "0"/"1" -> 0/1
      } else {
        y <- as.integer(y == lev[2L])            # 以第二个水平为阳性 -> 0/1
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
          sep_flag <<- TRUE; invokeRestart("muffleWarning")
        } else if (grepl("algorithm did not converge", msg)) {
          nonconv_flag <<- TRUE; invokeRestart("muffleWarning")
        }
      }
    )
    
    coef_vec <- tryCatch(stats::coef(glm_fit), error = function(e) rep(NA_real_, ncol(x) + 1))
    bad_coef <- any(!is.finite(coef_vec))
    
    # ---- 定义一个统一的回退函数 ----
    .fallback_ridge <- function() {
      if (!requireNamespace("glmnet", quietly = TRUE)) return(NULL)
      cvfit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = alpha_ridge, standardize = TRUE)
      beta  <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
      intercept <- as.numeric(beta[1, 1])
      coefs <- as.numeric(beta[-1, 1]); names(coefs) <- rownames(beta)[-1]
      coefs <- coefs[colnames(x)]; names(coefs) <- colnames(x)
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
      if (!requireNamespace("logistf", quietly = TRUE)) return(NULL)
      lf <- tryCatch(
        logistf::logistf(
          y ~ .,
          data = logistic_data,
          control   = logistf::logistf.control(maxit = firth_maxit),
          pl        = isTRUE(firth_ci),
          plcontrol = if (isTRUE(firth_ci)) logistf::logistpl.control(maxit = firth_pl_maxit) else NULL
        ),
        error = function(e) NULL
      )
      if (is.null(lf)) return(NULL)
      cf <- tryCatch(stats::coef(lf), error = function(e) NULL)
      if (is.null(cf)) cf <- tryCatch(lf$coefficients, error = function(e) NULL)
      if (is.null(cf)) return(NULL)
      varmat <- tryCatch(lf$var, error = function(e) NULL)
      se <- if (!is.null(varmat) && is.matrix(varmat) &&
                nrow(varmat) == length(cf) && ncol(varmat) == length(cf)) {
        sqrt(pmax(diag(varmat), 0))
      } else rep(NA_real_, length(cf))
      zval <- ifelse(is.finite(se) & se > 0, cf / se, NA_real_)
      pval <- ifelse(is.finite(zval), 2 * stats::pnorm(-abs(zval)), NA_real_)
      coef_summary <- cbind(
        Estimate     = as.numeric(cf),
        `Std. Error` = as.numeric(se),
        `z value`    = as.numeric(zval),
        `Pr(>|z|)`   = as.numeric(pval)
      )
      rownames(coef_summary) <- names(cf)
      
      intercept <- unname(cf["(Intercept)"]); if (!is.finite(intercept)) intercept <- 0
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
      if (!is.null(firth_res)) return(firth_res)
      ridge_res <- .fallback_ridge()
      if (!is.null(ridge_res)) return(ridge_res)
      # 实在不行才兜底返回 glm 结果
      intercept <- coef_vec[1]
      coefs <- coef_vec[-1]; coefs[!is.finite(coefs)] <- 0
      names(coefs) <- colnames(x)
      coef_summary <- tryCatch(summary(glm_fit)$coefficients, error = function(e) {
        cm <- cbind(
          Estimate     = c(intercept, coefs),
          `Std. Error` = NA_real_,
          `z value`    = NA_real_,
          `Pr(>|z|)`   = NA_real_
        )
        rownames(cm) <- c("(Intercept)", colnames(x)); cm
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
    coefs <- coef_vec[-1]; names(coefs) <- colnames(x)
    p_hat <- tryCatch(stats::predict(glm_fit, type = "response"), error = function(e) NULL)
    extreme_pred <- !is.null(p_hat) && mean(p_hat < pred_eps | p_hat > 1 - pred_eps) > pred_guard_pct
    too_large <- (is.finite(intercept) && abs(intercept) > coef_guard) ||
      any(is.finite(coefs) & abs(coefs) > coef_guard)
    
    if (extreme_pred || too_large) {
      firth_res <- .fallback_firth()
      if (!is.null(firth_res)) return(firth_res)
      ridge_res <- .fallback_ridge()
      if (!is.null(ridge_res)) return(ridge_res)
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
  
  # helpers for QDA and Gaussian NB native formulas
  .qda_to_logit <- function(X, y) {
    X <- as.matrix(X)
    y <- .y_to_fac(y)                 
    cls <- levels(y)
    if (length(cls) != 2) stop("QDA requires binary y.")
    x0 <- X[y == cls[1], , drop = FALSE]
    x1 <- X[y == cls[2], , drop = FALSE]
    mu0 <- colMeans(x0); mu1 <- colMeans(x1)
    S0 <- stats::cov(x0); S1 <- stats::cov(x1)
    n0 <- nrow(x0); n1 <- nrow(x1)
    pi0 <- n0 / (n0 + n1); pi1 <- 1 - pi0
    
    # regularize if needed
    reg <- function(S) {
      if (!all(is.finite(S)) || det(S) == 0) S <- S + diag(1e-6, ncol(S))
      S
    }
    S0 <- reg(S0); S1 <- reg(S1)
    
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
    mu0 <- colMeans(x0); mu1 <- colMeans(x1)
    v0 <- apply(x0, 2, stats::var); v1 <- apply(x1, 2, stats::var)
    v0[!is.finite(v0) | v0 <= 0] <- 1e-6
    v1[!is.finite(v1) | v1 <= 0] <- 1e-6
    n0 <- nrow(x0); n1 <- nrow(x1)
    pi0 <- n0 / (n0 + n1); pi1 <- 1 - pi0
    
    # logit = c + sum(b_i x_i) + sum(a_i x_i^2), a_i only on diagonal (no cross terms)
    a <- 0.5 * (1 / v0 - 1 / v1)
    b <- (mu1 / v1) - (mu0 / v0)
    c <- -0.5 * (sum((mu1^2) / v1) - sum((mu0^2) / v0)) -
      0.5 * sum(log(v1 / v0)) + log(pi1 / pi0)
    
    list(intercept = c, linear = b, quad_diag = a)
  }
  
  formula_source <- NULL
  formula_str <- NULL
  coef_summary <- NULL
  
  if (identical(model_type, "glmnet")) {
    a_use <- parsed$alpha
    if (is.null(a_use) || is.na(a_use)) a_use <- 1
    y_all_num <- .y_to_num(y_all)
    cvfit <- glmnet::cv.glmnet(x_all, y_all_num, family = "binomial",
                           alpha = a_use, standardize = TRUE)
    beta  <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
    intercept <- as.numeric(beta[1, 1])
    
    # 先用beta自带的行名（去掉截距），再与x_all列名对齐，避免长度不一致
    raw_names <- rownames(beta)[-1]
    raw_coefs <- as.numeric(beta[-1, 1])
    names(raw_coefs) <- raw_names
    
    # 按 x_all 的列顺序对齐（有些列可能被惩罚到0，但仍保留在coef里）
    aligned <- setNames(rep(0, ncol(x_all)), colnames(x_all))
    common  <- intersect(names(raw_coefs), colnames(x_all))
    aligned[common] <- raw_coefs[common]
    coefs <- .safe_set_names(aligned, colnames(x_all))
    
    coef_summary <- cbind(
      Estimate     = c(intercept, as.numeric(coefs)),
      `Std. Error` = NA_real_,
      `z value`    = NA_real_,
      `Pr(>|z|)`   = NA_real_
    )
    rownames(coef_summary) <- c("(Intercept)", colnames(x_all))
    
    formula_str <- build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
    cat(formula_str, "\n")
    formula_source <- "native_glmnet"
  } else if (identical(model_type, "lr") || identical(model_type, "stepwise_lr")) {
    # native logistic on selected features (with robust fallback)
    res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
    intercept <- res$intercept; coefs <- res$coefs; coef_summary <- res$coef_summary
    feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
    names(coefs) <- feat_names_out
    rownames(coef_summary) <- c("(Intercept)", names(coefs))
    formula_str <- build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
    cat(formula_str, "\n")
    formula_source <- paste0(model_type, "_native_logit")
    
  } else if (identical(model_type, "lda")) {
    if (length(unique(y_all)) != 2) {
      # fallback: 用逻辑回归/岭回归
      res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
      intercept   <- res$intercept
      coefs       <- .safe_set_names(res$coefs, colnames(x_all))
      coef_summary <- cbind(
        Estimate     = c(intercept, as.numeric(coefs)),
        `Std. Error` = NA_real_,
        `z value`    = NA_real_,
        `Pr(>|z|)`   = NA_real_
      )
      rownames(coef_summary) <- c("(Intercept)", colnames(x_all))
      formula_str <- build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
      cat(formula_str, "\n")
      formula_source <- paste0("lda_", res$source)
      
    } else {
      # 原生 LDA → 等价的logit系数
      cls <- sort(unique(y_all))
      x0 <- x_all[y_all == cls[1], , drop = FALSE]
      x1 <- x_all[y_all == cls[2], , drop = FALSE]
      mu0 <- colMeans(x0); mu1 <- colMeans(x1)
      S0 <- stats::cov(x0); S1 <- stats::cov(x1)
      n0 <- nrow(x0); n1 <- nrow(x1)
      Sp <- ((n0 - 1) * S0 + (n1 - 1) * S1) / (n0 + n1 - 2)
      if (any(!is.finite(Sp)) || det(Sp) == 0) {
        Sp <- Sp + diag(1e-6, ncol(Sp))
      }
      beta_vec <- tryCatch(solve(Sp, (mu1 - mu0)),
                           error = function(e) { MASS::ginv(Sp) %*% (mu1 - mu0) })
      beta_vec <- as.numeric(beta_vec)
      coefs <- .safe_set_names(beta_vec, colnames(x_all))
      
      pi1 <- n1 / (n0 + n1); pi0 <- 1 - pi1
      intercept <- as.numeric(-0.5 * crossprod((mu1 + mu0), as.numeric(coefs)) + log(pi1 / pi0))
      
      coef_summary <- cbind(
        Estimate     = c(intercept, as.numeric(coefs)),
        `Std. Error` = NA_real_,
        `z value`    = NA_real_,
        `Pr(>|z|)`   = NA_real_
      )
      rownames(coef_summary) <- c("(Intercept)", colnames(x_all))
      
      formula_str <- build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
      cat(formula_str, "\n")
      formula_source <- "native_lda_logit"
    }
    
  } else if (identical(model_type, "qda")) {
    q <- .qda_to_logit(x_all, y_all)
    intercept <- q$intercept
    lin <- q$linear; names(lin) <- colnames(x_all)
    quad <- q$quad; rownames(quad) <- colnames(x_all); colnames(quad) <- colnames(x_all)
    coef_summary <- cbind(Estimate = c(intercept, lin),
                          `Std. Error` = NA_real_,
                          `z value` = NA_real_,
                          `Pr(>|z|)` = NA_real_)
    rownames(coef_summary) <- c("(Intercept)", names(lin))
    formula_str <- build_quadratic_formula(intercept, lin, quad, colnames(x_all), prefix = "logit(p) = ")
    cat(formula_str, "\n")
    formula_source <- "native_qda_quadratic"
    
  } else if (identical(model_type, "nb")) {
    nbp <- .nb_gaussian_to_logit(x_all, y_all)
    intercept <- nbp$intercept
    lin <- nbp$linear; names(lin) <- colnames(x_all)
    quad_diag <- nbp$quad_diag
    Q <- diag(as.numeric(quad_diag), nrow = ncol(x_all))
    rownames(Q) <- colnames(x_all); colnames(Q) <- colnames(x_all)
    coef_summary <- cbind(Estimate = c(intercept, lin),
                          `Std. Error` = NA_real_,
                          `z value` = NA_real_,
                          `Pr(>|z|)` = NA_real_)
    rownames(coef_summary) <- c("(Intercept)", names(lin))
    formula_str <- build_quadratic_formula(intercept, lin, Q, colnames(x_all), prefix = "logit(p) = ")
    cat(formula_str, "\n")
    formula_source <- "native_nb_gaussian_quadratic"
    
  } else if (identical(model_type, "svm_linear")) {
    svm_fit <- e1071::svm(x = x_all, y = factor(y_all),
                          kernel = "linear", probability = TRUE,
                          cost = svm_cost, scale = TRUE)
    if (is.null(svm_fit$coefs) || is.null(svm_fit$SV)) {
      res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
      intercept <- res$intercept; coefs <- res$coefs; coef_summary <- res$coef_summary
      feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
      names(coefs) <- feat_names_out
      formula_str <- build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
      cat(formula_str, "\n")
      formula_source <- "svm_linear_fallback_logit"
    } else {
      wv <- as.numeric(crossprod(svm_fit$coefs, svm_fit$SV))
      b <- -as.numeric(svm_fit$rho)
      names(wv) <- colnames(x_all)
      margin_str <- build_linear_formula(b, wv, names(wv), prefix = "margin = ")
      cat(margin_str, "\n")
      coef_summary <- cbind(Estimate = c(b, wv),
                            `Std. Error` = NA_real_,
                            `z value` = NA_real_,
                            `Pr(>|z|)` = NA_real_)
      rownames(coef_summary) <- c("(Intercept)", names(wv))
      formula_str <- margin_str
      formula_source <- "native_svm_linear_margin"
    }
    
  } else if (model_type %in% c("mlp","rf","xgb","gbm","svm_rbf","dt","tree","knn")) {
    # black-box models: interpretable linear surrogate
    res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
    intercept <- res$intercept; coefs <- res$coefs; coef_summary <- res$coef_summary
    feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
    names(coefs) <- feat_names_out
    formula_str <- build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
    cat(formula_str, "\n")
    formula_source <- paste0(model_type, "_surrogate_", res$source)
    
  } else {
    # unknown: fallback to interpretable logistic
    res <- .fit_logit_or_ridge(x_all, y_all, alpha_ridge = 0)
    intercept <- res$intercept; coefs <- res$coefs; coef_summary <- res$coef_summary
    feat_names_out <- setdiff(rownames(coef_summary), "(Intercept)")
    names(coefs) <- feat_names_out
    formula_str <- build_linear_formula(intercept, coefs, names(coefs), prefix = "logit(p) = ")
    cat(formula_str, "\n")
    formula_source <- "unknown_surrogate_logit"
  }
  
  
  
  # ---- persist formula to text (snapshot + rolling log) ----
  ts_disp <- if (exists(".ts")) .ts() else format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ts_file <- gsub("[: ]", "-", ts_disp)
  
  file_ts <- file.path(plots_dir, sprintf("fh_formula_%s.txt", ts_file))
  log_file <- file.path(plots_dir, "fh_formula_log.txt")
  
  lines <- c(
    sprintf("[fh_hunter] timestamp: %s", ts_disp),
    sprintf("[fh_hunter] model_type: %s", if (exists("model_type")) as.character(model_type) else NA_character_),
    sprintf("[fh_hunter] formula_source: %s", if (exists("formula_source")) as.character(formula_source) else NA_character_),
    sprintf("[fh_hunter] n_features: %s", if (exists("coef_summary")) nrow(coef_summary) - 1L else NA_integer_),
    sprintf("[fh_hunter] formula: %s", if (exists("formula_str")) formula_str else NA_character_)
  )
  
  writeLines(lines, file_ts)  # snapshot file
  write(paste(lines, collapse = "\n"), file = log_file, append = TRUE)
  write("\n---\n", file = log_file, append = TRUE)
  
  message(sprintf("[%s] Saved formula to: %s", if (exists(".ts")) .ts() else ts_disp, file_ts))
  message(sprintf("[%s] Appended formula to: %s", if (exists(".ts")) .ts() else ts_disp, log_file))
  
  # ---- wrap up and return ----
  invisible(list(
    params = list(
      top_models_csv = top_models_csv,
      score_name = namesS[score_index],
      picked_row = pick_index,
      parsed = parsed,
      seed = seed,
      out_dir = out_dir,
      method_weights = method_weights,
      ci_level = ci_level,
      strict_min_freq = strict_min_freq,
      strict_min_effect = strict_min_effect,
      strict_ci_gate = strict_ci_gate,
      strict_knee = strict_knee,
      strict_nmax = strict_nmax,
      apply_selection_to_plots = apply_selection_to_plots,
      formula_source = formula_source
    ),
    importance_df = importance_df,
    top_list = top_list,
    final_top = final_top,
    composite_mat = composite_mat,
    glm_summary = list(
      coefficients = coef_summary,
      formula = formula_str
    )
  ))
}