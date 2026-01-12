#' Internal training utilities for FeatureHunter
#' @keywords internal
#' @importFrom stats predict sd median mad quantile reorder glm binomial
#' @importFrom utils txtProgressBar setTxtProgressBar

.train_single_model <- function(model_type, train_exp, train_y, seed, fit_params, run_id = 1) {
    # fit_params is a list containing merged parsed params and explicit arguments

    fit <- switch(model_type,
        mlp = {
            if (!requireNamespace("keras3", quietly = TRUE)) stop("Package 'keras3' needed for MLP model.")
            if (!requireNamespace("tensorflow", quietly = TRUE)) stop("Package 'tensorflow' needed for MLP model.")

            get_or <- function(name, default = NULL) {
                if (!is.null(fit_params[[name]])) fit_params[[name]] else default
            }

            seed_i <- as.integer(seed)
            set.seed(seed_i)
            if (reticulate::py_module_available("tensorflow")) {
                tensorflow::set_random_seed(seed_i)
                tensorflow::tf$random$set_seed(seed_i)
            }

            bn_on <- isTRUE(get_or("use_batchnorm", TRUE))
            l2_reg <- get_or("l2", 1e-4) # default from code
            if (is.null(l2_reg)) l2_reg <- 1e-4
            gnoise_sd <- get_or("gaussian_noise_sd", 0.0)
            if (is.null(gnoise_sd)) gnoise_sd <- 0.0
            act_hidden <- get_or("activation", "relu")
            if (is.null(act_hidden)) act_hidden <- "relu"

            min_lr_i <- get_or("min_lr", 1e-5)
            plat_fac <- get_or("plateau_factor", 0.5)
            plat_pat <- as.integer(get_or("plateau_patience", 6))
            es_pat <- as.integer(get_or("early_patience", 9))
            shuffle_i <- isTRUE(get_or("shuffle", TRUE))

            units_vec <- get_or("hidden_units", c(32L, 16L, 8L))
            if (is.null(units_vec)) units_vec <- c(32L, 16L, 8L)
            units_vec <- as.integer(units_vec[is.finite(units_vec) & units_vec > 0])

            # params from parsed/overrides
            learningratei <- get_or("lr")
            if (is.null(learningratei)) stop("MLP lr missing")
            batch_sizei <- get_or("bs")
            if (is.null(batch_sizei)) stop("MLP bs missing")
            epochselecti <- get_or("ep")
            if (is.null(epochselecti)) stop("MLP ep missing")
            dropoutratei <- get_or("dropout")
            if (is.null(dropoutratei)) stop("MLP dropout missing")

            y01 <- as.integer(train_y)
            pos_rate <- mean(y01 == 1)
            neg_rate <- 1 - pos_rate
            imb_th <- get_or("imbalance_thresh", 0.35)
            if (is.null(imb_th)) imb_th <- 0.35
            use_imbalance <- (min(pos_rate, neg_rate) < imb_th)

            class_weighti <- if (use_imbalance) list(`0` = 0.5 / max(neg_rate, 1e-12), `1` = 0.5 / max(pos_rate, 1e-12)) else NULL
            out_bias <- if (use_imbalance) keras3::initializer_constant(log(pos_rate / max(neg_rate, 1e-12))) else "zeros"

            use_l2 <- is.numeric(l2_reg) && length(l2_reg) == 1 && is.finite(l2_reg) && l2_reg > 0
            reg <- if (use_l2) keras3::regularizer_l2(l2_reg) else NULL

            p <- ncol(train_exp)
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
                verbose = 0
            )

            # Prepare training plots info if run 1
            plots_data <- NULL
            if (run_id == 1) {
                loss_vec <- as.numeric(history$metrics$loss)
                acc_vec <- if (!is.null(history$metrics$accuracy)) {
                    as.numeric(history$metrics$accuracy)
                } else if (!is.null(history$metrics$binary_accuracy)) {
                    as.numeric(history$metrics$binary_accuracy)
                } else {
                    NULL
                }
                n_ep <- length(loss_vec)
                if (!is.null(acc_vec)) n_ep <- min(n_ep, length(acc_vec))
                if (n_ep > 0) {
                    plots_data <- data.frame(
                        Epoch    = seq_len(n_ep),
                        Loss     = loss_vec[seq_len(n_ep)],
                        Accuracy = if (is.null(acc_vec)) NA_real_ else acc_vec[seq_len(n_ep)]
                    )
                }
            }
            list(model = model, plots_data = plots_data)
        },
        rf = {
            if (!requireNamespace("ranger", quietly = TRUE)) stop("Package 'ranger' needed for random forest.")
            m <- ranger::ranger(
                dependent.variable.name = "y",
                data = data.frame(y = factor(train_y), train_exp),
                probability = TRUE,
                importance = "impurity",
                num.trees = as.integer(fit_params$num.trees),
                mtry = as.integer(fit_params$mtry),
                seed = seed
            )
            list(model = m)
        },
        glmnet = {
            # glmnet is in Imports, no check needed (or do it to be safe if moved to Suggests?)
            # Moved to Imports in my replace list? No, I moved glmnet to IMPORTS?
            # Wait, I originally moved glmnet to Suggests in my mind but let me check my last replace call.
            # Step 108: "Imports: ... glmnet ...". Ah. I kept glmnet in Imports!
            # OK, then no check needed.
            # Actually, let's verify. The replacement said "glmnet" in "Imports".
            # But wait, looking at description content in Step 108 args:
            # "Imports: ... glmnet ..."
            # So checks are OPTIONAL but good practice if I decide to move it later.
            # I'll skip check for glmnet if it is in Imports.
            m <- glmnet::cv.glmnet(
                x = train_exp,
                y = train_y,
                family = "binomial",
                alpha = if (!is.null(fit_params$alpha)) fit_params$alpha else 1,
                type.measure = "class",
                standardize = TRUE
            )
            list(model = m)
        },
        xgb = {
            if (!requireNamespace("xgboost", quietly = TRUE)) stop("Package 'xgboost' needed for this model_type.")
            dtrain <- xgboost::xgb.DMatrix(data = train_exp, label = train_y)
            m <- xgboost::xgboost(
                data = dtrain,
                objective = "binary:logistic",
                nrounds = as.integer(fit_params$nrounds),
                max_depth = as.integer(fit_params$max_depth),
                eta = fit_params$eta,
                subsample = 0.8, colsample_bytree = 0.8,
                verbose = 0
            )
            list(model = m)
        },
        svm_linear = {
            if (!requireNamespace("e1071", quietly = TRUE)) stop("Package 'e1071' needed for SVM.")
            m <- e1071::svm(x = train_exp, y = factor(train_y), kernel = "linear", probability = TRUE, cost = fit_params$cost, scale = TRUE)
            list(model = m)
        },
        svm_rbf = {
            if (!requireNamespace("e1071", quietly = TRUE)) stop("Package 'e1071' needed for SVM.")
            m <- e1071::svm(x = train_exp, y = factor(train_y), kernel = "radial", probability = TRUE, cost = fit_params$cost, gamma = fit_params$gamma, scale = TRUE)
            list(model = m)
        },
        lda = {
            list(model = MASS::lda(x = train_exp, grouping = factor(train_y)))
        },
        qda = {
            list(model = MASS::qda(x = train_exp, grouping = factor(train_y)))
        },
        nb = {
            if (!requireNamespace("e1071", quietly = TRUE)) stop("Package 'e1071' needed for Naive Bayes.")
            list(model = e1071::naiveBayes(x = train_exp, y = factor(train_y)))
        },
        gbm = {
            if (!requireNamespace("gbm", quietly = TRUE)) stop("Package 'gbm' needed for GBM.")
            m <- gbm::gbm(
                formula = y ~ .,
                data = data.frame(y = train_y, train_exp),
                distribution = "bernoulli",
                n.trees = 200L, # 可暴露为参数
                interaction.depth = 3L, # 可暴露为参数
                shrinkage = 0.05, # 可暴露为参数
                n.minobsinnode = 10L,
                bag.fraction = 0.8,
                train.fraction = 1.0,
                keep.data = FALSE,
                verbose = FALSE
            )
            list(model = m)
        },
        lr = {
            m <- stats::glm(y ~ .,
                data = data.frame(y = train_y, train_exp),
                family = stats::binomial(link = "logit")
            )
            list(model = m)
        },
        stepwise_lr = {
            base_fit <- stats::glm(y ~ .,
                data = data.frame(y = train_y, train_exp),
                family = stats::binomial(link = "logit")
            )
            m <- if (requireNamespace("MASS", quietly = TRUE)) {
                suppressWarnings(MASS::stepAIC(base_fit, trace = FALSE))
            } else {
                base_fit
            }
            list(model = m)
        },
        dt = {
            m <- if (requireNamespace("rpart", quietly = TRUE)) {
                rpart::rpart(y ~ .,
                    data = data.frame(y = factor(train_y), train_exp),
                    method = "class", parms = list(split = "gini")
                )
            } else {
                stop("rpart not installed but model_type=dt/tree was selected.")
            }
            list(model = m)
        },
        tree = {
            m <- if (requireNamespace("rpart", quietly = TRUE)) {
                rpart::rpart(y ~ .,
                    data = data.frame(y = factor(train_y), train_exp),
                    method = "class", parms = list(split = "gini")
                )
            } else {
                stop("rpart not installed but model_type=dt/tree was selected.")
            }
            list(model = m)
        },
        knn = {
            # “训练”就是保存标准化后的训练集与 y，以及 k
            k_knn <- max(3L, min(15L, as.integer(round(sqrt(nrow(train_exp))))))
            cen <- colMeans(train_exp)
            scl <- apply(train_exp, 2, stats::sd)
            scl[!is.finite(scl) | scl == 0] <- 1
            Xs <- scale(train_exp, center = cen, scale = scl)
            m <- list(
                kind = "knn", X_train = Xs, y = factor(train_y), k = k_knn,
                center = cen, scale = scl
            )
            list(model = m)
        }
    )

    fit
}
