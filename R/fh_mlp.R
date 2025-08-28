#' @title MLP Grid Search Runner (regularized & stabilized)
#'
#' @description
#' Grid search MLP hyperparameters (learning rate, batch size, dropout, etc.),
#' pick thresholds on validation set B using an automatic method selector
#' (merged with fixed cutoff values), and evaluate F1 / Accuracy / Recall on
#' datasets A/B/C/D... Results are stored into
#' \code{all_result_FS}, \code{all_result_acc}, \code{all_result_recall},
#' \code{all_result_summary}.
#'
#' Includes stronger regularization and training stabilization:
#' - Standardize features with train-set mean/sd (shared to val/test)
#' - Batch Normalization (toggle via \code{use_batchnorm})
#' - L2 kernel regularization (\code{l2})
#' - Gaussian input noise (\code{gaussian_noise_sd})
#' - ReduceLROnPlateau + EarlyStopping
#'
#' @param train_exp        Training expression matrix
#' @param train_labels     Training labels (0/1 or binary factor)
#' @param test_exp         List: first element = validation set B, others = test sets C/D/…
#' @param epochselecti     Number of training epochs
#' @param batchsize_grid   Batch size grid; auto-generated if "auto"
#' @param learningrate_grid Learning rate grid
#' @param dropoutrate_grid  Dropout rate grid
#' @param cutoff_grid      Cutoff grid (merged with auto thresholds)
#' @param imbalance_thresh Imbalance threshold for deciding class weighting
#' @param auto_th_method   Thresholding method for the auto threshold on B:
#'   one of \code{"youden"}, \code{"f1"}, \code{"auto"} (default).
#'   If "auto", the method is decided by \code{decide_threshold_method()}.
#' @param early_patience   Patience for early stopping (default 8)
#' @param seed             Random seed
#' @param labels_list      Label list (for \code{get_labels_for_matrix()})
#' @param collector        (optional) Environment collector; if supplied, stores results into it
#'
#' @param standardize      Logical; standardize by train-set mean/sd (default TRUE)
#' @param hidden_units     Integer vector; sizes of hidden layers (default c(32,16,8))
#' @param activation       Activation for hidden layers (default "relu")
#' @param use_batchnorm    Logical; add BatchNorm after each Dense (default TRUE)
#' @param l2               Numeric; L2 kernel regularization (default 1e-4)
#' @param gaussian_noise_sd Numeric; std of input GaussianNoise (default 0)
#' @param min_lr           Numeric; floor learning rate for ReduceLROnPlateau (default 1e-5)
#' @param plateau_factor   Numeric; LR decay factor on plateau (default 0.5)
#' @param plateau_patience Integer; epochs to wait before LR reduce (default 4)
#'
#' @return Invisible NULL or a list of result containers if no collector is provided.
#'
#' @details
#' - Validation set B is always \code{test_exp[[1]]}.
#' - Auto thresholding calls \code{decide_threshold_method()} to choose between
#'   "youden" and "f1" when \code{auto_th_method="auto"}, then applies
#'   \code{choose_threshold()} on B to obtain the numeric cutoff.
#' - Each model key includes cutoff, learning rate, batch size, epoch, and dropout.
#'
#' @import keras3
#' @importFrom tensorflow tf
#' @export
fh_mlp <- function(
    train_exp,
    train_labels,
    test_exp,
    epochselecti = 50,
    batchsize_grid = "auto",
    learningrate_grid = c(1e-3, 5e-3, 1e-2, 5e-2),
    dropoutrate_grid = c(0.25, 0.5, 0.75),
    cutoff_grid = c(0.25, 0.5, 0.75),
    imbalance_thresh = 0.35,
    auto_th_method = "auto",
    early_patience = 8,
    seed = seed,
    labels_list = labels_list,
    collector = NULL,
    # New knobs for generalization
    standardize = FALSE,
    hidden_units = c(32L, 16L, 8L),
    activation = "relu",
    use_batchnorm = TRUE,
    l2 = 1e-4,
    gaussian_noise_sd = 0,
    min_lr = 1e-5,
    plateau_factor = 0.5,
    plateau_patience = 4
) {
  # ---- 1) Inputs & labels ----
  train_exp <- safe_matrix(train_exp)
  train_labels <- as_binary01(train_labels)
  input_dim <- ncol(train_exp)
  stopifnot(input_dim > 0L, nrow(train_exp) > 0L)
  stopifnot(nrow(train_exp) == length(train_labels))
  
  # Auto batchsize grid
  if (identical(batchsize_grid, "auto")) {
    step <- max(1, floor(nrow(train_exp) / 50))
    batchsize_grid <- seq(5, nrow(train_exp), by = step)
  }
  
  set.seed(seed)
  if (reticulate::py_module_available("tensorflow")) {
    tensorflow::set_random_seed(seed)
    tensorflow::tf$random$set_seed(as.integer(seed))
  }
  
  # ---- 2) Validation/Test sets ----
  if (length(test_exp) < 1L)
    stop("test_exp must include at least one dataset; test_exp[[1]] is used as validation set B.")
  x_val <- safe_matrix(test_exp[[1]])
  y_val <- as_binary01(get_labels_for_matrix(x_val, labels_list))
  
  test_exp_rest <- if (length(test_exp) >= 2L) test_exp[-1] else list()
  n_test <- length(test_exp_rest)
  
  # ---- 3) Optional standardization ----
  if (isTRUE(standardize)) {
    tr_mean <- colMeans(train_exp, na.rm = TRUE)
    tr_sd   <- apply(train_exp, 2, stats::sd)
    tr_sd[!is.finite(tr_sd) | tr_sd == 0] <- 1
    scale_mat <- function(M) {
      M <- safe_matrix(M)
      sweep(sweep(M, 2, tr_mean, "-"), 2, tr_sd, "/")
    }
    train_exp <- scale_mat(train_exp)
    x_val     <- scale_mat(x_val)
    if (n_test) test_exp_rest <- lapply(test_exp_rest, scale_mat)
  }
  
  # ---- 4) Imbalance handling ----
  use_imbalance <- need_imbalance_handling(train_labels, thresh = imbalance_thresh)
  cw <- if (use_imbalance) class_weights_balanced(train_labels) else NULL
  init_bias <- if (use_imbalance) init_bias_from_labels(train_labels) else 0
  
  auto_th_method <- match.arg(auto_th_method, c("youden", "f1", "auto"))
  
  # Helper: build model with current (dr, lr, etc.)
  build_mlp <- function(input_dim, dr, lr) {
    reg <- if (l2 > 0) keras3::regularizer_l2(l2) else NULL
    mdl <- keras3::keras_model_sequential()
    
    mdl$add(keras3::layer_input(shape = c(input_dim)))
    
    if (is.finite(gaussian_noise_sd) && gaussian_noise_sd > 0) {
      mdl$add(keras3::layer_gaussian_noise(stddev = gaussian_noise_sd))
    }
    
    # hidden blocks
    for (u in hidden_units) {
      mdl$add(keras3::layer_dense(units = as.integer(u),
                                  activation = activation,
                                  kernel_regularizer = reg))
      if (isTRUE(use_batchnorm)) {
        mdl$add(keras3::layer_batch_normalization())
      }
      mdl$add(keras3::layer_dropout(rate = dr))
    }
    
    # output
    out_bias <- if (use_imbalance) keras3::initializer_constant(init_bias) else "zeros"
    mdl$add(keras3::layer_dense(units = 1, activation = "sigmoid",
                                bias_initializer = out_bias))
    
    opt <- keras3::optimizer_adam(learning_rate = lr)
    mdl |>
      keras3::compile(
        optimizer = opt,
        loss = "binary_crossentropy",
        metrics = c(keras3::metric_auc(name = "auc_roc"),
                    keras3::metric_auc(curve = "PR", name = "prc"))
      )
    mdl
  }
  
  # Assemble datasets list (names must match downstream)
  mk_sets_list <- function() {
    tests_list <- if (n_test) {
      labs <- paste0("Dataset", LETTERS[3:(2 + n_test)], "(Test)")
      lapply(seq_len(n_test), function(i) {
        Xi <- safe_matrix(test_exp_rest[[i]])
        yi <- as_binary01(get_labels_for_matrix(Xi, labels_list))
        list(x = Xi, y = yi, name = labs[i])
      })
    } else list()
    c(
      list(list(x = train_exp, y = train_labels, name = "DatasetA(Train)")),
      list(list(x = x_val,    y = y_val,        name = "DatasetB(Val)")),
      tests_list
    )
  }
  
  # ---- 5) Grid search ----
  for (bs in batchsize_grid) {
    if (!is.finite(bs) || bs < 1 || bs > nrow(train_exp)) next
    
    for (lr in learningrate_grid) {
      for (dr in dropoutrate_grid) {
        
        model <- build_mlp(input_dim, dr = dr, lr = lr)
        
        # callbacks: EarlyStopping + ReduceLROnPlateau
        cbs <- list(
          keras3::callback_early_stopping(
            monitor = "val_prc", mode = "max",
            patience = early_patience, restore_best_weights = TRUE
          ),
          keras3::callback_reduce_lr_on_plateau(
            monitor = "val_prc", mode = "max",
            factor = plateau_factor, patience = plateau_patience,
            min_lr = min_lr, verbose = 0
          )
        )
        
        hist <- tryCatch(
          keras3::fit(
            model,
            x = train_exp, y = train_labels,
            epochs = epochselecti, batch_size = as.integer(bs),
            validation_data = list(x_val, y_val),
            class_weight = cw,
            shuffle = TRUE,
            callbacks = cbs,
            verbose = 0
          ),
          error = function(e) {
            message(sprintf(
              "[MLP][skip] fit error @ lr=%.4g, bs=%d, dr=%.2f | %s",
              lr, bs, dr, conditionMessage(e)
            ))
            NULL
          }
        )
        print(model)
        if (is.null(hist)) next
        
        # ---- 6) Auto threshold selection on B ----
        prob_val <- as.numeric(model$predict(x_val))
        if (!all(is.finite(prob_val))) {
          message("[MLP] non-finite probabilities, skip this combination.")
          next
        }
        method_chosen <- if (identical(auto_th_method, "auto")) {
          decide_threshold_method(prob_val, y_val)
        } else {
          auto_th_method
        }
        th_auto <- choose_threshold(prob_val, y_val, method = method_chosen)
        
        # ---- 7) Eval on A/B/C... at fixed and auto thresholds ----
        sets_list <- mk_sets_list()
        names_now <- vapply(sets_list, `[[`, "", "name")
        
        th_vec <- c(cutoff_grid, th_auto)
        th_lab <- c(paste0("cutoff:", cutoff_grid),
                    sprintf("cutoff:auto(%.3f,%s)", th_auto, method_chosen))
        
        for (k in seq_along(th_vec)) {
          th <- th_vec[k]; th_name <- th_lab[k]
          
          ev <- lapply(sets_list, function(d) {
            p  <- as.numeric(model$predict(d$x))
            pr <- as.integer(p >= th)
            list(
              raw = data.frame(predict_p = p, predict_result = pr, real_label = d$y),
              met = metrics_from_preds(pr, d$y)
            )
          })
          mets <- do.call(cbind, lapply(ev, `[[`, "met"))
          
          if (ncol(mets) != (n_test + 2L)) {
            stop(sprintf("Internal consistency error: expected %d columns, got %d. names: %s",
                         n_test + 2L, ncol(mets),
                         paste(names_now, collapse = " | ")))
          }
          
          colnames(mets) <- names_now
          result_FS     <- as.numeric(mets["F1",  ]);  names(result_FS)     <- names_now
          result_acc    <- as.numeric(mets["ACC", ]);  names(result_acc)    <- names_now
          result_recall <- as.numeric(mets["RECALL",]);names(result_recall) <- names_now
          
          key <- paste0("NN-MLP (", th_name,
                        ", lr:", lr,
                        ", bs:", bs,
                        ", ep:", epochselecti,
                        ", dropout:", dr, ")")
          
          if (is.environment(collector)) {
            if (is.null(collector$all_result_acc))     collector$all_result_acc     <- list()
            if (is.null(collector$all_result_recall))  collector$all_result_recall  <- list()
            if (is.null(collector$all_result_FS))      collector$all_result_FS      <- list()
            if (is.null(collector$all_result_summary)) collector$all_result_summary <- list()
            
            collector$all_result_acc[[key]]     <- result_acc
            collector$all_result_recall[[key]]  <- result_recall
            collector$all_result_FS[[key]]      <- result_FS
            collector$all_result_summary[[key]] <- lapply(ev, `[[`, "raw")
          } else {
            if (!exists("acc_local", inherits = FALSE)) {
              acc_local     <- list(); recall_local  <- list()
              fs_local      <- list(); summary_local <- list()
            }
            acc_local[[key]]     <- result_acc
            recall_local[[key]]  <- result_recall
            fs_local[[key]]      <- result_FS
            summary_local[[key]] <- lapply(ev, `[[`, "raw")
            
            if (!is.environment(collector)) {
              return(invisible(list(
                all_result_acc     = acc_local,
                all_result_recall  = recall_local,
                all_result_FS      = fs_local,
                all_result_summary = summary_local
              )))
            } else {
              return(invisible(NULL))
            }
          }
        }
      }
    }
  }
  invisible(NULL)
}