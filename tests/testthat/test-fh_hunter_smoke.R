test_that("fh_hunter runs on iris dataset (RF)", {
    data(iris)
    # Binary classification: Virginica vs others
    y <- ifelse(iris$Species == "virginica", 1, 0)
    X <- iris[, 1:4]

    # Create a dummy leaderboard csv
    tmp_csv <- tempfile(fileext = ".csv")
    write.csv(data.frame(Model = c("rf"), F_score = c(0.9)), tmp_csv, row.names = FALSE)

    # Run fh_hunter
    # We use a temporary directory for output
    tmp_dir <- tempdir()

    expect_error(
        {
            res <- fh_hunter(
                train_exp = X,
                train_labels = y,
                nshow = 1,
                top_models_csv = tmp_csv,
                num_runs = 2,
                num_coregene = 2,
                n_likes = 2,
                n_interest = 2,
                out_dir = tmp_dir,
                rf_num_trees = 10, # smooth run
                score_index = 2 # F_score is 2nd col in our dummy
            )
        },
        NA
    ) # Expect no error

    expect_true(is.list(res))
    expect_true("importance_df" %in% names(res))
    expect_true("composite_mat" %in% names(res))
})
