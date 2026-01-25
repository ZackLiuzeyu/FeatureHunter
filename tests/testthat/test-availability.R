test_that("FeatureHunter package loads and key functions exist", {
    expect_true(exists("fh_hunter"))
    expect_true(exists("fh_run_ml_models"))
})
