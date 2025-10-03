# R/globals.R  — silence R CMD check NOTES for globals and NSE columns

## 1) Create visible bindings for the global result containers
##    (this fixes: "no visible binding for '<<-' assignment to ...")
all_result_acc        <- NULL
all_result_recall     <- NULL
all_result_FS         <- NULL
all_result_summary    <- NULL
all_result_importance <- NULL

## 2) Declare symbols used in tidy eval / ggplot etc. (NSE column names)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # columns / symbols seen in plots & data wrangling
    "pred","real","actual","Freq","recall","precision","label","score",
    "x","y","curvetypes","ctype","curve_type","type","NA_actual_",
    # any other column names that appear in your code
    "Accuracy","Composite","Composite_sd","Epoch","Loss","Mean_Composite"
  ))
}


utils::globalVariables(c(
  "Frequency", "Gene", "Label", "Median", "Method",
  "UMAP1", "UMAP2", "Value", "tape"
))

utils::globalVariables(c("collector", "prob_val"))

utils::globalVariables(c("Composite_lo", "Composite_hi"))

utils::globalVariables(c("Group", "Sample", "Combo_z"))