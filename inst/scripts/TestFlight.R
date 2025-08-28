# This is an example script to demonstrate the main functions of FeatureHunter package
# welcome here and thanks for your interest in FeatureHunter
# on behalf of PNCS headquarters, we sincerely hope you will find FeatureHunter useful and user-friendly
# If you are in need of any help, we will be glad to assist you 
#Contact me via 
#Rednote Number 3683833756
#Instagram https://www.instagram.com/zackerweiss?igsh=MXFpYTN6NGU5OG56aQ%3D%3D&utm_source=qr
#Of course you can contact me on Github. I may not reply such a quick way but I will definitely reply you
#Enjoy FeatureHunter and have a nice day!

library(FeatureHunter)

#Input at least 4 datasets, 1 as training set called DatasetA, 1 as validation set DatasetB, 2 as test sets DatasetC and DatasetD, save as DatasetA.txt, DatasetB.txt, DatasetC.txt, DatasetD.txt
#Input at least 4 sample information tables, similar to my example document, named SampleA.txt, SampleB.txt, SampleC.txt, SampleD.txt,
#Input a common gene list file, named Common_genes.txt

#If more than 4 datasets are provided, the 3rd and subsequent datasets will be treated as test sets
dir_sample <- fh_example_path("raw_data_pncs") #this is the path to your data folder, modify it to your own data folder
# dir_sample <- "your_data_folder_path"
dat <- fh_load_data(
  dir            = dir_sample,
  train_exp_name = "DatasetA.txt",
  hub_file       = "Common_genes.txt",
  positive_label = "Disease",
  assign_to_global = TRUE,       # write the objects to global env
)

# 现在全局里就有：
# allfiles, labels_files, labels_list, exp_files, exp_list,
# com_genes, all_labels, x_train, y_train, x_tests


out <- fh_preprocess_data(
  exp_list       = exp_list,
  labels_list    = labels_list,
  exp_files      = exp_files,
  train_exp_name = "DatasetA.txt",
  ref_index      = 1,           #  use DatasetA( training set) as reference
  include_ref    = TRUE,        # training set also used for standardization
  align          = "exact"      
)

print(out$pos_rate_summary)  # postive rate summary
out$log_report               # log decision report
out$scale_report             # standardization report

exp_list   <- out$exp_list   
train_exp <- out$train_exp 
test_exp   <- out$test_exp

lasso <- fh_lasso_prescreen(
  x_train = x_train,
  y_train = y_train,
  x_tests = test_exp,
  alpha = 1,
  nfolds = 10,
  nlambda = 100,
  intercept = FALSE,
  impute_median = TRUE,
  selection = "conservative",   # choose from c("lenient","balanced","conservative","custom_fraction") SEE DETAILS BELOW
  fraction = 0.5,          # only used when selection = "custom_fraction"
  min_keep = 0,
  max_keep = Inf,
  plot = "both"       
)

# DO NOT RUN the following code chunk, it is only for your reference to understand how the lasso prescreening works
# lenient       = lambda_min,
# conservative  = lambda_1se,
# balanced      = exp(0.5 * (log(lambda_min) + log(lambda_1se))),  
# custom_fraction = {
#   if (!is.finite(fraction) || fraction < 0 || fraction > 1)
#     stop("fraction must be in [0,1] when selection='custom_fraction'.")
#   exp((1 - fraction) * log(lambda_min) + fraction * log(lambda_1se))
# }



train_expd_lasso <- lasso$x_train
test_exp_lasso   <- lasso$x_tests  
lassogene        <- lasso$genes    



all_result_summary <- list()
all_result_acc     <- list()
all_result_recall  <- list()
all_result_FS      <- list()


collector <- new.env(parent = emptyenv())
mlp_result_pro <- fh_mlp(
  x_train,
  y_train,
  test_exp = test_exp,
  epochselecti = 50,
  batchsize_grid =  c(32L, 48L, 64L, 96L, 128L),
  learningrate_grid = c(1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2),
  dropoutrate_grid  = c(0.10, 0.25, 0.40, 0.50),
  labels_list = labels_list,
  collector = collector ,
  hidden_units = c(32L, 16L, 8L),
  use_batchnorm      = TRUE,
  l2                 = 1e-4, 
  gaussian_noise_sd  = 0.01,
  min_lr             = 1e-5,
  plateau_factor     =  0.5,
  plateau_patience   = 4,
  early_patience     = 10,
  seed               = 42
)

snap_path <- file.path(getwd(), "collector_after_mlp.RData")
save(collector, file = snap_path)

# ########### machine learnings one by one(use when necessary) #########
# fh_logistic(
#   train_exp = train_exp,
#   train_labels = y_train,
#   test_exp = test_exp,  
#   labels_list = labels_list,
#   all_labels = all_labels,
#   auto_th_method = "auto" , 
#   collector = collector
# )
# 
# fh_lda(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   auto_th_method = "auto" , 
#   collector = collector
# )
# 
# fh_lasso_lda(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   lassogene    = lassogene,   
#   fold         = 10,      
#   auto_th_method = "auto" ,
#   collector = collector
# )
# 
# fh_qda(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   auto_th_method = "auto" ,
#   collector = collector
# )
# 
# fh_lasso_qda(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   lassogene    = lassogene,
#   fold         = 10,        
#   auto_th_method = "auto" ,
#   collector = collector
# )
# 
# fh_knn(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   knumber      = c(1,2,3,4,5),
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_lasso_knn(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   lassogene    = lassogene,
#   fold         = 10,
#   knumber      = c(1,2,3,4,5),
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# fh_tree(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# fh_lasso_tree(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   lassogene    = lassogene,
#   fold         = 10,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_rf(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   com_genes    = com_genes,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_lasso_rf(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   lassogene    = lassogene,
#   fold         = 10,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_xgboost(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   nround       = 100,
#   max_depth    = 6,
#   eta          = 0.5,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# fh_lasso_xgboost_default(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,
#   labels_list  = labels_list,
#   lassogene    = lassogene,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_lasso_xgboost_best(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,
#   labels_list  = labels_list,
#   fold         = 10,
#   lassogene    = lassogene,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_ridge_glmnet(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   fold         = 10,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_lasso(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   fold         = 10,
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# 
# fh_elastic_net(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   fold         = 10,
#   alpha_all    = seq(0, 1, 0.1),
#   auto_th_method = "auto",
#   collector = collector  
# )
# 
# fh_svm(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   kernel_all   = c("linear","polynomial","radial"),
#   auto_th_method = "auto",
#   collector = collector  
# )
# 
# 
# fh_svm_best(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   fold         = 10,
#   lassogene    = lassogene,
#   collector = collector  
# )
# 
# 
# fh_lasso_svm_best(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   fold         = 10,        
#   lassogene    = lassogene,
#   collector = collector  
# )
# 
# fh_gbm_default(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   auto_th_method = "auto",
#   collector = collector  
# )
# 
# fh_gbm_cv_best(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,
#   fold         = 10,
#   auto_th_method = "auto",
#   collector = collector  
# )
# 
# 
# fh_gbm_lasso_default(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  ...
#   labels_list  = labels_list,     
#   lassogene    = lassogene,    
#   fold         = 10,         
#   auto_th_method = "auto",
#   collector = collector
# )
# 
# fh_gbm_lasso_best(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  ...
#   labels_list  = labels_list,   
#   lassogene    = lassogene,     
#   fold         = 10,  
#   collector = collector
# )
# 
# 
# fh_stepwise_lr(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,  
#   labels_list  = labels_list,  
#   collector = collector
# )
# 
# 
# 
# fh_naive_bayes(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,
#   labels_list  = labels_list,  
#   collector = collector
# )
# 
# 
# fh_naive_bayes_lasso(
#   train_exp    = train_exp,
#   train_labels = y_train,
#   test_exp     = test_exp,
#   labels_list  = labels_list,
#   lassogene    = lassogene,
#   fold         = 10,  
#   collector = collector
# )
# 

############# machine learnings all together#####
fh_run_ml_models(
  train_exp    = train_exp,
  train_labels = y_train,
  test_exp     = test_exp,
  labels_list  = labels_list,
  lassogene    = lassogene,
  com_genes    = com_genes,
  all_labels   = all_labels,
  fold         = 10,
  knumber      = c(1,3,5),
  alpha_all    = seq(0.1,0.9,0.1),
  kernel_all   = c("linear","polynomial","radial"),
  nround       = 100,
  max_depth    = 6,
  eta          = 0.5,
  auto_th_method = "auto",
  cores        = 4,                  # auto select (parallel::detectCores()-1)
  which_models = "all"   
  # use which_models = c("fh_rf","fh_gbm_cv_best") when you don't prefer to run all of machine learning methods
)


all_result_summary <- as.list(collector$all_result_summary)
all_result_acc     <- as.list(collector$all_result_acc)
all_result_recall  <- as.list(collector$all_result_recall)
all_result_FS      <- as.list(collector$all_result_FS)


## Visualization of models across datasets
nshow <- 40 #heatmap lines you want to show

fh_all_heatmap(
  all_result_acc    = all_result_acc,
  all_result_recall = all_result_recall,
  all_result_FS     = all_result_FS,
  exp_files          = exp_files,  
  heatmapcolor      = c("#8ecae6","#219ebc","#ffc8dd","#ffafcc","#023047","#ffb703","#fb8500",
                        "#bde0fe","#a2d2ff","#bde0fe"),
  nshow             = nshow        
)

# ###optional
# ###train a single MLP model with multiple random seeds (optional)
# res_single <- fh_mlp_single(
#   train_exp    = x_train,
#   train_labels = y_train,
#   test_exp     = test_exp,
#   epochselecti = 50,
#   batchsize    = 64L,
#   learningrate = 1e-2,         
#   dropoutrate  = 0.40,
#   cutoff_grid  = sort(unique(c(
#     seq(0.40, 0.70, by = 0.05),
#     seq(0.70, 0.90, by = 0.01),
#     seq(0.80, 0.84, by = 0.001)
#   ))),
#   imbalance_thresh = 0.35,
#   auto_th_method   = "youden",
#   early_patience   = 10,
#   seeds        = c(42L,43L,44L,45L,46L,47L,48L,49L,50L,51L),
#   labels_list  = labels_list,
#   collector    = collector,
#   standardize  = FALSE,
#   hidden_units = c(32L,16L,8L),
#   activation   = "relu",
#   use_batchnorm = TRUE,
#   l2           = 1e-3,        
#   gaussian_noise_sd = 0.05,
#   min_lr       = 1e-6,
#   plateau_factor   = 0.5,
#   plateau_patience = 4
# )


fh_roc(
  all_result_summary      = all_result_summary,
  top_n       = 6,
  skip_first  = TRUE,
  namesroc    = c("DatasetA(Train)", "DatasetB(Val)", "DatasetC(Test)", "DatasetD(Test)"),
  out_dir     = getwd(),      #optional,default is current working directory
  out_modelnames_csv = "modelnames.csv",
  out_dir_roc    = "roc",
  out_dir_result = "result"
)


####parameters for strict filtering
# want to be more "strict": increase strict_min_freq to 0.70 or strict_min_effect to 0.60.
# want to be more "loose": decrease strict_min_freq to 0.60, keep strict_ci_gate = TRUE.
# if you want to apply strict filtering to all plots, set apply_selection_to_plots = TRUE.

# if you choose MLP model, it takes several minutes to wait for hunter results, other models are faster
res <- fh_hunter(
  train_exp = train_exp,
  train_labels = y_train,
  nshow = 40,
  namesS = c("Accuracy","Recall","F-score"),
  score_index = 3,# 1 refers to Accuracy, 2 to Recall, 3 to F-score ,actually the index of the metric you care most in the top nshow models you exported
  pick_index = 2,# the rank of the model you like in the dataset you specified above 
  #Ooops,the first model's threshold is too low(0.01) and the confusion curve is not good in A, so I choose the second model, the effect is actually similar
  seed = 42,
  num_runs = 10,
  num_coregene = 29,     #genes to show in the composite barplot
  n_likes = 6,           #genes for UMAP and Regression
  n_interest = 4,        #genes for stability heatmap
  
  ##parameters below are optional and have default values adjust when necessary
  perm_metric    = "f1",
  perm_nrep      = 10L,
  shap_nsim      = 30L,
  shap_subsample = 200L,
  shuffle            = TRUE,
  standardize        = TRUE,
  # parameters below are for MLP model, need to be consistent with previous MLP training parameters, otherwise the results will differ
  hidden_units       = c(32L,16L,8L),
  activation         = "relu",
  use_batchnorm      = TRUE,
  l2                 = 1e-4,
  gaussian_noise_sd  = 0.01,
  min_lr             = 1e-5,
  plateau_factor     = 0.5,
  plateau_patience   = 4L,
  early_patience     = 10L,
  imbalance_thresh   = 0.35,
  auto_th_method     = "auto",
  # knobs for strict filtering,DO NOT change the default values unless you know what you are doing
  method_weights     = c(model = 1, perm = 1, shap = 1),
  ci_level           = 0.99,
  strict_min_freq    = 0.75,
  strict_min_effect  = 0.55,
  strict_ci_gate     = TRUE,
  strict_knee        = TRUE,
  strict_nmax        = 10,
  apply_selection_to_plots = FALSE  #whether to activate strict filtering mentioned above
)






# top by F-score, pick the 1st ranked row (any model type)
out <- fh_plot_model(
  all_result_summary = all_result_summary,
  all_result_FS      = all_result_FS,
  nshow        = 40,
  namesS       = c("Accuracy","Recall","F-score"),
  score_index  = 3,
  pick_index   = 2,
  cm_palette_low  = "#20ACBD",
  cm_palette_mid  = "#FCECCF",
  cm_palette_high = "#F69896",
  base_size = 11
)

# # pick the 7th ranked row from an explicit CSV
# out2 <- fh_plot_model(
#   all_result_summary = all_result_summary,
#   all_result_FS      = all_result_FS,
#   top_models_csv = "heatmap/40_top_F-score_models.csv",
#   pick_index   = 7
# )
# 
# 





fh_plot_PRROC(
  n_model = 10,               # displaying top 10 models
  score_index = 3,   # use accuracy ranking
  out_pdf = "PR-ROC.pdf",
  nshow = nshow
)


#save the important objects to an RData file
save.image(file = "FeatureHunter_Workspace.RData")















