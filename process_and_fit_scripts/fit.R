setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(magrittr)
library(dplyr)
library(HOSEA) # requires v >=0.0.0.9008

dir_models = "./R_data/results/models/imputation3/"
imputation = "mice"
to_drop = c("chol", "rbc", "hgb")

# Training
nrounds=2000
print_every_n=10
early_stopping_rounds=100
  

complete_data = readRDS('R_data/processed_records/5-1_merged.rds')
master = complete_data$master
# df = complete_data$df

# HOSEA.fit
outcome = "ANY"


cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
cat("[HOSEA] Loading data\n")
df_list = list()

filepath = paste0("./R_data/imputed_records/5-1valid_", imputation, "_", tolower(outcome), ".rds")
print(filepath)
df_list$valid = readRDS(filepath)

cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")

filepath = paste0("./R_data/imputed_records/5-1test_", imputation, "_", tolower(outcome), ".rds")
print(filepath)
df_list$test = readRDS(filepath)

cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")

filepath = paste0("./R_data/imputed_records/5-1train_", imputation, "_", tolower(outcome), ".rds")
print(filepath)
df_list$train = readRDS(filepath)

cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
features = HOSEA:::get_feature_names(df_list$train, to_drop)

for(outcome in c("ANY", "EAC", "EGJAC")){
  
  
  # prepare outcome
  cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  cat("[HOSEA] Selecting outcome ", outcome, "\n")
  df_list1 = df_list %>% lapply(function(df) patch_outcome(df, master, outcome, drop=T))
  
  n0train = (df_list1$train$casecontrol==0)%>%sum
  n1train = (df_list1$train$casecontrol==1)%>%sum
  
  n0test = (df_list1$test$casecontrol==0)%>%sum
  n1test = (df_list1$test$casecontrol==1)%>%sum
  
  cat("[HOSEA] Summary of train/valid/test\n")
  HOSEA:::summary_df_list(df_list1) %>% print()
  
  
  
  watchlist = HOSEA:::prepare_watchlist(df_list1, features)
  params_xgb = HOSEA:::xgboost_options(
    eta=0.05
  )
  params_xgb$scale_pos_weight = (n1test/n1train) * (n0train/n0test)
  cat("[HOSEA] Starting XGBoost fitting\n")
  cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  set.seed(0)
  xgb_fit = xgboost::xgb.train(
    params=params_xgb,
    data=watchlist$train,
    watchlist=watchlist,
    nrounds=nrounds,
    verbose=1, 
    print_every_n=print_every_n,
    early_stopping_rounds=early_stopping_rounds
  )
  cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  out = list(
    xgb_fit=xgb_fit,
    imputation=imputation
  )
  class(out) = "HOSEA.fit"
  
  filepath = paste0(dir_models, "xgb_", imputation, "_", tolower(outcome), ".rds")
  print(filepath)
  saveRDS(out, filepath)
  
  filepath = paste0(dir_models, "xgb_", imputation, "_", tolower(outcome), ".model")
  print(filepath)
  xgboost::xgb.save(xgb_fit, filepath)
  
}

