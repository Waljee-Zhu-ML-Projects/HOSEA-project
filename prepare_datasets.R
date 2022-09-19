setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(magrittr)
library(dplyr)
library(HOSEA) # requires v >=0.0.0.9008
nc = 10
# imputer_path = "./R_data/results/imputers/srs.rds"
# imputer = readRDS(imputer_path)
dir_models = "./R_data/results/models/imputation/"
imputation = "mice"
to_drop = c("chol", "rbc", "hgb")

outcome = "ANY"

y0s = c(5, 5, 5, 5, 4, 3, 2, 4)
y1s = c(1, 3, 2, 4, 1, 1, 1, 2)

# no 4 did not work, increase sample size for that one

for(i in c(1,2,6,8)){
  y0 = y0s[i]; y1 = y1s[i]

  cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  cat("[HOSEA] Loading data", y0, y1, "\n")
  
  complete_data = readRDS(paste0('R_data/processed_records/', y0, '-', y1, '_merged.rds'))
  master = complete_data$master
  df = complete_data$df
  
  
  # HOSEA.fit
  set.seed(0)
  df %<>% HOSEA:::patch_outcome(master, outcome)
  cat("[HOSEA] Selected outcome ", outcome, "\n")
  features = HOSEA:::get_feature_names(df, to_drop)
  n0all = (df$casecontrol==0)%>%sum
  n1all = (df$casecontrol==1)%>%sum
  df_list = HOSEA:::stratified_split(df, c(train=2, test=1, valid=1))
  cat("[HOSEA] Split into train/valid/test\n")
  HOSEA:::summary_df_list(df_list) %>% print()
    
  # # imputation
  # set.seed(0)
  # df_list$train %<>% sample_n(10000)
  # 
  # cat("[HOSEA] Building imputation models ...\n")
  # cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  # nc = 4
  # cluster = parallel::makeCluster(nc)
  # set.seed(0)
  # imputer = switch(
  #   imputation,
  #   "mice"=HOSEA.mice.fit(df_list$train, features, n_rounds=5, cluster=cluster),
  #   "srs"=HOSEA.srs.fit(df_list$train, features, n_quantiles=10000)
  # )
  # parallel::stopCluster(cluster)
  # cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  # cat("[HOSEA] Saving Imputer\n")
  imputer_path = paste0("./R_data/results/imputers/mice_", y0, "-", y1, ".rds")
  # saveRDS(imputer, imputer_path)
  # cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
# End imputation
# df_list %<>% lapply(function(df) df %>% sample_n(10000))
  
  # Training
  nrounds=2000
  print_every_n=10
  early_stopping_rounds=100
  imputer = readRDS(imputer_path)
  
  missingness = HOSEA:::summarize_missingness(df_list$train)
  set.seed(0)
  df_list$train %<>% HOSEA:::balanced_resample()
  ids = lapply(df_list, function(df) df$id %>% unique())
  cat("[HOSEA] Resample cases in training set\n")
  HOSEA:::summary_df_list(df_list) %>% print()
  
  set.seed(0)
  cat("[HOSEA] Imputation starting (valid)\n")
  cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  cluster = parallel::makeCluster(nc)
  df_list$valid = impute(imputer, df_list$valid, 1, cluster=cluster)
  filepath = paste0("./R_data/imputed_records/", y0, '-', y1, "valid_", imputation, "_", tolower(outcome), ".rds")
  print(filepath)
  saveRDS(df_list$valid, filepath)
  parallel::stopCluster(cluster)
  
  set.seed(0)
  cat("[HOSEA] Imputation starting (test)\n")
  cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  cluster = parallel::makeCluster(nc) 
  df_list$test = impute(imputer, df_list$test, 1, cluster=cluster)
  filepath = paste0("./R_data/imputed_records/", y0, '-', y1, "test_", imputation, "_", tolower(outcome), ".rds")
  print(filepath)
  saveRDS(df_list$test, filepath)
  parallel::stopCluster(cluster)
  
  set.seed(0)
  cat("[HOSEA] Imputation starting (train)\n")
  cat("[HOSEA]", timestamp(prefix="", suffix="", quiet=T), "\n")
  cluster = parallel::makeCluster(nc)
  df_list$train = impute(imputer, df_list$train, 1, cluster=cluster)
  filepath = paste0("./R_data/imputed_records/", y0, '-', y1, "train_", imputation, "_", tolower(outcome), ".rds")
  print(filepath)
  saveRDS(df_list$train, filepath)
  parallel::stopCluster(cluster)
  
  
}
