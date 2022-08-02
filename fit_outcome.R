setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/evaluation_split.R')
library(HOSEA)

# import data
complete_data = readRDS('R_data/processed_records/5-1_merged.rds')
master = complete_data$master
dff = complete_data$df

# drop columns
dff %<>% select(-c(starts_with("chol"), starts_with("rbc"), starts_with("hgb")))

# get outcomes
outcomes = master %>% select(id, casecontrol, cancertype)
outcomes %<>% mutate(
  ANY=casecontrol,
  EAC=as.integer(cancertype=="EAC"),
  EGJAC=as.integer(cancertype=="EGJAC")
)

# parameters
set.seed(0)
n_quantiles = 10000
method = "resample"
# n = 1e6
nname = "all"

# logging
log = function(df) cat(paste("Full data set: ", nrow(df), "observations,", 
                             df$casecontrol%>%sum, "cases", 
                             (df$casecontrol==0)%>%sum, "controls"), fill=T)
# old settings
param_xg = list(
  max_depth = 5,
  subsample = 0.1,
  eta = 2,
  objective = 'binary:logistic',
  eval_metric = 'auc',
  nthread=-1
)


# prepare stuff
# param_xg = list(
#   max_depth = 4,
#   subsample = 0.2,
#   eta = .02,
#   objective = 'binary:logistic',
#   eval_metric = 'auc',
#   nthread=-1
# )

outcome_list = c("ANY", "EAC", "EGJAC")
outcome_list = c("EGJAC")
for(outcome in outcome_list){
  # select outcome
  outcomes_ = outcomes%>%select(id, !!outcome)
  df = dff %>% 
    left_join(outcomes_, by="id") %>%
    select(-casecontrol) %>% 
    rename(casecontrol=!!outcome)
    
  n0all = (df$casecontrol==0)%>%sum
  n1all = (df$casecontrol==1)%>%sum
  # train-test split
  # df = subsample_controls(df, n)
  set.seed(0)
  out = train_test_split(df=df, weights=c(3, 1))
  train = out[[1]]
  test = out[[2]]
  rm(df);gc()
  
  log(train)
  log(test)
  
  # validation split
  set.seed(0)
  out = train_test_split(df=train, weights=c(2, 1))
  rm(train);gc()
  train_n = out[[1]]
  valid_n = out[[2]]
  rm(out);gc()
  log(train_n)
  log(valid_n)
  n0 = (train_n$casecontrol==0)%>%sum
  n1 = (train_n$casecontrol==1)%>%sum
  timestamp()
  # produce training set
  train_ = balanced_resample(train_n)
  log(train_)
  # imputation
  missing_prop = train_n %>% mutate(across(everything(), is.na)) %>% colMeans()
  missing_rate = train_n %>% mutate(across(everything(), is.na)) %>% rowSums() %>% table()
  quantiles = compute_quantiles(train_n, n_quantiles)
  rm(train_n);gc()
  set.seed(0)
  train_ = impute_srs(train_, quantiles)
  test_ = impute_srs(test, quantiles)
  valid_ = impute_srs(valid_n, quantiles)
  log(train_);log(valid_);log(test_)
  dwatchlist = xgb_prep(train=train_,
                        test=test_,
                        valid=valid_)
  # dwatchlist$test = NULL
  gc()
  # fit
  param_xg$scale_pos_weight = (n1all/n0all)*(n0/n0)
  set.seed(0)
  timestamp()
  xgb_fit = xgb.train(param_xg,
                      dwatchlist$train,
                      nrounds=2000,
                      dwatchlist,
                      verbose=1,print_every_n=10,
                      early_stopping_rounds=100)
  timestamp()
  # save
  out = list(
    xgb_fit=xgb_fit,
    quantiles=quantiles,
    test_ids=test_$id,
    missing_prop=missing_prop,
    missing_rate=missing_rate
  )
  filepath = paste0("R_data/results/models/XGB_", nname, "_", outcome, ".rds")
  print(filepath)
  saveRDS(out, filepath)
  
  filepath = paste0("R_data/results/models/xgb_", tolower(outcome), ".model")
  print(filepath)
  xgboost::xgb.save(HOSEA::XGB_EGJAC$xgb_fit, filepath)

}
