source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')

# import data
complete_data = readRDS('R_data/complete_data_raw.rds')
complete_data$n_visits = NULL 
cc_test <- readRDS('R_data/cc_complete_data.rds')

# replciability
set.seed(0)
n_quantiles = 10000

# logging
log = function(df) cat(paste("Full data set: ", nrow(df), "observations,", 
                             df$CaseControl%>%sum, "cases", 
                             (df$CaseControl==0)%>%sum, "controls"), fill=T)

# subsampling of controls
n_controls = 1e6
sub_complete_data = subsample_controls(complete_data, n_controls)
log(sub_complete_data)
rm(complete_data)
gc()

# train-test split
out = train_test_split(df=sub_complete_data, weights=c(3, 1))
train = out[[1]]
test = out[[2]]

log(train)
log(test)

# loop over n from here
ns = c(1e4, 2e4) #, 5e4, 1e5, 2e5, 5e5)

ns = c(1e6)
# loop
for(n in ns){
  # prepare outpuut
  aucs = data.frame(matrix(0., 1, 6))
  colnames(aucs) = c("method", "n", "train", "valid", "test", "cc")
  sub_train = subsample_controls(train, n)
  log(sub_train)
  out = train_test_split(df=sub_train, weights=c(2, 1))
  train_n = out[[1]]
  valid_n = out[[2]]
  log(train_n)
  log(valid_n)
  n0 = (train_n$CaseControl==0)%>%sum
  n1 = (train_n$CaseControl==1)%>%sum
  
  # downsample
  train_downsample = subsample_controls(train_n, n1)
  log(train_downsample)
  
  quantiles_downsample = compute_quantiles(train_downsample, n_quantiles)
  set.seed(0)
  train_downsample = impute_srs(train_downsample, quantiles_downsample)
  test_downsample = impute_srs(valid_n, quantiles_downsample)
  valid_downsample = impute_srs(test, quantiles_downsample)
  
  # out = impute_missing_hosea(
  #   train=train_downsample, 
  #   valid=valid_n, 
  #   test=test, seed=0)$dfs
  # train_downsample = out$train
  # test_downsample = out$test
  # valid_downsample = out$valid
  
  dwatchlist_downsample = xgb_prep(train_downsample,
                                   test_downsample,
                                   valid_downsample,
                                   # dname='impsamp',
                                 cc=cc_test)
  
  # resampling
  train_resample = balanced_resample(train_n)
  log(train_resample)
  
  quantiles_resample = compute_quantiles(train_resample, n_quantiles)
  set.seed(0)
  train_resample = impute_srs(train_resample, quantiles_resample)
  test_resample = impute_srs(valid_n, quantiles_resample)
  valid_resample = impute_srs(test, quantiles_resample)
  
  # out = impute_missing_hosea(
  #   train=train_resample, 
  #   valid=valid_n, 
  #   test=test, seed=0)$dfs
  # train_resample = out$train
  # test_resample = out$test
  # valid_resample = out$valid
  
  dwatchlist_resample = xgb_prep(train_resample,
                                 test_resample,
                                 valid_resample,
                                  # dname='impsamp',
                                  cc=cc_test)
  
  # reweighing cases
  
  
  quantiles_weight = compute_quantiles(train_n, n_quantiles)
  set.seed(0)
  train_weight = impute_srs(train_resample, quantiles_weight)
  test_weight = impute_srs(valid_n, quantiles_weight)
  valid_weight = impute_srs(test, quantiles_weight)
  
  # out = impute_missing_hosea(
  #   train=train_n, 
  #   valid=valid_n, 
  #   test=test, seed=0)$dfs
  # train_weight = out$train
  # test_weight = out$test
  # valid_weight = out$valid
  

  dwatchlist = xgb_prep(train_weight,
                               test_weight,
                               valid_weight,
                               # dname='impsamp',
                               cc=cc_test)
  
  
  # XGBoost parameters
  param_xg = list(
    max_depth=5,
    subsample=min(0.5, sqrt(10000/n)),
    eta=.05,
    objective='binary:logistic',
    eval_metric='auc'
  )
  
  # fit
  set.seed(1)
  
  param_xg$scale_pos_weight = 1.
  xgb_fit = xgb.train(param_xg,
                      dwatchlist$train,
                      nrounds=5000,
                      dwatchlist,
                      verbose=1,print_every_n=10,
                      early_stopping_rounds=50)
  
  param_xg$scale_pos_weight = 1.
  xgb_fit_downsample = xgb.train(param_xg,
                               dwatchlist_downsample$train,
                               nrounds=5000,
                               dwatchlist_downsample,
                               verbose=1,print_every_n=10,
                               early_stopping_rounds=50)
  
  param_xg$scale_pos_weight = 1.
  xgb_fit_resample = xgb.train(param_xg,
                               dwatchlist_resample$train,
                               nrounds=5000,
                               dwatchlist_resample,
                               verbose=1,print_every_n=10,
                               early_stopping_rounds=50)
  
  param_xg$scale_pos_weight = n0/n1
  xgb_fit_weight = xgb.train(param_xg,
                             dwatchlist$train,
                             nrounds=5000,
                             dwatchlist,
                             verbose=1,print_every_n=10,
                             early_stopping_rounds=50)
  
  param_xg$scale_pos_weight = sqrt(n0/n1)
  xgb_fit_weight_sqrt = xgb.train(param_xg,
                                  dwatchlist$train,
                                  nrounds=5000,
                                  dwatchlist,
                                  verbose=1,print_every_n=10,
                                  early_stopping_rounds=50)
  
  param_xg$scale_pos_weight = (n0/n1)^2
  xgb_fit_weight_sq = xgb.train(param_xg,
                                  dwatchlist$train,
                                  nrounds=5000,
                                  dwatchlist,
                                  verbose=1,print_every_n=10,
                                  early_stopping_rounds=50)
  
  aucs = rbind(aucs, c("unweighted", n, best_auc(xgb_fit)))
  aucs = rbind(aucs, c("downsample", n, best_auc(xgb_fit_downsample)))
  aucs = rbind(aucs, c("resample", n, best_auc(xgb_fit_resample)))
  aucs = rbind(aucs, c("weighted", n, best_auc(xgb_fit_weight)))
  aucs = rbind(aucs, c("weighted_sqrt", n, best_auc(xgb_fit_weight_sqrt)))
  aucs = rbind(aucs, c("weighted_sq", n, best_auc(xgb_fit_weight_sq)))
  
  write.csv(aucs, paste0("R_data/results/best_aucs_large_", n, ".csv"))
}

aucs = read.csv("R_data/results/best_aucs_large.csv")
aucs$X = NULL
aucs = aucs[-1, ]
aucs = aucs[order(aucs$n), ]
aucs = aucs[order(aucs$method), ]
xtable::xtable(aucs, digits=3, include.rownames=F)

library(ggplot2)
ggplot(aucs, aes(n, cc, group=method, colour=method)) + 
  geom_line()
