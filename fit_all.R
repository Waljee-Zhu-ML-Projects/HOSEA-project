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
complete_data = readRDS('R_data/processed_records/5-1.rds')
master = complete_data$master
complete_data = complete_data$df


# replciability
set.seed(0)
n_quantiles = 10000

# logging
log = function(df) cat(paste("Full data set: ", nrow(df), "observations,", 
                             df$CaseControl%>%sum, "cases", 
                             (df$CaseControl==0)%>%sum, "controls"), fill=T)

n0all = (complete_data$CaseControl==0)%>%sum
n1all = (complete_data$CaseControl==1)%>%sum
# train-test split
# n_controls = 2e6
# complete_data = subsample_controls(complete_data, n_controls)
set.seed(0)
out = train_test_split(df=complete_data, weights=c(3, 1))
train = out[[1]]
test = out[[2]]
rm(complete_data)
gc()

log(train)
log(test)

# prepare stuff
methods = c(
  "resample"
  # "weighted",
  # "weighted_sqrt"
)
param_xg = list(
  max_depth = 5,
  subsample = 0.1,
  eta = .05,
  objective = 'binary:logistic',
  eval_metric = 'auc',
  nthread=-1
)

# loop
n = "all"
set.seed(0)
out = train_test_split(df=train, weights=c(2, 1))
rm(train);gc()
train_n = out[[1]]
valid_n = out[[2]]
rm(out);gc()
log(train_n)
log(valid_n)
n0 = (train_n$CaseControl==0)%>%sum
n1 = (train_n$CaseControl==1)%>%sum
method = "resample"
timestamp()
# produce training set
train_ = switch(method,
               "downsample" = subsample_controls(train_n, n1),
               "resample" = balanced_resample(train_n),
               "unweighted" = train_n
)
log(train_)
# imputation
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
dwatchlist$test = NULL
gc()
# fit
set.seed(0)
param_xg$scale_pos_weight = switch(
  method,
  "downsample" = 1.,
  "resample" = (n1all/n0all)*(n0/n0),
  "unweighted" = (n1all/n0all)*(n0/n1),
  "weighted" = n0/n1,
  "weighted_sqrt" = sqrt(n0/n1),
  "weighted_sq" = (n0/n1)^2
)
timestamp()
xgb_fit = xgb.train(param_xg,
                    dwatchlist$train,
                    nrounds=20000,
                    dwatchlist,
                    verbose=1,print_every_n=10,
                    early_stopping_rounds=100)
timestamp()
# save
out = list(
  xgb_fit=xgb_fit,
  quantiles=quantiles,
  test_ids=test_$ID
)
filepath = paste0("R_data/results/models/XGB_nALL_typeANY.rds")
saveRDS(out, filepath)