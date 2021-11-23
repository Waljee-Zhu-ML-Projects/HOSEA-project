setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/evaluation_split.R')

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
n_controls = 2e6
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

# prepare stuff
n = 1e5
method = "resample"

# base XGBoost parameters
param_xg_base = list(
  max_depth=5,
  subsample=0.2,
  eta=.01,
  objective='binary:logistic',
  eval_metric='auc',
  scale_pos_weight = 1.
)


# grids
values = list(
  # max_depth = c(3, 4, 5, 7, 10),
  # subsample = c(0.05, 0.1, 0.2, 0.3, 0.5),
  eta = c(0.0001, 0.0003, 0.0005, 0.001)
)

# datasets
sub_train = subsample_controls(train, n)
log(sub_train)
set.seed(0)
out = train_test_split(df=sub_train, weights=c(2, 1))
train_n = out[[1]]
valid_n = out[[2]]
n0 = (train_n$CaseControl==0)%>%sum
n1 = (train_n$CaseControl==1)%>%sum
# produce training set
train_ = switch(method,
                "downsample" = subsample_controls(train_n, n1),
                "resample" = balanced_resample(train_n),
                "unweighted" = train_n,
                "weighted" = balanced_resample(train_n),
                "weighted_sqrt" = balanced_resample(train_n),
                "weighted_sq" = balanced_resample(train_n)
)
# imputation
quantiles = compute_quantiles(train_n, n_quantiles)
set.seed(0)
train_ = impute_srs(train_, quantiles)
test_ = impute_srs(test, quantiles)
valid_ = impute_srs(valid_n, quantiles)
dwatchlist = xgb_prep(train=train_,
                      test=test_,
                      valid=valid_,
                      cc=cc_test)
dwatchlist$test = NULL
dwatchlist$cc = NULL


# prepare testing sets
test_sets = evaluation_split(test, test_)

# loop
for(param in names(values)){
  for(value in values[[param]]){
    cat("=======================================================================", fill=T)
    cat(paste("===", param, "=", value), fill=T)
    cat("=======================================================================", fill=T)
    param_xg = param_xg_base
    param_xg[[param]] = value
  
    # fit
    set.seed(0)
    param_xg$scale_pos_weight = switch(
      method,
      "downsample" = 1.,
      "resample" = 1.,
      "unweighted" = 1.,
      "weighted" = n0/n1,
      "weighted_sqrt" = sqrt(n0/n1),
      "weighted_sq" = (n0/n1)^2
    )
    xgb_fit = xgb.train(param_xg,
                        dwatchlist$train,
                        nrounds=20000,
                        dwatchlist,
                        verbose=1,print_every_n=200,
                        early_stopping_rounds=1000)
    # evaluate
    calibration_metrics = lapply(test_sets, function(df) calibration(xgb_fit, df))
    threshold_metrics = lapply(test_sets, function(df) 
      classification_metrics(xgb_fit, df, seq(0.01, 0.99, 0.01)))
    # save
    out = list(
      xgb_fit=xgb_fit,
      metrics=list(calibration=calibration_metrics, classification=threshold_metrics)
    )
    filepath = paste0("R_data/results/models_tuning_large/", param, "_", value, "_patchsmall.rds")
    saveRDS(out, filepath)
  }
}

test_sets
sapply(test_sets, function(df) sum(getinfo(df, "label")))

       