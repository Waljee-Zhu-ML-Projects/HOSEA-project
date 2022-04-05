setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/evaluation_split.R')
library(HOSEA)
library(dplyr)
library(magrittr)
# import data
complete_data = readRDS('R_data/processed_records/5-1.rds')
master = complete_data$master
complete_data = complete_data$df


# nas = complete_data %>% select(chol_mean, hdl_mean, ldl_mean) %>% is.na()
# patterns = apply(nas, 1, function(x) paste0(as.integer(x), collapse=""))
# table(patterns)


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
n_controls = 1e6
complete_data = subsample_controls(complete_data, n_controls)
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
  subsample = 0.5,
  eta = .2,
  objective = 'binary:logistic',
  eval_metric = 'auc',
  nthread=-1
)

# loop
n = "1M"
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

clinical = c("colonoscopy_n", "colonoscopy_maxdiff",
             "labs_fobt_n", "labs_fobt_maxdiff")
hgb = c("hgb_mean", "hgb_min", "hgb_max", "hgb_mindiff", "hgb_maxdiff", "hgb_tv")
hct = c("hct_mean", "hct_min", "hct_max", "hct_mindiff", "hct_maxdiff", "hct_tv")
rbc = c("rbc_mean", "rbc_min", "rbc_max", "rbc_mindiff", "rbc_maxdiff", "rbc_tv")
chol = c("chol_mean", "chol_min", "chol_max", "chol_mindiff", "chol_maxdiff", "chol_tv")
ldl = c("ldl_mean", "ldl_min", "ldl_max", "ldl_mindiff", "ldl_maxdiff", "ldl_tv")
hdl = c("hdl_mean", "hdl_min", "hdl_max", "hdl_mindiff", "hdl_maxdiff", "hdl_tv")


subsets = list(
  "all" = c(),
  "proposed" = c(clinical, hgb, rbc, chol),
  "no_clinical" = clinical,
  "no_rbc_hgb" = c(hgb, rbc),
  "no_rbc_hct" = c(hct, rbc),
  "no_hct_hgb" = c(hgb, hct),
  "no_chol" = c(chol),
  "no_ldl" = c(ldl),
  "no_hdl" = c(hdl),
  "no_hdl_ldl" = c(hdl, ldl)
)

for(i in seq(length(subsets))){
  mname = names(subsets)[i]
  subset = subsets[[i]]
  dwatchlist = xgb_prep(train=train_ %>% select(-subset),
                        test=test_ %>% select(-subset),
                        valid=valid_ %>% select(-subset))
  dwatchlist$test = NULL
  
  # fit
  set.seed(0)
  n0all = n_controls
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
  cat(xgb_fit$nfeatures)
  # save
  out = list(
    xgb_fit=xgb_fit,
    quantiles=quantiles,
    test_ids=test_$ID
  )
  filepath = paste0("R_data/results/models/XGB_n", n, "_typeANY_", mname, ".rds")
  saveRDS(out, filepath)
}


