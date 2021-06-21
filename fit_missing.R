# code to fit xgboost models to HOSEA data (subsample)

# starts from complete_data_raw (i.e. after link_summaries.R)
# split into test,train,validation
# impute independently with function
# prepare imputed data for xgboost
# fit models, inspect AUC, PDP plots, etc.

#### source helper functions ####
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/impute_models.R')
source('R_code/hosea-project/impute_missing.R')

#### import data ####
complete_data <- readRDS('R_data/subsample/sub_complete_data_raw.rds')

#### train/test/validation split #### 
# take a test set of 1000
set.seed(1994)
n <- nrow(complete_data)
ntest <- as.integer(.1*n)
itest <- sample(1:n,ntest) 
test <- rep(FALSE,n)
test[itest] <- TRUE
train <- !test

# select a validation set (25% of training data)
set.seed(1997)
ivalid <- sample(which(train),size=as.integer(.25*sum(train)),replace=FALSE)
valid <- rep(FALSE,n)
valid[ivalid] <- TRUE
traintrain <- as.logical(train*(!valid))

#### impute each chunk ####
train_data_impute <- impute_missing_hosea(complete_data[traintrain,],seed=1995)
test_data_impute <- impute_missing_hosea(complete_data[test,],seed=1996)
valid_data_impute <- impute_missing_hosea(complete_data[valid,],seed=1997)

#### format imputed data for xgboost ####
# with NAs
dwatchlist_na <- xgb_prep(train_data_impute,
                          test_data_impute,
                          valid_data_impute,
                          dname='clean')

# imputed with sampling
dwatchlist_samp <- xgb_prep(train_data_impute,
                          test_data_impute,
                          valid_data_impute,
                          dname='impsamp')

# imputed with regression
dwatchlist_reg <- xgb_prep(train_data_impute,
                          test_data_impute,
                          valid_data_impute,
                          dname='impreg')

#### fit with NAs ####

# global parameters for the 3 models
param_xg <- list(max_depth=5,eta=.05,objective='binary:logistic',
                 eval_metric='logloss')

# fit xgboost model
xgb_fit_na <- xgb.train(param_xg,
                     dwatchlist_na$train, # training set
                     nrounds=1000,
                     dwatchlist_na, # data watchlist
                     verbose=1,print_every_n=8,
                     early_stopping_rounds=10)

# print info about important variables
nsumm <- 20
xgb_summ_na <- xgb.importance(model=xgb_fit_na)
print(xgp_summ_na[1:nimp,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
xgb_pdp('ageatindex',xgb_fit_na,train_data_impute$clean)
# ...

# evaluate AUCs (as helper in xgb_utils.R)
xgb_auc_na <- xgb_auc(xgb_fit_na,dwatchlist_na)
print(xgb_auc_na)

#### fit with random sample imputation ####

# fit xgboost model
xgb_fit_samp <- xgb.train(param_xg,
                        dwatchlist_samp$train, # training set
                        nrounds=1000,
                        dwatchlist_samp, # data watchlist
                        verbose=1,print_every_n=8,
                        early_stopping_rounds=10)

# print info about important variables
xgb_summ_samp <- xgb.importance(model=xgb_fit_samp)
print(xgp_summ_samp[1:nimp,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
xgb_pdp('ageatindex',xgb_fit_samp,train_data_impute$impsamp)
# ...

# evaluate AUCs (as helper in xgb_utils.R)
xgb_auc_samp <- xgb_auc(xgb_fit_samp,dwatchlist_samp)
print(xgb_auc_samp)

#### fit with regression imputation ####

# fit xgboost model
xgb_fit_reg <- xgb.train(param_xg,
                        dwatchlist_reg$train, # training set
                        nrounds=1000,
                        dwatchlist_reg, # data watchlist
                        verbose=1,print_every_n=8,
                        early_stopping_rounds=10)

# print info about important variables
xgb_summ_reg <- xgb.importance(model=xgb_fit_reg)
print(xgp_summ_reg[1:nimp,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
xgb_pdp('ageatindex',xgb_fit_reg,train_data_impute$impreg)
# ...

# evaluate AUCs (as helper in xgb_utils.R)
xgb_auc_reg <- xgb_auc(xgb_fit_reg,dwatchlist_reg)
print(xgb_auc_reg)



