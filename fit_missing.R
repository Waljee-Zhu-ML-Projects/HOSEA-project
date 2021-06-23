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
valid_data_impute <- impute_missing_hosea(complete_data[valid,],seed=1998)

# save all imputed data
save(train_data_impute,test_data_impute,valid_data_impute,
     file='R_data/subsample/sub_complete_data_impute.RData')

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
param_xg <- list(max_depth=4,eta=.05,objective='binary:logistic',
                 eval_metric='logloss')

# fit xgboost model
xgb_fit_na <- xgb.train(param_xg,
                     dwatchlist_na$train, # training set
                     nrounds=1000,
                     dwatchlist_na, # data watchlist
                     verbose=1,print_every_n=8,
                     early_stopping_rounds=10)
# stops at 146 trees

# print info about important variables
nsumm <- 20
xgb_summ_na <- xgb.importance(model=xgb_fit_na)
print(xgb_summ_na[1:nsumm,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
# age
xgb_pdp('ageatindex',xgb_fit_na,train_data_impute$clean)
# smoking
xgb_pdp('smoke_current',xgb_fit_na,train_data_impute$clean)
# heart failure
xgb_pdp('CHF',xgb_fit_na,train_data_impute$clean)
# blood labs: HCT, MCH, RBC
xgb_pdp('hct_mean',xgb_fit_na,train_data_impute$clean)
xgb_pdp('mch_mean',xgb_fit_na,train_data_impute$clean)
xgb_pdp('rbc_mean',xgb_fit_na,train_data_impute$clean)

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
# stops at 109 trees

# print info about important variables
xgb_summ_samp <- xgb.importance(model=xgb_fit_samp)
print(xgb_summ_samp[1:nsumm,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
# age
xgb_pdp('ageatindex',xgb_fit_samp,train_data_impute$impsamp)
# weight
xgb_pdp('weight',xgb_fit_samp,train_data_impute$impsamp)
# labs: WBC, BUN, Na (min)
xgb_pdp('wbc_mean',xgb_fit_samp,train_data_impute$impsamp)
xgb_pdp('bun_mean',xgb_fit_samp,train_data_impute$impsamp)
xgb_pdp('na_min',xgb_fit_samp,train_data_impute$impsamp)

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
# stops at 244 trees

# print info about important variables
xgb_summ_reg <- xgb.importance(model=xgb_fit_reg)
print(xgb_summ_reg[1:nsumm,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
# age
xgb_pdp('ageatindex',xgb_fit_reg,train_data_impute$impreg)
xgb_pdp('hct_min_diff',xgb_fit_reg,train_data_impute$impreg) # probably an artifact of imputation?
xgb_pdp('lymph_mean',xgb_fit_reg,train_data_impute$impreg)

# evaluate AUCs (as helper in xgb_utils.R)
xgb_auc_reg <- xgb_auc(xgb_fit_reg,dwatchlist_reg)
print(xgb_auc_reg)



