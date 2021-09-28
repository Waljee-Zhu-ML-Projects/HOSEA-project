# code to fit xgboost models to HOSEA data (subsample)

# starts from complete_data_raw (i.e. after link_summaries.R)
# split into test,train,validation
# impute independently with function
# prepare imputed data for xgboost
# fit models, inspect AUC, PDP plots, etc.

#### source helper functions ####
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/impute_models.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/impute_multisamp.R')

#### import data ####
cc_test <- readRDS('R_data/cc_complete_data.rds')
load('R_data/subsample/sub_complete_data_impute.RData')

#### format imputed data for xgboost ####
# with NAs
dwatchlist_na = xgb_prep(train_data_impute,
                          test_data_impute,
                          valid_data_impute,
                          dname='clean',
                         cc=cc_test)

# imputed with sampling
dwatchlist_samp = xgb_prep(train_data_impute,
                            test_data_impute,
                            valid_data_impute,
                            dname='impsamp',
                           cc=cc_test)

# imputed with median
dwatchlist_med = xgb_prep(train_data_impute,
                           test_data_impute,
                           valid_data_impute,
                           dname='impmed',
                          cc=cc_test)

# imputed with regression
dwatchlist_reg = xgb_prep(train_data_impute,
                           test_data_impute,
                           valid_data_impute,
                           dname='impreg',
                          cc=cc_test)

#### XGBoost parameters ####
param_xg = list(
  max_depth=4,
  eta=.05,
  objective='binary:logistic',
  eval_metric='auc'
)

# fit xgboost model, treating NA as extra bin
xgb_fit_na <- xgb.train(param_xg,
                        dwatchlist_na$train,
                        nrounds=1000,
                        dwatchlist_na,
                        verbose=1,
                        print_every_n=8,
                        early_stopping_rounds=10)

# fit xgboost model with random sampling imputation
xgb_fit_samp <- xgb.train(param_xg,
                          dwatchlist_samp$train,
                          nrounds=1000,
                          dwatchlist_samp,
                          verbose=1,print_every_n=8,
                          early_stopping_rounds=10)

# fit xgboost model with median imputation
xgb_fit_med <- xgb.train(param_xg,
                         dwatchlist_med$train,
                         nrounds=1000,
                         dwatchlist_med,
                         verbose=1,print_every_n=8,
                         early_stopping_rounds=10)

# fit xgboost model with CART imputation
xgb_fit_reg <- xgb.train(param_xg,
                         dwatchlist_reg$train,
                         nrounds=1000,
                         dwatchlist_reg,
                         verbose=1,print_every_n=8,
                         early_stopping_rounds=10)

#### Get AUC table ####

