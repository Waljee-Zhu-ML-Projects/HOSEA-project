# code to compare model with different subsets of variables
# fitting with single random sample imputation

#### source helper functions ####
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_xgb.R')

#### import data ####

load('R_data/subsample/sub_complete_data_impute.RData')
# named lists of data frames: train/test/valid_data_impute
# load complete case test set
cc_test <- readRDS('R_data/cc_complete_data.rds')

#### subsets of variables ####

# always keep columns 1+2: ID and CaseControl indicator
demo_vars <- c(3:11,28:29)
charlson_vars <- 12:27
med_vars <- 34:43
lab_vars <- c(28:33,44:241)

dwatchlist_demo <- xgb_prep_sub(train_data_impute,
                                test_data_impute,
                                valid_data_impute,
                                'impsamp',
                                subset=demo_vars)

dwatchlist_charlson <- xgb_prep_sub(train_data_impute,
                                test_data_impute,
                                valid_data_impute,
                                'impsamp',
                                subset=charlson_vars)

dwatchlist_med <- xgb_prep_sub(train_data_impute,
                                    test_data_impute,
                                    valid_data_impute,
                                    'impsamp',
                                    subset=med_vars)

dwatchlist_lab <- xgb_prep_sub(train_data_impute,
                               test_data_impute,
                               valid_data_impute,
                               'impsamp',
                               subset=lab_vars)

#### fitting and AUCs ####

# global parameters for XGB
param_xg <- list(max_depth=4,eta=.05,objective='binary:logistic',
                 eval_metric='logloss')

# Demographic vars
# fit xgboost models, get AUCs 
xgb_fit_demo <- xgb.train(param_xg,
                        dwatchlist_demo$train, # training set
                        nrounds=1000,
                        dwatchlist_demo, # data watchlist
                        verbose=1,print_every_n=8,
                        early_stopping_rounds=10)

# evaluate AUCs
xgb_auc_demo <- xgb_auc(xgb_fit_demo,dwatchlist_demo)
print(xgb_auc_demo)
# evaluate AUC on complete cases
xgb_auc_demo_cc <- xgb_auc_external(xgb_fit_demo,cc_test[,c(1,2,demo_vars)])
print(xgb_auc_demo_cc) 

# Charlson vars
# fit xgboost models, get AUCs 
xgb_fit_charlson <- xgb.train(param_xg,
                          dwatchlist_charlson$train, # training set
                          nrounds=1000,
                          dwatchlist_charlson, # data watchlist
                          verbose=1,print_every_n=8,
                          early_stopping_rounds=10)

# evaluate AUCs
xgb_auc_charlson <- xgb_auc(xgb_fit_charlson,dwatchlist_charlson)
print(xgb_auc_charlson)
# evaluate AUC on complete cases
xgb_auc_charlson_cc <- xgb_auc_external(xgb_fit_charlson,cc_test[,c(1,2,charlson_vars)])
print(xgb_auc_charlson_cc) 

# Medication vars
# fit xgboost models, get AUCs 
xgb_fit_med <- xgb.train(param_xg,
                          dwatchlist_med$train, # training set
                          nrounds=1000,
                          dwatchlist_med, # data watchlist
                          verbose=1,print_every_n=8,
                          early_stopping_rounds=10)

# evaluate AUCs
xgb_auc_med <- xgb_auc(xgb_fit_med,dwatchlist_med)
print(xgb_auc_med)
# evaluate AUC on complete cases
xgb_auc_med_cc <- xgb_auc_external(xgb_fit_med,cc_test[,c(1,2,med_vars)])
print(xgb_auc_med_cc) 

# Blood lab vars
# fit xgboost models, get AUCs 
xgb_fit_lab <- xgb.train(param_xg,
                          dwatchlist_lab$train, # training set
                          nrounds=1000,
                          dwatchlist_lab, # data watchlist
                          verbose=1,print_every_n=8,
                          early_stopping_rounds=10)

# evaluate AUCs
xgb_auc_lab <- xgb_auc(xgb_fit_lab,dwatchlist_lab)
print(xgb_auc_lab)
# evaluate AUC on complete cases
xgb_auc_lab_cc <- xgb_auc_external(xgb_fit_lab,cc_test[,c(1,2,lab_vars)])
print(xgb_auc_lab_cc) 

