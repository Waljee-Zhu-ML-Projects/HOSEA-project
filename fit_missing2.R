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
source('R_code/hosea-project/impute_multisamp2.R')

#### import data ####
cc_test <- readRDS('R_data/cc_complete_data.rds')
load('R_data/subsample/sub_complete_data_impute.RData')


set.seed(1)

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
  subsample=0.5,
  eta=.05,
  objective='binary:logistic',
  eval_metric='auc'
)

# fit xgboost model, treating NA as extra bin
xgb_fit_na <- xgb.train(param_xg,
                        dwatchlist_na$train,
                        nrounds=10000,
                        dwatchlist_na,
                        verbose=1,
                        print_every_n=10,
                        early_stopping_rounds=50)

# fit xgboost model with random sampling imputation
xgb_fit_samp <- xgb.train(param_xg,
                          dwatchlist_samp$train,
                          nrounds=10000,
                          dwatchlist_samp,
                          verbose=1,print_every_n=10,
                          early_stopping_rounds=50)

# fit xgboost model with median imputation
xgb_fit_med <- xgb.train(param_xg,
                         dwatchlist_med$train,
                         nrounds=10000,
                         dwatchlist_med,
                         verbose=1,print_every_n=10,
                         early_stopping_rounds=50)

# fit xgboost model with CART imputation
xgb_fit_reg <- xgb.train(param_xg,
                         dwatchlist_reg$train,
                         nrounds=10000,
                         dwatchlist_reg,
                         verbose=1,print_every_n=10,
                         early_stopping_rounds=50)

# fit xgboost model with multiple samples imputation

xgb_multisamp = xgb_multisamp_prog(
  train=train_data_impute$clean,
  valid=valid_data_impute$clean,
  test=test_data_impute$clean,
  cc=cc_test,
  nreps=c(1, 2, 5, 10, 20, 50, 100),
  param_xg,
  mice_method="sample"
)

#### Get AUC table ####

best_auc = function(xgb_fit){
  i = xgb_fit$best_iteration
  best_auc = c(
    train = xgb_fit$evaluation_log$train_auc[i],
    test = xgb_fit$evaluation_log$test_auc[i],
    cc = xgb_fit$evaluation_log$cc_auc[i]
  )
  return(best_auc)
}

best_aucs = rbind(
  na = best_auc(xgb_fit_na),
  med = best_auc(xgb_fit_med),
  reg = best_auc(xgb_fit_reg),
  samp = best_auc(xgb_fit_samp),
  xgb_multisamp
)

write.csv(best_aucs, "R_data/results/best_auc_4subsample.csv")

best_aucs = read.csv("R_data/results/best_auc_4subsample.csv")
rownames(best_aucs) = best_aucs$X
best_aucs$X = NULL
xtable::xtable(best_aucs, digits=3)
