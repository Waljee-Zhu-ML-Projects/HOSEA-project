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
source('R_code/hosea-project/impute_multisamp.R')

#### import data ####

complete_data <- readRDS('R_data/subsample/sub_complete_data_raw.rds')
# load complete case test set
cc_test <- readRDS('R_data/subsample/cc_test.rds')
cc_test <- cc_test[,c(1:33,35,34,36:38,40,39,41:241)]
colnames(cc_test) <- colnames(complete_data)

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
# 10 cycles
train_data_impute <- impute_missing_hosea(complete_data[traintrain,],ncycles=10,seed=1995,hybrid_reg=TRUE)
test_data_impute <- impute_missing_hosea(complete_data[test,],ncycles=10,seed=1996,hybrid_reg=TRUE)
valid_data_impute <- impute_missing_hosea(complete_data[valid,],ncycles=10,seed=1998,hybrid_reg=TRUE)

# save all imputed data
save(train_data_impute,test_data_impute,valid_data_impute,
     file='R_data/subsample/sub_complete_data_impute_v03.RData')
# v02 with renormalized regression imputation
# can reload from here
#load('R_data/subsample/sub_complete_data_impute.RData')

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

# imputed with median
dwatchlist_med <- xgb_prep(train_data_impute,
                           test_data_impute,
                           valid_data_impute,
                           dname='impmed')

# imputed with regression
dwatchlist_reg <- xgb_prep(train_data_impute,
                          test_data_impute,
                          valid_data_impute,
                          dname='impreg')

#### fit with NAs ####

# global parameters for the 3 models
param_xg <- list(max_depth=4,eta=.05,objective='binary:logistic',
                 eval_metric='logloss')
# think about tuning max_depth,eta,nrounds with CV?

# fit xgboost model
xgb_fit_na <- xgb.train(param_xg,
                     dwatchlist_na$train, # training set
                     nrounds=1000,
                     dwatchlist_na, # data watchlist
                     verbose=1,print_every_n=8,
                     early_stopping_rounds=10)
# stops at 202 trees

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
# inflated .942 test AUC with missing values
# evaluate AUC on complete cases
xgb_auc_na_cc <- xgb_auc_external(xgb_fit_na,cc_test)
print(xgb_auc_na_cc) # 0.643

#### fit with random sample imputation ####

# baseline logistic
logistic_fit_samp <- fit_logistic(train_data = rbind(train_data_impute$impsamp[,-1],
                                                     valid_data_impute$impsamp[,-1]),
                                  test_data = test_data_impute$impsamp[,-1])
print(logistic_fit_samp$auc)
# honest .689 test AUC with logistic regression
logistic_auc_samp_cc <- logistic_auc_external(logistic_fit_samp$model,
                                              cc_test)
print(logistic_auc_samp_cc) # 0.724 

# fit xgboost model
xgb_fit_samp <- xgb.train(param_xg,
                        dwatchlist_samp$train, # training set
                        nrounds=1000,
                        dwatchlist_samp, # data watchlist
                        verbose=1,print_every_n=8,
                        early_stopping_rounds=10)
# stops at 169 trees

# print info about important variables
xgb_summ_samp <- xgb.importance(model=xgb_fit_samp)
print(xgb_summ_samp[200:219,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrixi
# age
xgb_pdp('ageatindex',xgb_fit_samp,train_data_impute$impsamp)
# weight
xgb_pdp('weight',xgb_fit_samp,train_data_impute$impsamp)
xgb_pdp('pud',xgb_fit_samp,train_data_impute$impsamp)
# labs: WBC, BUN, Na (min)
xgb_pdp('wbc_mean',xgb_fit_samp,train_data_impute$impsamp)
xgb_pdp('bun_mean',xgb_fit_samp,train_data_impute$impsamp)
xgb_pdp('na_min',xgb_fit_samp,train_data_impute$impsamp)

# evaluate AUCs (as helper in xgb_utils.R)
xgb_auc_samp <- xgb_auc(xgb_fit_samp,dwatchlist_samp)
print(xgb_auc_samp)
# honest .729 test AUC with sampling
# evaluate AUC on complete cases
xgb_auc_samp_cc <- xgb_auc_external(xgb_fit_samp,cc_test)
print(xgb_auc_samp_cc) # 0.756

#### fit with median imputation ####

# baseline logistic
logistic_fit_med <- fit_logistic(train_data = rbind(train_data_impute$impmed[,-1],
                                                     valid_data_impute$impmed[,-1]),
                                  test_data = test_data_impute$impmed[,-1])
print(logistic_fit_med$auc)
# .810 test AUC with logistic regression + median imputation
logistic_auc_med_cc <- logistic_auc_external(logistic_fit_med$model,
                                              cc_test)
print(logistic_auc_med_cc) # 0.585 significantly worse, negative effects of Charlson?


# fit xgboost model
xgb_fit_med <- xgb.train(param_xg,
                          dwatchlist_med$train, # training set
                          nrounds=1000,
                          dwatchlist_med, # data watchlist
                          verbose=1,print_every_n=8,
                          early_stopping_rounds=10)
# stops at 319 trees

# print info about important variables
xgb_summ_med <- xgb.importance(model=xgb_fit_med)
print(xgb_summ_med[1:nsumm,])

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
# age
xgb_pdp('ageatindex',xgb_fit_med,train_data_impute$impmed)
# labs
xgb_pdp('hct_mean',xgb_fit_med,train_data_impute$impmed) # spiked at median
xgb_pdp('mch_mean',xgb_fit_med,train_data_impute$impmed)
xgb_pdp('A1c_max',xgb_fit_med,train_data_impute$impmed)

# evaluate AUCs (as helper in xgb_utils.R)
xgb_auc_med <- xgb_auc(xgb_fit_med,dwatchlist_med)
print(xgb_auc_med)
# .857 test AUC with median
# evaluate AUC on complete cases
xgb_auc_med_cc <- xgb_auc_external(xgb_fit_med,cc_test)
print(xgb_auc_med_cc) # 0.623

#### fit with regression imputation ####

# baseline logistic
logistic_fit_reg <- fit_logistic(train_data = rbind(train_data_impute$impreg[,-1],
                                                     valid_data_impute$impreg[,-1]),
                                  test_data = test_data_impute$impreg[,-1])
print(logistic_fit_reg$auc)
# .760 test AUC - better than baseline, worse than median, still using charlson?
# charlson gone with distn matching, still poor test performance?
logistic_auc_reg_cc <- logistic_auc_external(logistic_fit_reg$model,
                                             cc_test)
print(logistic_auc_reg_cc) # 0.557 significantly worse, again negative effects of Charlson?
# 0.567 with distn matching, unclear why that performs so poorly
# 0.657 with distn matching + hybrid (may have been fitting spurious patterns to longitudinal vars)

# fit xgboost model
xgb_fit_reg <- xgb.train(param_xg,
                        dwatchlist_reg$train, # training set
                        nrounds=1000,
                        dwatchlist_reg, # data watchlist
                        verbose=1,print_every_n=8,
                        early_stopping_rounds=10)
# stops at 362 trees
# 228 with distn matching
# 147 with distn matching + hybrid

# print info about important variables
xgb_summ_reg <- xgb.importance(model=xgb_fit_reg)
print(xgb_summ_reg[1:nsumm,])
# still some weird patterns with hct_max, etc.

# pdp plots (as helper in xgb_utils.R)
# need to pass raw training data matrix
# age
xgb_pdp('ageatindex',xgb_fit_reg,train_data_impute$impreg)
xgb_pdp('A1c_mean',xgb_fit_reg,train_data_impute$impreg)
xgb_pdp('hct_max',xgb_fit_reg,train_data_impute$impreg)
xgb_pdp('RD',xgb_fit_reg,train_data_impute$impreg)

# evaluate AUCs (as helper in xgb_utils.R)
xgb_auc_reg <- xgb_auc(xgb_fit_reg,dwatchlist_reg)
print(xgb_auc_reg)
# .815 test AUC (possibly inflated)
# .725 after distn matching
# evaluate AUC on complete cases
xgb_auc_reg_cc <- xgb_auc_external(xgb_fit_reg,cc_test)
print(xgb_auc_reg_cc) # 0.705, significantly worse but approaching
# random sample
# 0.757 with distribution matching, performance similar to random sample
# still 0.757 with distn matching + hybrid, again similar to RS

#### fit with multiple sample imputation ####

# # with 10 reps
# set.seed(2003)
# xgb_fit_multisamp10 <- xgb_multisamp(train_data_impute$clean,
#                                      test_data_impute$clean,
#                                      valid_data_impute$clean,
#                                      cc_test,
#                                      nreps=10,
#                                      param_xg=param_xg)
# 
# # with 20 reps
# set.seed(2004)
# xgb_fit_multisamp20 <- xgb_multisamp(train_data_impute$clean,
#                                      test_data_impute$clean,
#                                      valid_data_impute$clean,
#                                      cc_test,
#                                      nreps=20,
#                                      param_xg=param_xg)
# 
# # with 30 reps
# set.seed(2005)
# xgb_fit_multisamp30 <- xgb_multisamp(train_data_impute$clean,
#                                      test_data_impute$clean,
#                                      valid_data_impute$clean,
#                                      cc_test,
#                                      nreps=30,
#                                      param_xg=param_xg)

#### fit with multiple sample progressive #### 

# with 10 reps
set.seed(2006)
xgb_fit_multisampprog10 <- xgb_multisamp_prog(train_data_impute$clean,
                                     test_data_impute$clean,
                                     valid_data_impute$clean,
                                     cc_test,
                                     nreps=10,nrounds=50,
                                     param_xg=param_xg)
# see trace of training process by plotting model evaluation_log
matplot(1:500,xgb_fit_multisampprog10$model$evaluation_log[,-1],
        xlab='Trees',ylab='Loss',main='Training for 10 reps, 50 trees each',
        type='p',pch=1,lty=1)
for(ii in 1:10){
  abline(v=50*ii,lty=2)
}

# print AUCs
print(xgb_fit_multisampprog10$aucs)
# test AUC on complete records is 0.781

# with 20 reps
set.seed(2007)
xgb_fit_multisampprog20 <- xgb_multisamp_prog(train_data_impute$clean,
                                              test_data_impute$clean,
                                              valid_data_impute$clean,
                                              cc_test,
                                              nreps=20,nrounds=25,
                                              param_xg=param_xg)
# see trace of training process by plotting model evaluation_log
matplot(1:500,xgb_fit_multisampprog20$model$evaluation_log[,-1],
        xlab='Trees',ylab='Loss',main='Training for 20 reps, 25 trees each',
        type='p',pch=1,lty=1)
for(ii in 1:20){
  abline(v=25*ii,lty=2)
}

# print AUCs
print(xgb_fit_multisampprog20$aucs)
# test AUC on complete records is 0.778


# with 30 reps
set.seed(2008)
xgb_fit_multisampprog30 <- xgb_multisamp_prog(train_data_impute$clean,
                                              test_data_impute$clean,
                                              valid_data_impute$clean,
                                              cc_test,
                                              nreps=30,nrounds=20,
                                              param_xg=param_xg)
# see trace of training process by plotting model evaluation_log
matplot(1:600,xgb_fit_multisampprog30$model$evaluation_log[,-1],
        xlab='Trees',ylab='Loss',main='Training for 30 reps, 20 trees each',
        type='p',pch=1,lty=1)
for(ii in 1:30){
  abline(v=20*ii,lty=2)
}

# print AUCs
print(xgb_fit_multisampprog30$aucs)
# test AUC on complete records is 0.778

# with 100 reps
set.seed(2009)
xgb_fit_multisampprog100 <- xgb_multisamp_prog(train_data_impute$clean,
                                              test_data_impute$clean,
                                              valid_data_impute$clean,
                                              cc_test,
                                              nreps=100,nrounds=5,
                                              param_xg=param_xg)
# see trace of training process by plotting model evaluation_log
matplot(1:500,xgb_fit_multisampprog100$model$evaluation_log[,-1],
        xlab='Trees',ylab='Loss',main='Training for 100 reps, 5 trees each',
        type='p',pch=1,lty=1)
for(ii in 1:100){
  abline(v=5*ii,lty=2)
}

# variable importance
xgb_summ_multisamp <- xgb.importance(model=xgb_fit_multisampprog100$model)
print(xgb_summ_multisamp)
# look at partial dependence plots
xgb_pdp('hct_mean',xgb_fit_multisampprog100$model,train_data_impute$impsamp)
xgb_pdp('bun_mean',xgb_fit_multisampprog100$model,train_data_impute$impsamp)
xgb_pdp('mch_tv',xgb_fit_multisampprog100$model,train_data_impute$impsamp)


# print AUCs
print(xgb_fit_multisampprog100$aucs)
# test AUC on complete records is 0.806

# save all the multiple sample imputation models/results
# save(xgb_fit_multisamp10,
#      xgb_fit_multisamp20,
#      xgb_fit_multisamp30,
#      xgb_fit_multisampprog10,
#      xgb_fit_multisampprog20,
#      xgb_fit_multisampprog30,
#      file='R_data/models/multisamp_models.RData')

