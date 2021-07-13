# helper script to run multiple imputation plus xgboost

# inputs: 'clean' data split into train, test and validation; complete test set; 
# xgb parameters

# output: final xgb predictions on the 4 data sets 

# note: to save/predict one would save all the xgb models and run a test example through
# all of them, with newly sampled imputed values each time 

impute_multisamp <- function(data){
  # variable blocks
  # demographic (age,bmi,weight numeric, SmokeStatus multiclass, remaining 2-class)
  demo_vars_uni <- c('ageatindex','Gender',
                     'agentorange',
                     'GerdAtIndex',
                     'CHF','CTD','DEM','DIAB_C','HIV','MLD','MSLD','PARA','RD',
                     'cd','copd','diab_nc','mi','pud','pvd')
  demo_vars_mult <- list(c('bmi','weight'),
                         c('Asian','Black','HawaiianPacific','IndianAlaskan'),
                         c('smoke_current','smoke_former'))
  # blood lab measurements (all numeric)
  lab_vars <- list(c('A1c'),
                   c('bun','calc','chlor','co2','creat','gluc','k','na'),
                   c('baso','eos','hct','hgb','lymph','mch','mchc','mcv','mono','mpv','neut','platelet','rbc','rdw','wbc'),
                   c('CRP'),
                   c('alkphos','alt','ast','totprot'),
                   c('chol','hdl','ldl','trig'))
  # blood lab longitudinal summaries
  lab_suffix <- c('_mean',
                  '_max','_min',
                  '_max_diff','_min_diff','_tv')
  
  # impute
  for(lab in lab_vars){
    for(v in lab){
      varnames <- paste0(v,lab_suffix)
      data[,varnames] <- fill_by_sample_mat(data[,varnames])
    }
  }
  
  # univariate demographic variables
  for(varname in demo_vars_uni){
    data[[varname]] <- fill_by_sample(data[[varname]])
  }
  # multivariate demographic variables
  for(varnames in demo_vars_mult){
    data[,varnames] <- fill_by_sample_mat(data[,varnames])
  }
  return(data)
}

xgb_multisamp <- function(train,test,valid,cc_test,
                          nreps=1,
                          param_xg){
  
  # initialize predictions
  ptrain <- rep(0,nrow(train))
  ptest <- rep(0,nrow(test))
  pvalid <- rep(0,nrow(valid))
  pcctest <- rep(0,nrow(cc_test))
  
  for(nn in 1:nreps){
    print(paste0('For rep ',nn,':'))
    # impute
    temp_train <- temp_test <- temp_valid <- list()
    temp_train[['x']] <- impute_multisamp(train)
    temp_test[['x']] <- impute_multisamp(test)
    temp_valid[['x']] <- impute_multisamp(valid)
    print('Impute')
    # prep data
    temp_data <- xgb_prep(train,test,valid,dname='x')
    temp_cctestdata <- xgb.DMatrix(as.matrix(cc_test[-c(1,2)]),
                        label=cc_test$CaseControl)
    print('Prepare data')
    # fit
    temp_model <- xgb.train(param_xg,
                            temp_data$train, # training set
                            nrounds=500,
                            temp_data, # data watchlist
                            verbose=0,
                            early_stopping_rounds=10)
    print('Fit XGB')
    # print ntrees
    print(paste0(temp_model$best_ntreelimit,' trees'))
    # update predictions 
    ptrain <- ptrain + (1/nreps)*predict(temp_model,newdata=temp_data$train,
                                         ntreelimit=temp_model$best_ntreelimit,
                                         outputmargin=TRUE)
    ptest <- ptest + (1/nreps)*predict(temp_model,newdata=temp_data$test,
                                         ntreelimit=temp_model$best_ntreelimit,
                                         outputmargin=TRUE)
    pvalid <- pvalid + (1/nreps)*predict(temp_model,newdata=temp_data$valid,
                                         ntreelimit=temp_model$best_ntreelimit,
                                         outputmargin=TRUE)
    pcctest <- pcctest + (1/nreps)*predict(temp_model,newdata=temp_cctestdata,
                                           ntreelimit=temp_model$best_ntreelimit,
                                           outputmargin=TRUE)
    print('Update predictions')
  }
  
  # calculate final AUCs
  auc_train <- ci.auc(response=train$CaseControl,
                      predictor=expit(ptrain),
                      conf.level=.95,
                      method='delong')
  auc_valid <- ci.auc(response=test$CaseControl,
                      predictor=expit(pvalid),
                      conf.level=.95,
                      method='delong')
  auc_test <- ci.auc(response=valid$CaseControl,
                     predictor=expit(ptest),
                     conf.level=.95,
                     method='delong')
  auc_cctest <- ci.auc(response=cc_test$CaseControl,
                     predictor=expit(pcctest),
                     conf.level=.95,
                     method='delong')
  # store in a matrix
  auc_mat <- rbind(auc_train,auc_valid,auc_test,auc_cctest)
  rownames(auc_mat) <- c('Training','Validation','Test','Test (complete)')
  colnames(auc_mat) <- c('LB','Est','UB')
  
  return(list(ptrain=ptrain,
              ptest=ptest,
              pvalid=pvalid,
              pcctest=pcctest,
              aucs=auc_mat))
}