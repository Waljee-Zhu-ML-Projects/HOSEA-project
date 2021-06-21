# function to take training, test and validation sets and 
# format for xgb fitting

xgb_prep <- function(train,test,valid,dname){
  # xgb formatting for each set
  dtrain <- xgb.DMatrix(as.matrix(train[[dname]][-c(1,2)]),
                           label=train[[dname]]$CaseControl)
  dvalid <- xgb.DMatrix(as.matrix(valid[[dname]][-c(1,2)]),
                           label=valid[[dname]]$CaseControl)
  dtest <- xgb.DMatrix(as.matrix(test[[dname]][-c(1,2)]),
                          label=test[[dname]]$CaseControl)
  # combine as a watchlist
  dwatchlist <- list(train=dtrain,test=dtest,valid=dvalid)
  # return
  return(dwatchlist)
}

# function to create a partial dependence plot for an xgboost model
xgp_pdp <- function(varname,
                    xgb_model,
                    dtrain){
  partial_var <- pdp::partial(xgb_model,
                              ntreelimit=xgb_model$best_ntreelimit,
                              train=dtrain[,-c(1,2)],
                              pred.var=varname,
                              grid.resolution=NULL,
                              quantiles=T,
                              probs=1:99/100,
                              plot=F)
  partial_var$yhat <- expit(partial_var$yhat)
  pdp::plotPartial(partial_var,rug=T,train=dtrain[,-1],
                   main=paste0('Partial dependence, ',varname),
                   ylab='prob')
}

# function to evaluate the training and test AUCs with confidence intervals
xgb_auc <- function(xgb_model,
                    xgb_data){
  # predict outsome
  ptrain <- predict(xgb_model,newdata=xgb_data$train,
                             ntreelimit=xgb_model$best_ntreelimit,
                             outputmargin=FALSE)
  pvalid <- predict(xgb_model,newdata=xgb_data$valid,
                             ntreelimit=xgb_model$best_ntreelimit,
                             outputmargin=FALSE)
  ptest <- predict(xgb_model,newdata=xgb_data$valid,
                            ntreelimit=xgb_fit$best_ntreelimit,
                            outputmargin=FALSE)
  # AUCs
  auc_train <- ci.auc(response=getinfo(xgb_data$train,name='label'),
                      predictor=ptrain,
                      conf.level=.95,
                      method='delong')
  auc_valid <- ci.auc(response=getinfo(xgb_data$valid,name='label'),
                      predictor=pvalid,
                      conf.level=.95,
                      method='delong')
  auc_test <- ci.auc(response=getinfo(xgb_data$test,name='label'),
                      predictor=ptest,
                      conf.level=.95,
                      method='delong')
  # store in a matrix
  auc_mat <- rbind(auc_train,auc_valid,auc_test)
  rownames(auc_mat) <- c('Training','Validation','Test')
  colnames(auc_mat) <- c('Est','LB','UB')
  return(auc_mat)
}

