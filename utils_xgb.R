# function to take training, test and validation sets and 
# format for xgb fitting

xgb_prep <- function(train,test,valid,dname,cc=NULL,weight=NULL){
  # xgb formatting for each set
  dtrain <- xgb.DMatrix(as.matrix(train%>%select(-c(ID,CaseControl))),
                           label=train$CaseControl)
  dvalid <- xgb.DMatrix(as.matrix(valid%>%select(-c(ID,CaseControl))),
                           label=valid$CaseControl)
  dtest <- xgb.DMatrix(as.matrix(test%>%select(-c(ID,CaseControl))),
                          label=test$CaseControl)
  if(!missing(cc)) dcc = xgb.DMatrix(as.matrix(cc%>%select(-c(ID,CaseControl))),
                                    label=cc$CaseControl)
  if(!missing(weight)) setinfo(dtrain, "weight", weight)
  # combine as a watchlist
  dwatchlist <- list(train=dtrain,test=dtest,valid=dvalid)
  if(!missing(cc)) dwatchlist = list(train=dtrain,test=dtest,cc=dcc,valid=dvalid)
  # return
  return(dwatchlist)
}

# same formatiing but drop columns
xgb_prep_drop <- function(train,test,valid,dname,cc=NULL,drop){
  # xgb formatting for each set
  dtrain <- xgb.DMatrix(as.matrix(train[[dname]] %>% select(-one_of(c("ID","CaseControl",drop)))),
                        label=train[[dname]]$CaseControl)
  dvalid <- xgb.DMatrix(as.matrix(valid[[dname]] %>% select(-one_of(c("ID","CaseControl",drop)))),
                        label=valid[[dname]]$CaseControl)
  dtest <- xgb.DMatrix(as.matrix(test[[dname]] %>% select(-one_of(c("ID","CaseControl",drop)))),
                       label=test[[dname]]$CaseControl)
  if(!missing(cc)) dcc = xgb.DMatrix(as.matrix(cc %>% select(-one_of(c("ID","CaseControl",drop)))),
                                     label=cc$CaseControl)
  # combine as a watchlist
  dwatchlist <- list(train=dtrain,test=dtest,valid=dvalid)
  if(!missing(cc)) dwatchlist = list(train=dtrain,test=dtest,cc=dcc,valid=dvalid)
  # return
  return(dwatchlist)
}

# same formatting for only a subset of predictors
xgb_prep_sub <- function(train,test,valid,dname,cc=NULL,subset){
  # xgb formatting for each set
  dtrain <- xgb.DMatrix(as.matrix(train[[dname]][subset]),
                        label=train[[dname]]$CaseControl)
  dvalid <- xgb.DMatrix(as.matrix(valid[[dname]][subset]),
                        label=valid[[dname]]$CaseControl)
  dtest <- xgb.DMatrix(as.matrix(test[[dname]][subset]),
                       label=test[[dname]]$CaseControl)
  if(!missing(cc)) dcc = xgb.DMatrix(as.matrix(cc[subset]),
                                     label=cc$CaseControl)
  # combine as a watchlist
  dwatchlist <- list(train=dtrain,test=dtest,valid=dvalid)
  if(!missing(cc)) dwatchlist = list(train=dtrain,test=dtest,cc=dcc,valid=dvalid)
  # return
  return(dwatchlist)
}

# function to create a partial dependence plot for an xgboost model
xgb_pdp <- function(varname,
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
  ptest <- predict(xgb_model,newdata=xgb_data$test,
                            ntreelimit=xgb_model$best_ntreelimit,
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
  colnames(auc_mat) <- c('LB','Est','UB')
  return(auc_mat)
}

# function to evaluate the training and test AUCs with confidence intervals
xgb_auc_external <- function(xgb_model,
                    new_data){
  # convert to xgb object
  dnew <- xgb.DMatrix(as.matrix(new_data[-c(1,2)]),
                        label=new_data$CaseControl)
  # predict outsome
  ptest <- predict(xgb_model,newdata=dnew,
                    ntreelimit=xgb_model$best_ntreelimit,
                    outputmargin=FALSE)
  # AUCs
  auc_test <- ci.auc(response=new_data$CaseControl,
                     predictor=ptest,
                     conf.level=.95,
                     method='delong')
  # store in a matrix
  auc_mat <- matrix(auc_test,nrow=1)
  rownames(auc_mat) <- c('New Data')
  colnames(auc_mat) <- c('LB','Est','UB')
  return(auc_mat)
}

# get AUCs
best_auc = function(xgb_fit){
  i = xgb_fit$best_iteration
  best_auc = c(
    train = xgb_fit$evaluation_log$train_auc[i],
    valid = xgb_fit$evaluation_log$valid_auc[i],
    test = xgb_fit$evaluation_log$test_auc[i],
    cc = xgb_fit$evaluation_log$cc_auc[i]
  )
  return(best_auc)
}

