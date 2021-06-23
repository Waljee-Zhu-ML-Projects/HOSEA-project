mmean <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    mean(x,na.rm=TRUE)
  }
}

mmax <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    max(x,na.rm=TRUE)
  }
}

mmin <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    min(x,na.rm=TRUE)
  }
}

ssum <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    sum(x,na.rm=TRUE)
  }
}

# apply median imputation
fill_by_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}
# apply zero imputation
fill_by_zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}
# apply random sample imputation
fill_by_sample <- function(x) {
  imiss <- is.na(x)
  nmiss <- sum(imiss)
  x[imiss] <- sample(x[!imiss],nmiss,replace=TRUE)
  return(x)
}

expit <- function(x){
  exp(x) / (1+exp(x))
}

# logistic fitting and AUC recovery
fit_logistic <- function(train_data,test_data){
  # fit model
  model <- glm(CaseControl ~ .,
               data=train_data,
               family=binomial)
  # get predictions
  ptrain <- predict(model,
                   newdata=train_data,
                   type='response')
  ptest <- predict(model,
                    newdata=test_data,
                    type='response')
  
  auc_train <- ci.auc(response=train_data$CaseControl,
                      predictor=ptrain,
                      conf.level=.95,
                      method='delong')
  auc_test <- ci.auc(response=test_data$CaseControl,
                     predictor=ptest,
                     conf.level=.95,
                     method='delong')
  # store in a matrix
  auc_mat <- rbind(auc_train,auc_test)
  rownames(auc_mat) <- c('Training','Test')
  colnames(auc_mat) <- c('LB','Est','UB')
  # return
  return(list(model=model,
              preds=preds,
              auc=auc_mat))
}


