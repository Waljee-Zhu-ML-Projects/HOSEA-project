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

# apply median imputation (y is target)
fill_by_median <- function(x, y=x) {
  y[is.na(y)] <- median(x, na.rm = TRUE)
  return(y)
}
# apply zero imputation
fill_by_zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}
# apply random sample imputation (y is target)
fill_by_sample <- function(x, y=x) {
  imiss_y <- is.na(y)
  nmiss <- sum(imiss_y)
  if(nmiss < 1) return(y)
  imiss_x <- is.na(x)
  y[imiss_y] <- sample(x[!imiss_x],nmiss,replace=TRUE)
  return(y)
}

fill_by_sample_mat <- function(x) {
  xmat <- as.matrix(x)
  isamp <- which(apply(xmat,1,function(y){all(!is.na(y))}))
  xsamp <- xmat[sample(isamp,nrow(xmat),replace=TRUE),]
  xmat[is.na(xmat)] <- xsamp[is.na(xmat)]
  return(as.data.frame(xmat))
}

# where first column is fully observed reference
fill_by_nn_mat <- function(x){
  xmat <- as.matrix(x)
  isamp <- which(apply(xmat,1,function(y){all(!is.na(y))}))
  for(ii in 1:nrow(xmat)){
    ireplace <- is.na(xmat[ii,])
    if(any(ireplace)){
      replacerow <- isamp[which.min(abs(xmat[isamp,1]-xmat[ii,1]))]
      xmat[ii,ireplace] <- xmat[replacerow,ireplace]
    }
  }
  return(as.data.frame(xmat))
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
              auc=auc_mat))
}

logistic_auc_external <- function(logistic_model,
                                  new_data){
  # predictions
  ptest <- predict(logistic_model,
                   newdata=new_data,
                   type='response')
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

# row-normalized and column normalized tables
table_rn <- function(x,y){
  temp <- table(x,y)
  temp_n <- temp / rowSums(temp)
  return(temp_n)
}

table_cn <- function(x,y){
  temp <- table(x,y)
  temp_n <- t(t(temp) / rowSums(t(temp)))
  return(temp_n)
}

# converting to a vector to its corresponding normal quantiles
quantile_normalize <- function(vec){
  obs <- !is.na(vec)
  temp <- rep(NA,length(vec))
  temp[obs] <- qnorm(rank(vec[obs])/(1+length(vec[obs])))
  return(temp)
}
# inverting
quantile_unnormalize <- function(vec,wrt){
  quantile(wrt,pnorm(vec),na.rm=TRUE,names=FALSE)
}





