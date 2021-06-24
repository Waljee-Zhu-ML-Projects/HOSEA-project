# imputation scripts for linear/multiclass regression

# linear regression for numeric/non-negative/binary cases
impute_reg <- function(df,response,imiss,
                               nneg=FALSE,
                               binary=FALSE){
  # store vector to return
  out <- response
  # response
  if(nneg){
    yobs <- log(1e-2+response)
  }
  else{
    yobs <- response
  }
  # add response
  df[['y']] <- yobs
  # fit model
  temp_model <- lm(y~.,data=df,subset=!imiss)
  # impute
  suppressWarnings(yimp <- predict(temp_model,newdata=df[imiss,]))
  # replace
  if(binary){
    out[imiss] <- as.integer(yimp > .5)
  }
  else{
    if(nneg){
      out[imiss] <- pmax(exp(yimp)-1e-2,0)
    }
    else{
      out[imiss] <- yimp
    }
  }
  return(out)
}

# convergence issues for current data
impute_logit <- function(df,response,imiss){
  # store vector to return
  out <- response
  # add response
  df[['y']] <- response
  # fit model
  temp_model <- glm(y~.,data=df,family=binomial,subset=!imiss)
  # impute
  suppressWarnings(yimp <- predict(temp_model,newdata=df[imiss,],type='response'))
  # replace
  out[imiss] <- as.integer(yimp > .5)
  return(out)
}

# extension to multiclass variables
impute_reg_multi <- function(df,responses,imiss){
  # matrix to return
  out <- responses
  # impute each class
  classes <- colnames(responses)
  yimp <- matrix(NA,sum(imiss),length(classes))
  for(ii in 1:length(classes)){
     # run imputation
     yimp[,ii] <- impute_reg(df,responses[[classes[ii]]],imiss)[imiss]
  }
  # replace
  out[imiss,] <- 0
  wmiss <- which(imiss)
  for(ii in 1:sum(imiss)){
    if(any(yimp[ii,]>.5)){
      out[wmiss[ii],which.max(yimp[ii,])] <- 1
    }
  }
  return(out)
}



