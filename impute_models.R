# imputation scripts for linear,logistic,multinomial regression
# always exclude a smoke status variable 

# linear regression for numeric variables
impute_reg_numeric <- function(df,varname,imiss,
                               nneg=FALSE,verbose=FALSE){
  col <- df[[varname]]
  # fit a linear regression model
  if(nneg){
    model <- lm(as.formula(paste0('log(',varname,')~.')),
                data=df,
                subset=!imiss)
  }
  else{
    model <- lm(as.formula(paste0(varname,'~.')),
              data=df,
              subset=!imiss)
  }
  if(verbose){
    print('Fit model')
  }
  # calculate fitted values
  fitted <- predict(model,newdata=df[imiss,])
  if(verbose){
    print('Predict missing')
  }
  if(nneg){
    col[imiss] <- exp(fitted)
  }
  else{
    col[imiss] <- fitted
  }
  return(col)
}

# logistic regression for binary variables
impute_reg_binary <- function(df,varname,imiss,verbose=FALSE){
  col <- df[[varname]]
  # fit a linear regression model
  model <- glm(as.formula(paste0(varname,'~.')),
               family=binomial,
               data=df,
               subset=!imiss)
  if(verbose){
    print('Fit model')
  }
  # calculate fitted values
  fitted <- as.integer(predict(model,newdata=df[imiss,],type='response')>.5)
  if(verbose){
    print('Predict missing')
  }
  col[imiss] <- fitted
  return(col)
}

# multinomial regression for multiclass variables
# logistic regression for binary variables
impute_reg_multi <- function(df,varnames,imiss,verbose=FALSE){
  cols <- as.matrix(select(df,all_of(varnames)))
  cols_extend <- cbind(cols,1-rowSums(cols))
  
  preds <- as.matrix(select(df,-all_of(varnames)))
  # fit a linear regression model
  model <- multinom(cols_extend~preds,
                    subset=!imiss)
  if(verbose){
    print('Fit model')
  }
  # calculate fitted values
  fitted <- predict(model,newdata=preds)
  if(verbose){
    print('Predict missing')
  }
  # fill in columns
  wmiss <- which(imiss)
  cols[imiss,] <- 0
  for(ii in 1:sum(imiss)){
    if(as.integer(fitted[ii]) <= ncol(cols)){
      cols[wmiss[ii],fitted[ii]] <- 1
    }
  }
  if(verbose){
    print('Fill values')
  }
  return(cols)
}

# new simpler routine purely with linear algebra...
impute_simple_numeric <- function(df,response,imiss,
                               nneg=FALSE,
                               binary=FALSE,
                               ridge=0){
  # store vector to return
  out <- response
  # predictors
  Xobs <- cbind(1,as.matrix(df[!imiss,]))
  Xmiss <- cbind(1,as.matrix(df[imiss,]))
  # response
  if(nneg){
    yobs <- log(1e-2+response[!imiss])
  }
  else{
    yobs <- response[!imiss]
  }
  # impute
  yimp <- c(Xmiss %*% solve(crossprod(Xobs) + ridge*diag(ncol(Xobs)),crossprod(Xobs,yobs)))
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

impute_simple_multi <- function(df,responses,imiss,
                                ridge=0){
  # matrix to return
  out <- responses
  # predictors
  Xobs <- cbind(1,as.matrix(df[!imiss,]))
  Xmiss <- cbind(1,as.matrix(df[imiss,]))
  # response
  yobs <- as.matrix(responses[!imiss,])
  # impute
  yimp <- Xmiss %*% solve(crossprod(Xobs) + ridge*diag(ncol(Xobs)),crossprod(Xobs,yobs))
  # replace
  out[imiss,] <- 0
  wmiss <- which(imiss)
  for(ii in 1:sum(imiss)){
    if(any(yimp[ii,])>.5){
      out[wmiss[ii],which.max(yimp[ii,])] <- 1
    }
  }
  return(out)
}



