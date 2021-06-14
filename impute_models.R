# imputation scripts for linear,logistic,multinomial regression
# always exclude a smoke status variable 

# linear regression for numeric variables
impute_reg_numeric <- function(df,varname,imiss,verbose=FALSE){
  col <- df[[varname]]
  # fit a linear regression model
  model <- lm(as.formula(paste0(varname,'~.')),
              data=df,
              subset=!imiss)
  if(verbose){
    print('Fit model')
  }
  # calculate fitted values
  fitted <- predict(model,newdata=df[imiss,])
  if(verbose){
    print('Predict missing')
  }
  col[imiss] <- fitted
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
  cols <- select(df,all_of(varnames))
  cols_extend <- cbind(cols,1-rowSums(cols))
  
  preds <- select(df,-all_of(varnames))
  # fit a linear regression model
  model <- multinom(cols_extend~preds,
                    subset=!imiss)
  if(verbose){
    print('Fit model')
  }
  # calculate fitted values
  fitted <- predict(model,newdata=df[imiss,])
  if(verbose){
    print('Predict missing')
  }
  # fill in columns
  wmiss <- which(imiss)
  cols[imiss,] <- 0
  for(ii in 1:sum(imiss)){
    if(fitted[ii] < ncols(cols)){
      cols[wmiss[ii],fitted[ii]] <- 1
    }
  }
  if(verbose){
    print('Fill values')
  }
  return(cols)
}

# need similar routines for xgboost 




