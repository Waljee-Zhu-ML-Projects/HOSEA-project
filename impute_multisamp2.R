source("impute_missing.R")
source("utils_xgb")

xgb_multisamp_prog = function(
  train, test, valid, cc, 
  nreps, param_xg, mice_method="sample", ...
){
  # ============================================================================
  cat('\n=== IMPUTATION ===', fill=T)
  
  # split x/y
  train_y = select(train,all_of(c('ID','CaseControl')))
  train_x = select(train,-all_of(c('ID','CaseControl')))
  valid_y = select(valid,all_of(c('ID','CaseControl')))
  valid_x = select(valid,-all_of(c('ID','CaseControl')))
  test_y  = select(test,all_of(c('ID','CaseControl')))
  test_x  = select(test,-all_of(c('ID','CaseControl')))
  
  # merged to pass to mice
  merged_x = bind_rows(train_x, valid_x, test_x)
  train_n  = nrow(train_x)
  valid_n  = nrow(valid_x)
  test_n   = nrow(test_x)
  ignore   = rep(T, train_n + test_n + valid_n)
  ignore[1:train_n] = F
  
  # MICE call
  mice_result = mice::mice(merged_x, m=nreps, method=mice_method,
                           maxit=ifelse(mice_method=="sample", 1, 5),
                           visitSequence="monotone", ignore=ignore, ...)
  
  # ============================================================================
  cat('\n=== XGBoost ===', fill=T)
  
  xgb_model = NULL
  dcc    = xgb.DMatrix(select(cc,-all_of(c('ID','CaseControl'))), label=cc$CaseControl)
  
  for(i in seq(nreps)){
    cat(paste0('\n=== Round ', i, ' ==='), fill=T)
    # recover imputed train/test/valid X
    merged_imputed_x = mice::complete(mice_result, i)
    train_x = merged_imputed_x[1:train_n, ]
    valid_x = merged_imputed_x[(train_n+1):(train_n+valid_n), ]
    # create watchlist
    dtrain = xgb.DMatrix(train_x, label=train_y)
    dvalid = xgb.DMatrix(valid_x, label=valid_y)
    dwatchlist = list(train=dtrain, cc=dcc, valid=dvalid)
    # fit model
    xgb_model = xgb.train(
      param_xg,
      dwatchlist$train,
      nrounds=10000/nreps,
      dwatchlist,
      verbose=1,print_every_n=10,
      early_stopping_rounds=10,
      xgb_model=xgb_model
    )
  }
  # aucs
  i = xgb_model$best_iteration
  best_aucs = c(
    train = xgb_model$evaluation_log$train_auc[i],
    test  = 0.,
    cc    = xgb_model$evaluation_log$cc_auc[i],
  )
  
  # ============================================================================
  cat('\n=== Evaluation on test sets ===', fill=T)
  ptest <- matrix(0, nrow(test), nreps)
  for(i in seq(nreps)){
    merged_imputed_x = mice::complete(mice_result, i)
    test_x = merged_imputed_x[(train_n+valid_n+1):(train_n+valid_n+test_n), ]
    dtest = xgb.DMatrix(test_x, label=test_y)
    ptest[, i] = predict(mice_result, newdata=dtest)
  }
  ptest = apply(ptest, 1, mean)
  best_aucs["test"] = auc(test_y, ptest)
  
  return(best_aucs)
  
}









