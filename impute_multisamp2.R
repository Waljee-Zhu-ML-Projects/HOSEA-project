
xgb_multisamp_prog = function(
  train, test, valid, cc, 
  nreps, param_xg, mice_method="sample", ...
){
  # ============================================================================
  cat('\n=== IMPUTATION ===', fill=T)
  if(is.numeric(nreps)) nreps = c(nreps)
  
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
  mice_result = mice::mice(merged_x, m=max(nreps), method=mice_method,
                           maxit=ifelse(mice_method=="sample", 1, 5),
                           visitSequence="monotone", ignore=ignore, ...)
  
  # ============================================================================
  
  best_aucs = matrix(0, length(nreps), 3)
  rownames(best_aucs) = paste0("multsamp", nreps)
  colnames(best_aucs) = c("train", "test", "cc")
  dcc = xgb.DMatrix(as.matrix(select(cc,-all_of(c('ID','CaseControl')))), 
                    label=cc$CaseControl)
  
  for(m in nreps){
    cat(paste0('\n=== XGBoost(', m, ') ==='), fill=T)
    xgb_model = NULL
    for(i in seq(m)){
      cat(paste0('\n=== Round ', i, '/', m, ' ==='), fill=T)
      # recover imputed train/test/valid X
      merged_imputed_x = mice::complete(mice_result, i)
      train_x = merged_imputed_x[1:train_n, ]
      valid_x = merged_imputed_x[(train_n+1):(train_n+valid_n), ]
      # create watchlist
      dtrain = xgb.DMatrix(as.matrix(train_x), label=train_y$CaseControl)
      dvalid = xgb.DMatrix(as.matrix(valid_x), label=valid_y$CaseControl)
      dwatchlist = list(train=dtrain, cc=dcc, valid=dvalid)
      # fit model
      xgb_model = xgb.train(
        params=param_xg,
        data=dwatchlist$train,
        nrounds=5000/m,
        watchlist=dwatchlist,
        verbose=1,
        print_every_n=50,
        #early_stopping_rounds=50,
        xgb_model=xgb_model
      )
    }
    # aucs
    best_aucs[paste0("multsamp", m), "train"] = tail(xgb_model$evaluation_log$train_auc, n=1)
    best_aucs[paste0("multsamp", m), "cc"]    = tail(xgb_model$evaluation_log$cc_auc, n=1)
    
    # ============================================================================
    cat('\n=== Evaluation on test sets ===', fill=T)
    ptest = matrix(0, nrow(test), m)
    for(i in seq(m)){
      merged_imputed_x = mice::complete(mice_result, i)
      test_x = merged_imputed_x[(train_n+valid_n+1):(train_n+valid_n+test_n), ]
      dtest = xgb.DMatrix(as.matrix(test_x), label=test_y$CaseControl)
      ptest[, i] = predict(xgb_model, newdata=dtest)
    }
    ptest = apply(ptest, 1, mean)
    best_aucs[paste0("multsamp", m), "test"] = auc(test_y$CaseControl, ptest)
  }
  return(best_aucs)
  
}









