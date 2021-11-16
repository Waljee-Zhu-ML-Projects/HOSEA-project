# for development
# xgb_fit = xgb_fit_resample
# df = dwatchlist_resample$valid
# threshold = seq(0.01, 0.99, 0.01)

classification_metrics = function(xgb_fit, df, threshold=0.5){
  threshold = c(threshold)
  threshold = sort(threshold)
  proba = predict(xgb_fit, newdata=df)
  y = xgboost::getinfo(df, "label")
  yhat = lapply(threshold, function(tr) as.numeric(proba>tr))
  # metrics depending on threshold
  cols = c("tpr", "tnr", "ppv", "npv", "detection_prevalance", 
           "accuracy", "balanced_accuracy", "F1_score", "mcc",
           "jaccard_index")
  metrics = data.frame(matrix(NA, nrow=length(threshold), ncol=length(cols)))
  rownames(metrics) = threshold
  colnames(metrics) = cols
  for(i in seq(length(threshold))){
    yy = yhat[[i]]
    N = length(y)
    p = sum(yy)
    n = sum(1-yy)
    tp = sum((yy*y))
    tn = sum((1-yy)*(1-y))
    fp = sum(yy*(1-y))
    fn = sum((1-yy)*y)
    metrics$ppv[i]                  = tp/p
    metrics$npv[i]                  = tn/n
    metrics$tpr[i]                  = tp/(tp+fn)
    metrics$tnr[i]                  = tn/(fp+tn)
    metrics$detection_prevalance[i] = p/N
    metrics$F1_score[i]             = 2 *tp/(2*tp+fp+fn)
    metrics$accuracy[i]             = (tp+tn)/N
    metrics$balanced_accuracy[i]    = (metrics$tpr[i] + metrics$tnr[i]) / 2
    metrics$mcc[i]                  = (tp*tn-fp*fn) / sqrt(p*n*(tp+fn)*(tn+fp))
    metrics$jaccard_index[i]        = (tp+tn) / (2*N-tp-tn)
    
  }
  
  return(metrics)
}

calibration = function(xgb_fit, df){
  proba = predict(xgb_fit, newdata=df)
  y = xgboost::getinfo(df, "label")
  # log loss, auroc, auprc, brier
  fg = proba[y==1]; bg = proba[y==0]
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  pr = PRROC::pr.curve(fg, bg ,curve=TRUE)
  return(list(roc=roc, pr=pr))
}
