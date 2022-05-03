# depends on global variables 'complete_data'
# complete_data <- readRDS('R_data/complete_data_raw.rds')

subsample_controls = function(df=complete_data, n_controls=1e6){
  i_case = which(df$casecontrol==1)
  i_control = sample(which(df$casecontrol==0),n_controls)
  i = c(i_case, i_control)  
  sub_master = df[i,]
  sub_ids = sub_master$id
  sub = filter(df, id %in% sub_ids)
  return(sub)
}

subsample = function(df=complete_data, n_samples=10000, balanced=TRUE){
  if(balanced){
     n0 = floor(n_samples/2)
     n1 = n_samples - n0
     i_case = sample(which(df$casecontrol==1),n1)
     i_control = sample(which(df$casecontrol==0),n0)
     i = c(i_case, i_control)
  }else{
    i = sample(df$id, n_samples)
  }
  sub_master = df[i,]
  sub_ids = sub_master$id
  sub = filter(df, id %in% sub_ids)
  return(sub)
}

balanced_resample = function(df){
  # assumes more controls than cases
  i1 = which(df$casecontrol==1)
  n0 = sum(df$casecontrol==0)
  n1 = length(i1)
  n = max(0, n0-n1)
  df1 = df[sample(i1, n, replace=T),]
  return(rbind(df, df1))
}

train_test_split = function(df=complete_data, weights=c(1), stratified=TRUE){
  k = length(weights)
  if(stratified){
    dfs = list()
    for(v in c(0,1)){
      i0 = sample(which(df$casecontrol==v)) 
      n0 = length(i0)
      ns0 = round(n0 * weights / sum(weights))
      ns0[1] = n0 - sum(ns0[-1])
      bounds = cumsum(c(0, ns0))
      dfs[[v+1]] = lapply(seq(k), function(i) {
        df[i0[(bounds[i]+1):bounds[i+1]], ]
        })
    }
    # merge
    out = list()
    for(i in seq(k)) out[[i]] = rbind(dfs[[1]][[i]], dfs[[2]][[i]])
  }else{
    out = list()
    i0 = sample(df$id) 
    n0 = length(i0)
    ns0 = round(n0 * weights / sum(weights))
    ns0[1] = n0 - sum(ns0[-1])
    bounds = cumsum(c(0, ns0))
    out = lapply(seq(k), function(i) {
      filter(df, id %in% i0[(bounds[i]+1):bounds[i+1]])
    })
  }
  
  lapply(out, function(df) df$casecontrol %>% mean)
  return(out)
}