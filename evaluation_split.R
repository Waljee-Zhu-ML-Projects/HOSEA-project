# for development
# df = test
# bins=c(0, 0.05, 0.1, 0.3, 1.0)
# quantile(na_prop)

evaluation_split = function(df, imputed, bins=c(0, 0.05, 0.1, 0.3, 1.0)){
  ids = df$ID
  na_prop = rowMeans(is.na(select(df, -c(ID, CaseControl))))
  cc = na_prop == 0.
  dfs = list(all=imputed, cc=imputed[cc, ])
  for(i in seq(length(bins)-1)){
    id = (na_prop>bins[i]) & (na_prop<=bins[i+1])
    dfs[[paste0("all_", bins[i]*100, "_", bins[i+1]*100)]] = imputed[id, ]
  }
  out = list()
  for(dfname in names(dfs)){
    df0 = dfs[[dfname]]
    out[[dfname]] = xgb.DMatrix(as.matrix(df0[-c(1,2)]),
                                label=df0$CaseControl)
  }
  return(out)
}


split_by_vargoups = function(test, test_){
  lab_vars0 = unlist(lab_vars)
  lab_vars0 = paste0(sort(rep(lab_vars0, length(lab_suffix))), lab_suffix)
  vargroups = list(lab_vars=lab_vars0, charlson_vars=charlson_vars, demo_vars=demo_vars,
                   other_vars=other_vars, smoke_vars=smoke_vars)
  vargroups_bins = list(lab_vars=c(0.00, 0.06, 0.14, 0.33, 1.00, 1.01), 
                        charlson_vars=c(0.00, 0.01, 1.01), 
                        demo_vars=c(0.00, 0.01, 1.01),
                        other_vars=c(0.00, 1.01), 
                        smoke_vars=c(0.00, 0.50, 1.01))
  dfs = list()
  for(vargroup in names(vargroups)){
    cat(vargroup, fill=T)
    vars = vargroups[[vargroup]]
    dfraw = test %>% select(one_of(vars))
    na_prop = rowMeans(is.na(dfraw))
    bins = vargroups_bins[[vargroup]]
    for(i in seq(length(bins)-1)){
      id = (na_prop>=bins[i]) & (na_prop<bins[i+1])
      dfs[[paste0(vargroup, "_", bins[i]*100, "_", bins[i+1]*100)]] = test_[id, ]
    }
  }
  out = lapply(dfs, function(df) xgb.DMatrix(as.matrix(df[-c(1,2)]),
                                             label=df$CaseControl))
  return(out)
}
