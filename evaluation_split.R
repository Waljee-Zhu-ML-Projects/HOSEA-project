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
    dfs[[paste0("(", bins[i]*100, ", ", bins[i+1]*100, "]% missing")]] = imputed[id, ]
  }
  out = list()
  for(dfname in names(dfs)){
    df0 = dfs[[dfname]]
    out[[dfname]] = xgb.DMatrix(as.matrix(df0[-c(1,2)]),
                                label=df0$CaseControl)
  }
  return(out)
}
