
# # import data
# complete_data = readRDS('R_data/complete_data_raw.rds')
# complete_data$n_visits = NULL 
# 
# # subset for development
# df = complete_data[seq(100000), ]

# compute quantiles
compute_quantiles = function(df, n_quantiles){
  qs = seq(n_quantiles)/(n_quantiles + 1)
  quantiles_df = df %>%
    summarise_all(quantile, na.rm=T, probs=qs, type=1)
  quantiles_df$ID = NULL
  quantiles_df$CaseControl = NULL
  return(quantiles_df)
}

# impute
impute_srs = function(df, quantiles_df){
  n_quantiles = nrow(quantiles_df)
  for(col in colnames(quantiles_df)){
    ids = which(is.na(df[[col]]))
    n = length(ids)
    if(n>0){
      is = sample.int(n=n_quantiles, size=n, replace=T)
      values = quantiles_df[is, col]
      df[ids, col] = values
    }
  }
  return(df)
}

