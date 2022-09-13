setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(tidyr)
library(HOSEA)
library(ggplot2)
theme_set(theme_minimal())
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_data = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
dir_models = "R_data/results/models/"

# =========================================================
# read in models
models = list(
  pre_egjac=readRDS(paste0(dir_models, "old/XGB_all_EAC.rds")),
  post_egjac=readRDS(paste0(dir_models, "XGB_all_EAC.rds"))
)

# =========================================================
# read in data

# data before the last update had min-max and mindiff-maxdiff switched in lab variables
switch_minmax = function(df){
  cols = df$df %>% colnames()
  cols = sapply(cols, function(col){
    if(grepl("h2r", col, fixed=TRUE)) return(col)
    if(grepl("ppi", col, fixed=TRUE)) return(col)
    if(grepl("min", col, fixed=TRUE)) return(stringr::str_replace(col, "min", "max"))
    if(grepl("max", col, fixed=TRUE)) return(stringr::str_replace(col, "max", "min"))
    return(col)
  })
  colnames(df$df) = cols
  return(df)
}

dfs = list(
  post_egjac=readRDS(paste0(dir_data, "5-1_merged.rds")),
  pre_egjac=readRDS(paste0(dir_data, "5-1_merged_old.rds"))
)


# =========================================================
# get predicted risks
outcome = "EAC"

probas = lapply(names(dfs), function(name){
  df = dfs[[name]]$df
  xgb_fit = models[[name]]$xgb_fit
  quantiles = models[[name]]$quantiles
  test_ids = models[[name]]$test_ids
  
  master = dfs[[name]]$master
  df1 = df %>% filter(id %in% test_ids)
  set.seed(0)
  df1 %<>% impute_srs(quantiles)
  
  outcomes = master %>% select(id, cancertype)
  outcomes %<>% mutate(
    ANY=ifelse(cancertype=="", 0, 1),
    EAC=ifelse(cancertype=="EAC", 1, 0),
    EGJAC=ifelse(cancertype=="EGJAC", 1, 0)
  )
  df1 %<>% left_join(outcomes, by="id")
  df1 %<>% mutate(casecontrol = get(outcome))
  xgb_df = xgb.DMatrix(as.matrix(df1 %>% select(xgb_fit$feature_names)),
                   label=df1$casecontrol)
  proba = predict(xgb_fit, newdata=xgb_df)
  out = data.frame(id=df1$id, name=proba)
  colnames(out) = c("id", name)
  return(out)
})

probas %<>% purrr::reduce(full_join, by="id")

probas %<>%
  mutate(
    diff=post_egjac - pre_egjac,
    ratio=post_egjac / pre_egjac,
    logratio=log(post_egjac) - log(pre_egjac) 
  )

# =========================================================
# inspect differences

ggplot(
  data=probas %>% sample_n(10000),
  mapping=aes(x=pre_egjac, y=post_egjac)
) + geom_point(alpha=0.1) + 
  geom_abline(intercept=0, slope=1) + 
  scale_x_log10() + scale_y_log10()

ggplot(
  data=probas,
  mapping=aes(x=logratio)
) + geom_histogram(bins=1000)
