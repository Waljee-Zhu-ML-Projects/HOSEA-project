setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())
library(tidyr)
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
model_path = "R_data/results/models/XGB_all_ANY.rds"

# =========================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1_test_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df
# subset to test set
df %<>% filter(id %in% test_ids)
# imputation
set.seed(0)
df = impute_srs(df, quantiles)


# =========================================================
# add correct staging
staging = read.csv("./R_data/staging.csv") %>% tibble::tibble()
colnames(staging) %<>% tolower()

patch_staging = function(df, staging){
  df %<>% left_join(staging, 
                    by=c("stagegroupclinical", "clinicalt", "clinicaln", "clinicalm"))
  df %<>% mutate(
    nccn_stage_2017=ifelse(
      (casecontrol==1) & (nccn_stage_2017==""), 
      "missing", 
      nccn_stage_2017
      ))
  return(df)
}

master %<>% patch_staging(staging)

df %<>% left_join(master %>% select(id, nccn_stage_2017), by="id")
  
# =========================================================
# function to get ROC from a df
get_roc = function(df){
  # ensure correct column ordering for xgb model
  df %<>% select(c(id, casecontrol, xgb_fit$feature_names))
  y = df$casecontrol
  # convert to xgb format
  df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                   label=df$casecontrol)
  # get predicted risk and ROC curve
  proba = predict(xgb_fit, newdata=df)
  fg = proba[y==1]; bg = proba[y==0]
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  roc$curve = roc$curve[seq(1, nrow(roc$curve), 
                            by=ceiling(nrow(df)/1000)), ]
  
  proc = pROC::roc(controls=bg, cases=fg)
  roc$ci = pROC::ci(proc, of="auc")
  roc$display.ci = paste0(
    round(roc$au, 3), " [",
    round(roc$ci[1], 3), ",",
    round(roc$ci[3], 3), "]"
  )
  roc$display = round(roc$au, 3)
  return(roc)
}




# =========================================================
# risk distribution per stage

# function to get risk histogram
get_risk_hist = function(df, breaks){
  dff = df %>% select(c(id, casecontrol, xgb_fit$feature_names))
  dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                    label=dff$casecontrol)
  proba = predict(xgb_fit, newdata=dff)
  proba_hist = proba %>% cut(breaks) %>% table()
  return((proba_hist / nrow(df)) %>% cumsum())
}

breaks = seq(0, 100000, 10) / 100000

stages = c("1", "2", "3", "3 or 4", "4", "missing", "u")

proba_hist_df = sapply(stages, function(stage){
  get_risk_hist(df %>% filter(nccn_stage_2017 == !!stage), breaks)
}) %>% data.frame()

stages_names = c("I", "II", "III", "III/IV", "IV", "missing", "unknown")
colnames(proba_hist_df) = stages_names

proba_hist_df$mid = seq(0, 100000-10, 10) + 5


proba_hist_df_long = proba_hist_df %>% pivot_longer(
  cols=stages_names
)
colnames(proba_hist_df_long) = c("mid", "stage", "freq")

filepath = paste0(dir_figures, "risk_stage.pdf")
g = ggplot(data=proba_hist_df_long %>% filter(stage %in% c("I", "II", "III", "IV")), 
           aes(x=mid, y=freq, color=stage)) + 
  geom_line() + scale_x_continuous(trans="log10") +
  xlab("Predicted risk (/100,000)") + ylab("Cumulative freq.") +
  ggtitle("Risk distribution per stage (test cases only)")
ggsave(filepath, g, width=8, height=4)




# =========================================================
# AUROCs

dff = df %>% select(c(id, casecontrol, xgb_fit$feature_names))
# convert to xgb format
xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                 label=df$casecontrol)
# get predicted risk and ROC curve
proba = predict(xgb_fit, newdata=xgb_df)


# which subsets to get ROC
stages = c("1", "2", "3", "3 or 4", "4", "missing", "u")
subsets = list(
  "all" = stages,
  "I" = c("1"),
  "II" = c("2"),
  "III" = c("3"),
  "IV" = c("4"),
  "I+" = c("1", "2", "3", "3 or 4", "4"),
  "II+" = c("2", "3", "3 or 4", "4"),
  "III+" = c("3", "3 or 4", "4"),
  "IV+" = c("4")
)

aurocs = sapply(names(subsets), function(name){
  which = subsets[[name]]
  y = df$nccn_stage_2017 %in% which
  fg = proba[y==1]; bg = proba[y==0]
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  
  proc = pROC::roc(controls=bg, cases=fg)
  roc$ci = pROC::ci(proc, of="auc")
  roc$display.ci = paste0(
    round(roc$au, 3), " [",
    round(roc$ci[1], 3), ",",
    round(roc$ci[3], 3), "]"
  )
  roc$display = round(roc$au, 3)
  
  return(c(auc=roc$auc, n=sum(y), display.ci=roc$display.ci))
})


xtable::xtable(t(aurocs) %>% data.frame() %>% select(n, display.ci), digits=3)


# EAC EGJAC 
# 8430  2965 


