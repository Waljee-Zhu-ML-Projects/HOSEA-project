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
  ANY=list(
    original=XGB_ANY,
    pre_egjac=readRDS(paste0(dir_models, "old/XGB_all_ANY.rds")),
    post_egjac=readRDS(paste0(dir_models, "XGB_all_ANY.rds"))
  ),
  EAC=list(
    original=XGB_EAC,
    pre_egjac=readRDS(paste0(dir_models, "old/XGB_all_EAC.rds")),
    post_egjac=readRDS(paste0(dir_models, "XGB_all_EAC.rds"))
  ),
  EGJAC=list(
    original=XGB_EGJAC,
    pre_egjac=readRDS(paste0(dir_models, "old/XGB_all_EGJAC.rds")),
    post_egjac=readRDS(paste0(dir_models, "XGB_all_EGJAC.rds"))
  )
)

# =========================================================
# read in data
dfs = list(
  pre_egjac=readRDS(paste0(dir_data, "5-1_merged_old.rds")),
  post_egjac=readRDS(paste0(dir_data, "5-1_merged.rds"))
)

# =========================================================
# function to get ROC from a df
get_roc = function(df, xgb_fit, proba=NULL){
  y = df$casecontrol
  if(is.null(proba)){
    # ensure correct column ordering for xgb model
    df %<>% select(c(id, casecontrol, xgb_fit$feature_names))
    # convert to xgb format
    df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                     label=df$casecontrol)
    # get predicted risk and ROC curve
    proba = predict(xgb_fit, newdata=df)
  }
  
  fg = proba[y==1]; bg = proba[y==0]
  
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  if(nrow(roc$curve) > 10000) {
    roc$curve = roc$curve[seq(1, nrow(roc$curve), by=ceiling(nrow(df)/1000)), ]
  }
  roc$curve %<>% data.frame()
  colnames(roc$curve) = c("fpr", "recall", "tr")
  
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
# comparison

preprocess_for_comparison = function(df, test_ids){
  dff = df %>% filter(id %in% test_ids)
  dff$white = !(dff$asian | dff$black | dff$hawaiianpacific | dff$indianalaskan)
  dff$smoke_ever = dff$smoke_current | dff$smoke_former
  complete_cases = !(
    is.na(dff$gender) |
      is.na(dff$age) |
      is.na(dff$smoke_former) |
      is.na(dff$smoke_current) |
      is.na(dff$bmi) |
      is.na(dff$white) | # will be NA iff same for other columns
      is.na(dff$gerd) |
      (is.na(dff$h2r_max) &
         is.na(dff$ppi_max))
  )
  test_complete = dff %>% subset(complete_cases)
  test_complete$casecontrol %>% table()
  
  # test_complete is in the format for our xgboost model (w/ extra variables)
  # the following keeps only the required info for kunzmann or hunt
  
  # variables
  dffk = test_complete %>% select(c(
    "id", "casecontrol", "smoke_current", "smoke_former", "gender", 
    "gerd", "white", "smoke_ever", "bmi", "age")
  )
  
  # check
  dffk %>% filter(casecontrol==1) %>% select(c(
    "id", "casecontrol", "smoke_current", "smoke_former", "gender", 
    "gerd", "white", "smoke_ever", "bmi", "age"
  )) %>% is.na() %>% colSums()
  
  # compute various bins
  dffk$k_age_bin = 
    relevel(cut(test_complete$age, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
  dffk$k_bmi_bin = cut(test_complete$bmi, breaks=c(0, 25, 30, 35, 100))
  dffk$h_age_bin = cut(test_complete$age, breaks=c(0, 50, 60, 70, 100))
  dffk$h_bmi_bin = cut(test_complete$bmi, breaks=c(0, 30, 100))
  dffk$h_smoke_any = pmax(
    test_complete$smoke_former,
    test_complete$smoke_current
  )
  dffk$k_ec = pmax(
    test_complete$h2r_max>0, 
    test_complete$ppi_max>0, 
    test_complete$gerd, na.rm=T
  ) > 0
  
  return(dffk)
}


kunzmann_score = function(dff){
  score = 0
  score = score + (dff$k_age_bin == "(55,60]") * 1.5
  score = score + (dff$k_age_bin == "(60,65]") * 2.5
  score = score + (dff$k_age_bin == "(65,100]") * 3.5
  score = score + dff$gender * 4
  score = score + (dff$k_bmi_bin == "(25,30]") * 1
  score = score + (dff$k_bmi_bin == "(30,35]") * 1.5
  score = score + (dff$k_bmi_bin == "(35,100]") * 2.5
  score = score + dff$smoke_former * 2.5
  score = score + dff$smoke_current * 3.5
  score = score + dff$k_ec * 1.5
  return(score)
}

hunt_score = function(dff){
  score = 3.6
  score = score * ifelse(dff$gender, 1.9, 1.)
  score = score * ifelse(dff$h_age_bin == "(50,60]", 2.1, 1.)
  score = score * ifelse(dff$h_age_bin == "(60,70]", 3.2, 1.)
  score = score * ifelse(dff$h_age_bin == "(70,100]", 3.1, 1.)
  score = score * ifelse(dff$h_bmi_bin == "(30,100]", 1.8, 1.)
  score = score * ifelse(dff$gerd, 3.7, 1.)
  score = score * ifelse(dff$h_smoke_any, 2.1, 1.)
  return(score/100000)
}

# =========================================================
# get ROCs
cols = c("outcome", "model", "data", "auc")
rocs = data.frame(matrix(NA, 3*0, ncol=length(cols)))
colnames(rocs) = cols

for(dfname in names(dfs)){
for(outcome in names(models)){
omodels = models[[outcome]]
for(model in names(omodels)){
  xgb_fit = omodels[[model]]$xgb_fit
  quantiles = omodels[[model]]$quantiles
  test_ids = omodels[[model]]$test_ids
  df = dfs[[dfname]]$df
  master = dfs[[dfname]]$master
  df0 = preprocess_for_comparison(df, test_ids)
  df1 = df %>% filter(id %in% (df0 %>% pull(id)))
  df1 %<>% impute_srs(quantiles)
  
  outcomes = master %>% select(id, cancertype)
  outcomes %<>% mutate(
    ANY=ifelse(cancertype=="", 0, 1),
    EAC=ifelse(cancertype=="EAC", 1, 0),
    EGJAC=ifelse(cancertype=="EGJAC", 1, 0)
  )
  df0 %<>% left_join(outcomes, by="id")
  df0 %<>% mutate(casecontrol = get(outcome))
  df1 %<>% left_join(outcomes, by="id")
  df1 %<>% mutate(casecontrol = get(outcome))
  
  roc = get_roc(df1, xgb_fit)
  row = c(outcome, model, dfname, roc$display.ci)
  print(row)
  rocs[nrow(rocs) + 1, ] = row
}
  # Kunzmann
  scores = kunzmann_score(df0)
  roc = get_roc(df0, xgb_fit, scores)
  row = c(outcome, "Kunzmann", dfname, roc$display.ci)
  print(row)
  rocs[nrow(rocs) + 1, ] = row
  # HUNT
  scores = hunt_score(df0)
  roc = get_roc(df0, xgb_fit, scores)
  row = c(outcome, "HUNT", dfname, roc$display.ci)
  print(row)
  rocs[nrow(rocs) + 1, ] = row
}}

rocs_wide = rocs %>% pivot_wider(names_from=outcome, values_from=auc)
rocs_wide %<>% select(data, model, everything())
xtable::xtable(rocs_wide)

# comparing test sets
test_ids = lapply(models$ANY, function(model){ model$test_ids})

identical = matrix(NA, 3, 3, dimnames=list(names(test_ids), names(test_ids)))
different = matrix(NA, 3, 3, dimnames=list(names(test_ids), names(test_ids)))

for(m0 in names(test_ids)){
  for(m1 in names(test_ids)){
    t0 = test_ids[[m0]]
    t1 = test_ids[[m1]]
    int = intersect(t0, t1) %>% length
    n = length(t0)
    dif = n - int
    identical[m0, m1] = int
    different[m0, m1] = dif
  }
}

# comparing test sets 
