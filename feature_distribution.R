setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
dir_model = "R_data/results/models/"

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df
outcomes = master %>% select(id, cancertype)
df %<>% left_join(outcomes, by="id")
df %<>% mutate(cancertype=ifelse(cancertype=="", "Control", cancertype))


# test cases
model = readRDS(paste0(dir_model, "XGB_all_ANY.rds"))
test_ids = model$test_ids
quantiles = model$quantiles
xgb_fit = model$xgb_fit

df %<>% filter(id %in% test_ids)

set.seed(0)
test_imputed = impute_srs(df, quantiles)

# get complete cases
complete_cases = !(
  is.na(df$gender) |
    is.na(df$age) |
    is.na(df$smoke_former) |
    is.na(df$smoke_current) |
    is.na(df$bmi) |
    is.na(df$black) | # will be NA iff same for other columns
    is.na(df$gerd) |
    (is.na(df$h2r_max) &
       is.na(df$ppi_max))
)

# get new datasets
test_complete = df %>% filter(complete_cases)

# =========================================================
# compute new variables 

compute_new_vars = function(df){
  df$white = !(df$asian | df$black | df$hawaiianpacific | df$indianalaskan)
  df$smoke_ever = df$smoke_current | df$smoke_former
  
  # compute various bins
  df$k_age_bin = 
    relevel(cut(df$age, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
  df$k_bmi_bin = cut(df$bmi, breaks=c(0, 25, 30, 35, 100))
  df$h_age_bin = cut(df$age, breaks=c(0, 50, 60, 70, 100))
  df$h_bmi_bin = cut(df$bmi, breaks=c(0, 30, 100))
  df$h_smoke_any = pmax(
    df$smoke_former,
    df$smoke_current
  )
  df$k_ec = pmax(
    df$h2r_max>0, 
    df$ppi_max>0, 
    df$gerd, na.rm=T
  ) > 0
  return(df)
}



# =========================================================
# create the representative sample
representative_sample = function(cdf){
  ids_cases_female = cdf%>%filter(gender==0,casecontrol==1)%>%pull(id)
  ids_cases_male = cdf%>%filter(gender==1,casecontrol==1)%>%pull(id)
  ids_controls_female = cdf%>%filter(gender==0,casecontrol==0)%>%pull(id)
  ids_controls_male = cdf%>%filter(gender==1,casecontrol==0)%>%pull(id)
  
  n_cases_female = length(ids_cases_female)
  n_controls_female = length(ids_controls_female)
  n_cases_male = round(8.33*n_cases_female)
  n_controls_male = n_controls_female + n_cases_female - n_cases_male
  
  set.seed(0)
  ids_cases_male = sample(ids_cases_male, n_cases_male, replace=F)
  ids_controls_male = sample(ids_controls_male, n_controls_male, replace=F)
  
  ids = c(
    ids_cases_female,
    ids_cases_male,
    ids_controls_female,
    ids_controls_male
  )
  
  dff = cdf %>% filter(id %in% ids)
  return(dff)
}

# =========================================================
set.seed(0)
test_complete %<>% representative_sample()
test_representative = test_imputed %>% representative_sample()

test_complete %<>% compute_new_vars()
test_representative %<>% compute_new_vars()
test_imputed %<>% compute_new_vars()


# ==============================================================================
# compute summaries

summarizes_features = function(test_df){
  out = test_df %>% group_by(cancertype) %>% summarise(
    n=n(),
    smoke_current=mean(smoke_current, na.rm=T),
    smoke_former=mean(smoke_former, na.rm=T),
    smoke_ever=mean(smoke_ever, na.rm=T),
    male=mean(gender, na.rm=T),
    gerd=mean(gerd, na.rm=T),
    gerd_or_med=mean(k_ec, na.rm=T),
    white=mean(white, na.rm=T),
    black=mean(black, na.rm=T),
    asian=mean(asian, na.rm=T),
    hawaiianpacific=mean(hawaiianpacific, na.rm=T),
    indianalaskan=mean(indianalaskan, na.rm=T),
    bmi=mean(bmi, na.rm=T),
    bmi_0_25=mean(k_bmi_bin=="(0,25]", na.rm=T),
    bmi_25_30=mean(k_bmi_bin=="(25,30]", na.rm=T),
    bmi_30_35=mean(k_bmi_bin=="(30,35]", na.rm=T),
    bmi_35_100=mean(k_bmi_bin=="(35,100]", na.rm=T),
    bmi_0_30=mean(h_bmi_bin=="(0,30]", na.rm=T),
    bmi_30_100=mean(h_bmi_bin=="(30,100]", na.rm=T),
    age=mean(age, na.rm=T),
    age_0_25=mean(k_age_bin=="(0,50]", na.rm=T),
    age_25_30=mean(k_age_bin=="(50,55]", na.rm=T),
    age_30_35=mean(k_age_bin=="(55,60]", na.rm=T),
    age_35_100=mean(k_age_bin=="(60,65]", na.rm=T),
    age_35_100=mean(k_age_bin=="(65,100]", na.rm=T),
    age_0_50=mean(h_age_bin=="(0,50]", na.rm=T),
    age_50_60=mean(h_age_bin=="(50,60]", na.rm=T),
    age_60_70=mean(h_age_bin=="(60,70]", na.rm=T),
    age_70_100=mean(h_age_bin=="(70,100]", na.rm=T)
  ) %>% tibble::column_to_rownames("cancertype") %>% t()
  return(out)
}

summaries = lapply(list(imputed=test_imputed, representative=test_representative, complete=test_complete),
                   summarizes_features)

lapply(summaries, function(s) xtable::xtable(s %>% round(2)))


df$k_age_bin = 
  relevel(cut(df$age, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
df$k_bmi_bin = cut(df$bmi, breaks=c(0, 25, 30, 35, 100))
df$h_age_bin = cut(df$age, breaks=c(0, 50, 60, 70, 100))
df$h_bmi_bin = cut(df$bmi, breaks=c(0, 30, 100))


# ==============================================================================
# ROCs M/F

add_ci = function(roc){
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

get_roc = function(df){
  y = df$casecontrol
  xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                       label=df$casecontrol)
  proba = predict(xgb_fit, newdata=xgb_df)
  fg = proba[y==1]; bg = proba[y==0]
  x_roc = PRROC::roc.curve(fg, bg ,curve=TRUE) %>% add_ci()
  return(x_roc)
}

rocs = lapply(
  list(
    all=test_imputed,
    males=test_imputed %>% filter(gender==1),
    females=test_imputed %>% filter(gender==0)
  ),
  get_roc
)

sapply(rocs, function(r) r$display)

