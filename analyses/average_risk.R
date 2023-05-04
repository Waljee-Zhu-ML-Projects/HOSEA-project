# ==============================================================================
# CALILBRATION CURVE AND THRESHOLD PERFORMANCE TABLE
# Author: Simon Fontaine (simfont@umich.edu)
# ------------------------------------------------------------------------------




# ==============================================================================
# REQUIRED PACKAGES
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(HOSEA)
theme_set(theme_minimal())
# ------------------------------------------------------------------------------




# ==============================================================================
# PATHS
imputation = "srs"
setwd('/nfs/turbo/umms-awaljee-secure/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/average_risk/")
dir_tables = paste0("./R_code/hosea-project/tables/", imputation, "/average_risk/")
imputed_data = paste0("5-1test_", imputation, "_any.rds")
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------



# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------




# ==============================================================================
# PARAMETERS
seed = 0
missing_which = "all"
representative = F
# ------------------------------------------------------------------------------



# ==============================================================================
# PREPARE DATA
# get complete cases for Kunzmann, HUNT and guidelines
if(missing_which == "complete"){
  ids = complete_for_comparison(raw_df$df)
}
if(missing_which == "all"){
  ids = imputed_df %>% pull(id)
}
if(missing_which == "incomplete"){
  complete_ids = complete_for_comparison(raw_df$df)
  ids = imputed_df %>% filter(!(id %in% complete_ids)) %>% pull(id)
}
# working df
imputed_wdf = imputed_df %>% filter(id %in% ids)
# imputed_wdf %<>% patch_outcome(master=raw_df$master, outcome=outcome)
# representative sample
if(representative){
  set.seed(seed)
  imputed_wdf %<>% representative_sample()
}
# ------------------------------------------------------------------------------




# ==============================================================================
# GET SCORES
models = load_models(
  files_meta=list(
    ANY=paste0("xgb_", imputation, "_any.meta"), 
    EAC=paste0("xgb_", imputation, "_eac.meta"), 
    EGJAC=paste0("xgb_", imputation, "_egjac.meta")
  ),
  files_models=list(
    ANY=paste0("xgb_", imputation, "_any.model"), 
    EAC=paste0("xgb_", imputation, "_eac.model"), 
    EGJAC=paste0("xgb_", imputation, "_egjac.model")
  )
)
proba = predict.HOSEA(imputed_wdf, imputer=NULL, models=models)
# ------------------------------------------------------------------------------






# ==============================================================================
# AVERAGE PREDICTED RISK PER AGE GROUP & SEX
df = proba %>% left_join(imputed_df %>% select(id, age, gender), by="id")
df %<>% mutate(age_bin=cut(age, c(0, 30, 40, 50, 60, 70, 80, 100)))
avg_risk = df %>% group_by(age_bin, gender) %>% summarize(
  ANY=mean(ANY)*100000,
  EAC=mean(EAC)*100000,
  EGJAC=mean(EGJAC)*100000
)
write.csv(avg_risk, paste0(dir_tables, "avg_risk_by_age_bin_sex.csv"))
# ------------------------------------------------------------------------------




# ==============================================================================
# POPULATION ATTRIBUTABLE RISK PER AGE GROUP & SEX
df = proba %>% left_join(imputed_df %>% select(id, age, gender), by="id")
df %<>% mutate(age_bin=cut(age, c(0, 30, 40, 50, 60, 70, 80, 100)))
par_denum = df %>% group_by(age_bin, gender) %>% summarize(
  ANY=1/mean(1/ANY)*100000,
  EAC=1/mean(1/EAC)*100000,
  EGJAC=1/mean(1/EGJAC)*100000
)
write.csv(par_denum, paste0(dir_tables, "par_denum_by_age_bin_sex.csv"))
# ------------------------------------------------------------------------------
      
      