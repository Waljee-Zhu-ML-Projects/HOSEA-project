# ==============================================================================
# INSPECTING FEATURE DISTRIBUTION AND IMPUTATION
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
setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = "./R_code/hosea-project/figures/feature_comparison/"
dir_tables = "./R_code/hosea-project/tables/feature_comparison/"
imputed_data = "5-1test_mice_any.rds"
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------




# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
master = raw_df$master
raw_df = raw_df$df
raw_df %<>% filter(id %in% imputed_df$id)
diff_df = HOSEA:::mask_observed(raw_df, imputed_df)

# # for development
# raw_df %<>% sample_n(20000)
# imputed_df %<>% sample_n(20000)
# diff_df %<>% sample_n(20000)
# ------------------------------------------------------------------------------




# ==============================================================================
# COMPARE FULL DFs
comparison = HOSEA:::compare_dfs(
  df_list=list(
    observed_cases=raw_df%>%filter(casecontrol==1),
    imputed_cases=diff_df%>%filter(casecontrol==1),
    observed_controls=raw_df%>%filter(casecontrol==0),
    imputed_controls=diff_df%>%filter(casecontrol==0)
))

write.csv(comparison$fdist, paste0(dir_tables, "fdist_mice_cases.csv"))
write.csv(comparison$coherence, paste0(dir_tables, "coherence_mice_cases.csv"))
# ------------------------------------------------------------------------------

