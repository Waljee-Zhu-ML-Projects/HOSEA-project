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
setwd('/nfs/turbo/umms-awaljee-secure/umms-awaljee-HOSEA/Peter files')
dir_raw_data = "./R_data/processed_records/"
dir_tables = "./R_code/hosea-project/tables/srs/features/"
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------




# ==============================================================================
# READ IN DATA
raw_df = readRDS(paste0(dir_raw_data, raw_data))
raw_df$master %<>% patch_staging()
raw_df$df %<>% left_join(raw_df$master %>% select(id, cancertype, nccn_stage_2017), by="id")
raw_df$df%<>%mutate(cancertype=ifelse(cancertype=="", "Control", cancertype))
# ------------------------------------------------------------------------------




# ==============================================================================
# COMORBIDITIES
vars = HOSEA:::charlson_vars
out = lapply(vars, function(vname){
  df = raw_df$df %>% select(cancertype, !!vname)
  tab = table(df[[vname]], df$cancertype, useNA="always")
  N = as.matrix(tab[, 1:3])
  prop = N / matrix(colSums(N), 3, 3, T)
  colnames(N) = paste0("N_", colnames(N))
  colnames(prop) = paste0("prop_", colnames(prop))
  vout = bind_cols(N, prop)
  vout$value = c("0", "1", "missing")
  vout
})
names(out) = vars
comorbidities = bind_rows(out, .id="comorbidity")
write.csv(comorbidities, paste0(dir_tables, "comorbidities.csv"))
# ------------------------------------------------------------------------------




# ==============================================================================
# MEDICATION
vars = HOSEA:::med_vars

out = lapply(vars, function(vname){
  df = raw_df$df %>% select(cancertype, starts_with(!!vname))
  vout = df %>% group_by(cancertype) %>% summarize(
    missing=sum(is.na(.data[[paste0(vname, "_max")]])),
    observed=sum(!is.na(.data[[paste0(vname, "_max")]]))
  )
  N = vout %>% tibble::column_to_rownames("cancertype") %>% t
  prop = N / matrix(colSums(N), 2, 3, T)
  colnames(N) = paste0("N_", colnames(N))
  colnames(prop) = paste0("prop_", colnames(prop))
  vout = bind_cols(N, prop)
  vout$value = c("missing", "observed")
  vout
})

names(out) = vars
medication = bind_rows(out, .id="medication")
write.csv(medication, paste0(dir_tables, "medication.csv"))
# ------------------------------------------------------------------------------



# ==============================================================================
# LAB
vars = HOSEA:::lab_vars

out = lapply(vars, function(vname){
  df = raw_df$df %>% select(cancertype, starts_with(!!vname))
  vout = df %>% group_by(cancertype) %>% summarize(
    missing=sum(is.na(.data[[paste0(vname, "_max")]])),
    observed=sum(!is.na(.data[[paste0(vname, "_max")]]))
  )
  N = vout %>% tibble::column_to_rownames("cancertype") %>% t
  prop = N / matrix(colSums(N), 2, 3, T)
  colnames(N) = paste0("N_", colnames(N))
  colnames(prop) = paste0("prop_", colnames(prop))
  vout = bind_cols(N, prop)
  vout$value = c("missing", "observed")
  vout
})

names(out) = vars
labs = bind_rows(out, .id="labs")
write.csv(labs, paste0(dir_tables, "labs.csv"))
# ------------------------------------------------------------------------------

