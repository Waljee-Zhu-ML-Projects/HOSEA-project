# ==============================================================================
# COMPARISON OF PERFORMANCE BY CANCERSTAGES
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
dir_figures = "./R_code/hosea-project/figures/cancerstage/"
dir_tables = "./R_code/hosea-project/tables/cancerstage/"
path_staging = "./R_data/staging.csv"
imputed_data = "test_mice_any.rds"
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------




# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------




# ==============================================================================
# PARAMETERS
# outcome = "ANY"
missing_which = "all"
representative = F # F: uses everything, T: downsamples males so get a more representative sample 
seed = 0
stages = c("1", "2", "3", "3 or 4", "4", "missing", "u")
subsets = list(
  "Any" = stages,
  "I" = c("1"),
  "II" = c("2"),
  "III" = c("3"),
  "IV" = c("4"),
  "I+" = c("1", "2", "3", "3 or 4", "4"),
  "II+" = c("2", "3", "3 or 4", "4"),
  "III+" = c("3", "3 or 4", "4"),
  "IV+" = c("4")
)
# ------------------------------------------------------------------------------

for(outcome in c("ANY", "EAC", "EGJAC")){


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
# compute new columns
comparison_df = imputed_df %>% compute_columns_for_comparison()
# change outcome
comparison_df %<>% patch_outcome(master=raw_df$master, outcome=outcome)
# working df
imputed_wdf = imputed_df %>% filter(id %in% ids)
imputed_wdf %<>% patch_outcome(raw_df$master, outcome=outcome)
comparison_wdf = comparison_df %>% filter(id %in% ids)
# representative sample
if(representative){
set.seed(seed)
imputed_wdf %<>% representative_sample()
comparison_wdf %<>% filter(id %in% (imputed_wdf %>% pull(id)))
}

n_cases = imputed_wdf$casecontrol %>% sum()
n_patients = imputed_wdf %>% nrow()
# ------------------------------------------------------------------------------




# ==============================================================================
# ADD STAGING
master= raw_df$master
master %<>% patch_staging(path_staging)
# ------------------------------------------------------------------------------




# ==============================================================================
# GET SCORES
scores = predict.HOSEA(imputed_wdf, imputer=NULL) %>% 
                      select(id, !!outcome) %>% rename(HOSEA=!!outcome)
scores %<>% 
  left_join(imputed_wdf %>% select(id, casecontrol), by="id") %>%
  left_join(master %>% select(id, nccn_stage_2017), by="id")
# ------------------------------------------------------------------------------




# ==============================================================================
# GET ROC

rocs = lapply(names(subsets), function(name){
  cat(name, "\n")
  which = subsets[[name]]
  # scores_ = scores %>% mutate(
  #   y=ifelse(nccn_stage_2017 %in% which, 1, 0)
  # ) %>% select(-nccn_stage_2017, -casecontrol)
  scores_ = scores %>% filter(nccn_stage_2017 %in% c(NA, which)) %>% 
    select(-nccn_stage_2017) %>% rename(y=casecontrol)
  out = roc(scores_)$HOSEA
  out$n_cases = scores_$y %>% sum()
  return(out)
})
# ------------------------------------------------------------------------------




# ==============================================================================
# WRITE TO TABLE
tab_df = data.frame(
  Stage=names(subsets),
  "Test. AUC"=sapply(rocs, function(x) x$display.ci),
  "Nb. cases"=sapply(rocs, function(x) x$n_cases)
)

colnames(tab_df) = c("Stage", "Test. AUC", "Nb. cases")

filepath = paste0(dir_tables, "auc_", outcome, "_", 
                  missing_which, 
                  ifelse(representative, "_representative", ""),
                  ".tex")



addspace = c(1, 5)

textable = xtable::xtable(tab_df, align="llrr", digits=0)
textable %<>% print(
  table.placement="ht",
  booktabs=T,
  include.rownames=F,
  add.to.row=list(pos=as.list(addspace), command=rep("\\addlinespace\n", length(addspace)))
)
cat(textable, file=filepath)
# ------------------------------------------------------------------------------


}
