# ==============================================================================
# ICD10 COHORT ANALYSIS
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
dir_icd10_processed = "./R_data/processed_records/icd10cohort/" 
dir_icd10_imputed = "./R_data/imputed_records/icd10cohort/" 
dir_og_processed = "./R_data/processed_records/"
dir_og_imputed = "./R_data/imputed_records/"
dir_figures = "./R_code/hosea-project/figures/srs/icd10/"
dir_tables = "./R_code/hosea-project/tables/srs/icd10/"
# ------------------------------------------------------------------------------







# ==============================================================================
# READ IN DATA
y0=3; y1=1
icd10_raw = readRDS(paste0(dir_icd10_processed, y0, "-", y1, "filtered.rds"))
icd10_imputed = readRDS(paste0(dir_icd10_imputed, y0, "-", y1, ".rds"))
icd10_master = readRDS(paste0(dir_icd10_processed, y0, "-", y1, ".rds"))$master
test_imputed = readRDS(paste0(dir_og_imputed, y0, "-", y1, "test_srs_any.rds"))
original = readRDS(paste0(dir_og_processed, y0, "-", y1, "_merged.rds"))$df
master = readRDS(paste0(dir_og_processed, y0, "-", y1, "_merged.rds"))$master
test_raw = original %>% filter(id %in% test_imputed$id)
# ------------------------------------------------------------------------------






# ==============================================================================
# COMPARE FEATURE DISTRIBUTIONS
comparison = HOSEA:::compare_dfs(
  df_list=list(
    new_cohort_raw=icd10_raw,
    new_cohort_imputed=icd10_imputed,
    test_cohort_raw=test_raw,
    test_cohort_imputed=test_imputed
  ))

write.csv(comparison$fdist, paste0(dir_tables, "fdist.csv"))
write.csv(comparison$coherence, paste0(dir_tables, "coherence.csv"))
# ------------------------------------------------------------------------------






# ==============================================================================
# GET SCORES
scores_icd10 = predict.HOSEA(icd10_imputed, imputer=NULL)
scores_test = predict.HOSEA(test_imputed, imputer=NULL)
# merge in outcome (ANY)
scores_icd10 %<>% 
  left_join(icd10_imputed %>% select(casecontrol, id), by="id") %>%
  rename(y=casecontrol)
scores_test %<>% 
  left_join(test_imputed %>% select(casecontrol, id), by="id") %>%
  rename(y=casecontrol)
# ------------------------------------------------------------------------------


for(mname in c("ANY", "EAC", "EGJAC")){
  




# ==============================================================================
# GET ROC
score_icd10 = scores_icd10 %>%
  select(id, !!mname, y) %>%
  rename(casecontrol=y, score=!!mname) %>%
  patch_outcome(master=icd10_master, outcome=mname) %>%
  rename(y=casecontrol)
score_test = scores_test %>%
  select(id, !!mname, y) %>%
  rename(casecontrol=y, score=!!mname) %>%
  patch_outcome(master=master, outcome=mname) %>%
  rename(y=casecontrol)

# # representative sample
# set.seed(0)
# # merge in gender
# score_icd10 %<>% left_join(icd10_imputed %>% select(id, gender), by="id")
# score_test %<>% left_join(test_imputed %>% select(id, gender), by="id")
# # subsample males
# nf_icd10 = (1-score_icd10 %>% pull(gender) )%>% sum()
# score_icd10 = bind_rows(
#   score_icd10 %>% filter(gender==0),
#   score_icd10 %>% filter(gender==1) %>% sample_n(nf_icd10)
# ) %>% select(-gender)
# 
# nf_test = (1-score_test %>% pull(gender) )%>% sum()
# score_test = bind_rows(
#   score_test %>% filter(gender==0),
#   score_test %>% filter(gender==1) %>% sample_n(nf_test)
# ) %>% select(-gender)


roc_icd10 = roc(score_icd10)
roc_test = roc(score_test)  

cat(mname, "\n")
cat("ICD 10\n")
print(score_icd10$y %>% table())
cat("Test\n")
print(score_test$y %>% table())
# ------------------------------------------------------------------------------






# ==============================================================================
# PLOT ROC
  df_curves = bind_rows(
    data.frame(
      roc_icd10[["score"]]$curve[, 1:2], 
      Cohort=paste0("ICD 10 (", roc_icd10[["score"]]$display.ci, ")")
    ),
    data.frame(
      roc_test[["score"]]$curve[, 1:2], 
      Cohort=paste0("Test (", roc_test[["score"]]$display.ci, ")")
    )
  )
# filepath = paste0(dir_figures, "roc_", mname, "_representative.pdf")
filepath = paste0(dir_figures, "roc_", mname, ".pdf")
  gtitle = paste0("Cancer type: ", mname, "\n")
  g = ggplot(data=df_curves, aes(x=fpr, y=recall, color=Cohort)) + geom_line() +
    theme(aspect.ratio=1) +
    xlab("1 - Specificity") + ylab("Sensitivity") + 
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    ggtitle(gtitle)
  g
  ggsave(filepath, g, width=8, height=6, bg="white")
  ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=8, height=6, bg="white")
# ------------------------------------------------------------------------------


}




  
# ==============================================================================
# PLOT RISK CURVES
for(mname in c("ANY", "EAC", "EGJAC")){
  df = bind_rows(
    scores_icd10 %>% select(!!mname) %>% mutate(Cohort="ICD 10"),
    scores_test %>% select(!!mname) %>% mutate(Cohort="Test")
  ) %>% rename(score=!!mname)
  filepath = paste0(dir_figures, "score_dist_", mname, ".pdf")
  gtitle = paste0("Cancer type: ", mname, "\n")
  g = ggplot(data=df, aes(x=score*100000, color=Cohort)) + 
    geom_density() + 
    xlab("Predicted risk (/100,000)") + ylab("Density") + 
    ggtitle(gtitle) + scale_x_log10()
  g
  ggsave(filepath, g, width=8, height=4, bg="white")
  ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=8, height=4, bg="white")
}
# ------------------------------------------------------------------------------