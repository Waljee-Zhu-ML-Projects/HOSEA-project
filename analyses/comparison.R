# ==============================================================================
# COMPARISON OF XGBOOST MODELS TO KUNZMANN, HUNT & GUIDELINES
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
dir_figures = "./R_code/hosea-project/figures/comparison/"
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
outcome = "EAC"
missing_which = "incomplete"
representative = F # F: uses everything, T: downsamples males so get a more representative sample 
seed = 0
# ------------------------------------------------------------------------------


# to get everything
for(outcome in c("ANY", "EAC", "EGJAC")){
for(missing_which in c("all", "complete", "incomplete")){
for(representative in c(T, F)){



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
comparison_df %<>% patch_outcome(raw_df$master, outcome=outcome)
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
# GET SCORES
scores = list()
scores[["HOSEA"]] = predict.HOSEA(imputed_wdf, imputer=NULL) %>% 
  select(id, !!outcome) %>% rename(HOSEA=!!outcome)
scores[["Kunzmann"]] = kunzmann_score(comparison_wdf)
scores[["HUNT"]] = hunt_score(comparison_wdf)
scores[["Guidelines"]] = guidelines(comparison_wdf)
score_df = scores %>% purrr::reduce(full_join, by="id")
score_df %<>% left_join(imputed_wdf %>% select(id, casecontrol) %>% rename(y=casecontrol), by="id")
# ------------------------------------------------------------------------------




# ==============================================================================
# GET ROC
rocs = roc(score_df)
# ------------------------------------------------------------------------------




# ==============================================================================
# PLOT
which_point = scores[["Guidelines"]] %>% select(-id) %>% colnames()
which_curve = c("HOSEA", "Kunzmann", "HUNT")

# build dfs
df_curves = sapply(which_curve, function(name){
  curve_df = data.frame(rocs[[name]]$curve[, 1:2], paste0(name, " (", rocs[[name]]$display.ci, ")"))
  colnames(curve_df) = c("fpr", "recall", "Method")
  return(curve_df)
}, simplify=F) %>% bind_rows()

df_points = sapply(which_point, function(name) rocs[[name]]$curve[2, 1:2], simplify=F) %>% bind_rows(.id="label")
df_points$Method = "Guideline"


filepath = paste0(dir_figures, outcome, "_", 
                  missing_which, 
                  ifelse(representative, "_representative", ""),
                  ".pdf")
g = ggplot(data=df_curves, aes(x=fpr, y=recall, color=Method)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  geom_point(data=df_points) + 
  ggtitle(paste0("Cancer type: ", outcome, "\n",
                 "Dataset: test, ", ifelse(complete, "complete", "imputed"), 
                              ifelse(representative, ", representative", ""), "\n",
                 "Cases: ", n_cases, "/", n_patients))
g
ggsave(filepath, g, width=8, height=6)
# ------------------------------------------------------------------------------

}}}
