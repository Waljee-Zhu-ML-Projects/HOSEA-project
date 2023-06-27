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
imputation = "srs"
# ------------------------------------------------------------------------------




# ==============================================================================
# PATHS
# setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/comparison_age/")
dir_rocs = dir_figures
imputed_data = paste0("5-1test_", imputation, "_any.rds")
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------



# ==============================================================================
# MODELS
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
# ------------------------------------------------------------------------------



# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------




# ==============================================================================
# PARAMETERS
outcome = "ANY"
missing_which = "complete"
representative = F # F: uses everything, T: downsamples males so get a more representative sample
seed = 0
below50 = F
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
# subset to age
if(below50){
  imputed_wdf %<>% filter(age < 50)
  comparison_wdf %<>% filter(age < 50)
}else{
  imputed_wdf %<>% filter(age >= 50)
  comparison_wdf %<>% filter(age >= 50)
}

n_cases = imputed_wdf$casecontrol %>% sum()
n_patients = imputed_wdf %>% nrow()
# ------------------------------------------------------------------------------




# ==============================================================================
# GET SCORES
scores = list()
scores[["HOSEA"]] = predict.HOSEA(imputed_wdf, models, imputer=NULL) %>% 
  select(id, !!outcome) %>% rename("K-ECAN"=!!outcome)
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
which_curve = c("K-ECAN", "Kunzmann", "HUNT")

# build dfs
df_curves = sapply(which_curve, function(name){
  curve_df = data.frame(rocs[[name]]$curve[, 1:2], paste0(name, " (", rocs[[name]]$display.ci, ")"))
  colnames(curve_df) = c("fpr", "recall", "Method")
  return(curve_df)
}, simplify=F) %>% bind_rows()

df_points = sapply(which_point, function(name) rocs[[name]]$curve[2, 1:2], simplify=F) %>% bind_rows(.id="label")
df_points$Method = "Guideline"
df_points %<>% arrange(fpr)

xorder = order(df_points$fpr)
xs = df_points$fpr[xorder]
ys = df_points$recall[xorder]

df_points$xlab = df_points$fpr
df_points$ylab = df_points$recall + 0.1*(-1)^(seq_along(which_point)) +
  0.05 * (-1)^ceiling(seq_along(which_point)/2)
df_points$ylab = ifelse(df_points$ylab > 1, 1-(df_points$ylab-1), df_points$ylab)


filepath = paste0(dir_figures, outcome, "_", 
                  missing_which, 
                  ifelse(below50, "_lt50", "ge50"),
                  ifelse(representative, "_representative", ""),
                  ".pdf")
g = ggplot(data=df_curves, aes(x=fpr, y=recall, color=Method)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  geom_point(data=df_points) + 
  ggtitle(paste0("Cancer type: ", outcome, "\n",
                 "Dataset: test, ", missing_which, 
                 ifelse(below50, ", <50", ">=50"),
                              ifelse(representative, ", representative", ""), "\n",
                 "Cases: ", n_cases, "/", n_patients))
g = g +
  # geom_segment(data=df_points, aes(x=fpr, xend=xlab, y=recall, yend=ylab)) + 
  ggrepel::geom_label_repel(
    data=df_points, 
    mapping=aes(label=label, x=fpr, y=recall), 
    size=3, 
    show.legend=FALSE,
    max.overlaps=20,
    direction="both",
    min.segment.length=0,
    nudge_x=0.1, 
    nudge_y=-0.1
    )
g
ggsave(filepath, g, width=8, height=6, bg="white")
ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=8, height=6, bg="white")
# ------------------------------------------------------------------------------





