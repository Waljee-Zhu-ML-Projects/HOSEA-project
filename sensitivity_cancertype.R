setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(ggridges)
library(HOSEA)
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')

# =========================================================
# paths and parameters
data_path = "R_data/processed_records/5-1.rds"
outcome_path = "unzipped_data/cancertype.rds"
dir_figures = "R_code/hosea-project/figures/sensitivity_cancertype/"
dir_model = "R_data/results/models/"

# =========================================================
# read in models
model_names = c(
  ANY="finalMP_resample_nall_d.rds",
  EAC="finalMP_resample_nall_d_EAC.rds",
  EGJAC="finalMP_resample_nall_d_EGJAC.rds"
)
models = lapply(model_names, function(file){
  results = readRDS(paste0(dir_model, file))
  return(results)
})

# =========================================================
# read in data
df = readRDS(data_path)
# subset to test set
df %<>% filter(ID %in% test_ids)
# ensure correct column ordering for xgb model
df %<>% select(c(ID, CaseControl, xgb_fit$feature_names))


# =========================================================
# get outcomes (NB see also HOSEA::patch_outcome())
outcomes_raw = readRDS(outcome_path)
ids = df %>% select(ID)
outcomes = ids %>% left_join(outcomes_raw %>% select(ID, CancerType), by="ID")
outcomes %<>% mutate(CancerType=ifelse(is.na(CancerType), "None", CancerType))
outcomes %<>% mutate(
  ANY=ifelse(CancerType=="None", 0, 1),
  EAC=ifelse(CancerType=="EAC", 1, 0),
  EGJAC=ifelse(CancerType=="EGJAC", 1, 0)
)
# check
outcomes %>% select(-c(ID, CancerType)) %>% colSums()
outcome_names = c("ANY", "EAC", "EGJAC")
# merge into df
df %<>% left_join(outcomes, by="ID")

# =========================================================
# get ROC curves

rocs = list()
for(outcome in outcome_names){
  rocs[[outcome]] = list()
  # move outcome into CaseControl
  df %<>% mutate(CaseControl:=!!sym(outcome))
  for(name in names(models)){
    dff = df %>% filter(ID %in% models[[name]]$test_ids)
    set.seed(0)
    dff %<>% impute_srs(models[[name]]$quantiles)
    dff %<>% select(c(ID, CaseControl, models[[name]]$xgb_fit$feature_names))
    y = dff$CaseControl
    dff = xgb.DMatrix(as.matrix(dff[-c(1,2)]),
                     label=dff$CaseControl)
    # get predicted risk and ROC curve
    proba = predict(models[[name]]$xgb_fit, newdata=dff)
    fg = proba[y==1]; bg = proba[y==0]
    roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
    roc$curve = roc$curve[seq(1, nrow(roc$curve), by=1000), ]
    rocs[[outcome]][[name]] = roc
  }
}


# =========================================================
# postprocessing
aucs = data.frame(lapply(rocs, function(roc) sapply(roc, function(r) r$auc)))
xtable::xtable(aucs, digits=3)
# columns = outcome, rows = models
curves = list()
for(outcome in outcome_names){
  for(name in names(models)){
    curve = data.frame(rocs[[outcome]][[name]]$curve)
    colnames(curve) = c("fpr", "recall", "threshold")
    curve$outcome = outcome
    curve$model = name
    curves[[paste0(outcome, "_", name)]] = curve
  }
}
curves %<>% bind_rows()


# =========================================================
# plot
filepath = paste0(dir_figures, "roc_curves_outcome_model.pdf")
g = ggplot(data=curves, 
           aes(x=fpr, y=recall, color=outcome, linetype=model)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  ggtitle("Sensitivity analysis: Cancer type")
ggsave(filepath, g, width=9, height=7)
