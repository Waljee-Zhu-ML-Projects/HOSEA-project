setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
# model_path = "R_data/results/models/XGB_nALL_typeANY.rds"
model_path = "R_data/results/models/XGB_n1M_typeANY_all.rds"

# =========================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()
rm(xgb_fit)

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1.rds")
df = readRDS(file_path)
master = df$master
df = df$df
# subset to test set
df %<>% filter(ID %in% test_ids)
# imputation
set.seed(0)
df = impute_srs(df, quantiles)

# =========================================================
# subsets
clinical = c("colonoscopy_n", "colonoscopy_maxdiff",
             "labs_fobt_n", "labs_fobt_maxdiff")
hgb = c("hgb_mean", "hgb_min", "hgb_max", "hgb_mindiff", "hgb_maxdiff", "hgb_tv")
hct = c("hct_mean", "hct_min", "hct_max", "hct_mindiff", "hct_maxdiff", "hct_tv")
rbc = c("rbc_mean", "rbc_min", "rbc_max", "rbc_mindiff", "rbc_maxdiff", "rbc_tv")
chol = c("chol_mean", "chol_min", "chol_max", "chol_mindiff", "chol_maxdiff", "chol_tv")
ldl = c("ldl_mean", "ldl_min", "ldl_max", "ldl_mindiff", "ldl_maxdiff", "ldl_tv")
hdl = c("hdl_mean", "hdl_min", "hdl_max", "hdl_mindiff", "hdl_maxdiff", "hdl_tv")

subsets = list(
  "all" = c(),
  "proposed" = c(clinical, hgb, rbc, chol),
  "no_clinical" = clinical,
  "no_rbc_hgb" = c(hgb, rbc),
  "no_rbc_hct" = c(hct, rbc),
  "no_hct_hgb" = c(hgb, hct),
  "no_chol" = c(chol),
  "no_ldl" = c(ldl),
  "no_hdl" = c(hdl),
  "no_hdl_ldl" = c(hdl, ldl)
)


# =========================================================
# function to get ROC from a df
get_roc = function(df, xgb_fit){
  y = df$CaseControl
  # convert to xgb format
  dff = df %>% select(xgb_fit$feature_names)
  xgb_df = xgb.DMatrix(as.matrix(dff),
                   label=y)
  # get predicted risk and ROC curve
  proba = predict(xgb_fit, newdata=xgb_df)
  print(head(proba))
  fg = proba[y==1]; bg = proba[y==0]
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  roc$curve = roc$curve[seq(1, nrow(roc$curve), 
                            by=ceiling(nrow(df)/1000)), ]
  return(roc)
}

rocs = lapply(names(subsets), function(mname){
  filepath = paste0("R_data/results/models/XGB_n1M_typeANY_", mname, ".rds")
  print(filepath)
  xgb_fit = readRDS(filepath)$xgb_fit
  print(paste(xgb_fit$best_score, xgb_fit$best_ntreelimit))
  roc = get_roc(df, xgb_fit)
  print(paste(mname, xgb_fit$nfeatures, roc$auc))
  return(roc)
})
names(rocs) = names(subsets)


# =========================================================
# post processing
test_aucs = sapply(rocs, function(roc) roc$auc)
valid_auc = sapply(names(subsets), function(mname){
  filepath = paste0("R_data/results/models/XGB_n1M_typeANY_", mname, ".rds")
  xgb_fit = readRDS(filepath)$xgb_fit
  return(xgb_fit$best_score)
})
tab = rbind(valid_auc, test_aucs) %>% t()
xtable::xtable(tab, digits=3)
curves = lapply(seq_along(rocs), function(i){
  curve = data.frame(rocs[[i]]$curve)
  colnames(curve) = c("fpr", "recall", "threshold")
  nm = names(rocs)[i]
  curve$window = nm
  curve$label = paste0(nm, " (AUC=", round(aucs[nm], 3), ")")
  curve
})
curves %<>% bind_rows()

# =========================================================
# Gender
filepath = paste0(dir_figures, "roc_subsets.pdf")
g = ggplot(data=curves, aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Sex") +
  ggtitle("Sensitivity analysis: Variables")
ggsave(filepath, g, width=6, height=4)