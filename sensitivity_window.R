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
model_path = "R_data/results/models/XGB_nALL_typeANY.rds"
model_path = "R_data/results/models/XGB_n7M_typeANY.rds"

# =========================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()

# =========================================================
# get results for all windows
y0s = c(5, 5, 5, 5, 4, 3, 2, 4)
y1s = c(4, 3, 2, 1, 1, 1, 1, 2)
rocs = list()

for(i in seq(8)){
  start = y0s[i]; end = y1s[i]
  cat(paste(i, start, end, "\n"))
  window = paste0("[", y0s[i], "-", y1s[i], "]")
  file_path = paste0(dir_path, start, "-", end, ".rds")
  # read in data
  df = readRDS(file_path)$df
  # subset to test set
  df %<>% filter(ID %in% test_ids)
  # imputation
  set.seed(0)
  df = impute_srs(df, quantiles)
  # ensure correct column ordering for xgb model
  df %<>% select(c(ID, CaseControl, xgb_fit$feature_names))
  y = df$CaseControl
  # convert to xgb format
  df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                        label=df$CaseControl)
  # get predicted risk and ROC curve
  proba = predict(xgb_fit, newdata=df)
  fg = proba[y==1]; bg = proba[y==0]
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  roc$curve = roc$curve[seq(1, nrow(roc$curve), by=1000), ]
  rocs[[window]] = roc
}

# =========================================================
# post processing
aucs = sapply(rocs, function(roc) roc$auc)
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
# plot ROC curves

# # 5-x curves
# filepath = paste0(dir_figures, "roc_window_5-x.pdf")
# g = ggplot(data=curves %>% filter(window %in% names(rocs)[1:4]), 
#            aes(x=fpr, y=recall, color=label)) + 
#   geom_line() +
#   theme(aspect.ratio=1) +
#   xlab("1 - Specificity") + ylab("Sensitivity") + 
#   geom_abline(intercept=0, slope=1, linetype="dotted") +
#   labs(color="Window") +
#   ggtitle("Sensitivity analysis: [5-x] prediction window")
# ggsave(filepath, g, width=6, height=5)

# this one is better to show time to prediction
# 2 curves
filepath = paste0(dir_figures, "roc_window_2.pdf")
g = ggplot(data=curves %>% filter(window %in% names(rocs)[c(2, 8, 6)]), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted")+
  labs(color="Window") +
  ggtitle("Sensitivity analysis: [x to x-2] prediction window")
ggsave(filepath, g, width=6, height=5)

# x-1 curves
filepath = paste0(dir_figures, "roc_window_x-1.pdf")
g = ggplot(data=curves %>% filter(window %in% names(rocs)[4:7]), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Window") +
  ggtitle("Sensitivity analysis: [x-1] prediction window")
ggsave(filepath, g, width=6, height=5)

