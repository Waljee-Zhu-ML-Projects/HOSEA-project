setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(magrittr)
library(ggplot2)
library(HOSEA)

# ==============================================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
model_path = "R_data/results/models/XGB_nALL_typeANY.rds"

# ==============================================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()

# ==============================================================================
# histogram Nb end vs enddate
dir = "unzipped_data/"
df = load_sas(paste0(dir, "sample.sas7bdat"), "sample")
start = -5; end = -1
master = df %>% transmute(
  ID=ID,
  case=!is.na(datedx),
  start=IndexDate + start * 365 + 1,
  end=IndexDate + end * 365 + 1
)
maxenddate = 17229
freq = master %>% group_by(end, case) %>% summarize(count=n())
g = ggplot(freq, aes(x=end, y=count, color=case)) + 
  geom_bar(stat="identity") +
  scale_y_continuous(trans="log10") + 
  geom_vline(xintercept=365*seq(-10, 3) + maxenddate, 
             color="black", linetype="dashed")
ggsave(g, file=paste0(dir_figures, "histogram.pdf"),
       width=6, height=4)

# ==============================================================================
# get ROCs
rocs = list()

# =========================================================
# function to get ROC from a df
get_roc = function(df){
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
  roc$curve = roc$curve[seq(1, nrow(roc$curve), 
                            by=ceiling(nrow(df)/1000)), ]
  return(roc)
}

# comparison to similar windows
y0s = c(3, 2)
y1s = c(1, 1)

for(i in seq(2)){
  start = y0s[i]; end = y1s[i]
  window = paste0("ANY [", y0s[i], "-", y1s[i], "]")
  file_path = paste0(dir_path, start, "-", end, ".rds")
  # read in data
  df = readRDS(file_path)
  # subset to test set
  df %<>% filter(ID %in% test_ids)
  # imputation
  set.seed(0)
  df = impute_srs(df, quantiles)
  rocs[[window]] = get_roc(df)
}

# icd10 
data = readRDS(paste0(dir_path, "5-1_icd10only.rds"))
df_icd10 = data$df
meta_icd10 = data$master
meta_icd10 = meta_icd10 %>% 
  mutate(window_months=(end-start)/365) %>%
  filter(window_months>0)
df_icd10 %<>% filter(ID %in% meta_icd10$ID)

maxs = c(0.5, 1, 1.5, 2)
mins = c(0, 0.5, 1, 1.5)

for(i in seq(4)){
  ma = maxs[i]; mi = mins[i]
  window = paste0("ICD10 [", ma, "-", mi, "]")
  which = (meta_icd10$window_months>mi) & (meta_icd10$window_months<=ma)
  df = df_icd10 %>% filter(which)
  # subset to test set
  df %<>% filter(ID %in% test_ids)
  # imputation
  set.seed(0)
  df = impute_srs(df, quantiles)
  rocs[[window]] = get_roc(df)
}


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

# plot ROC curves

# 5-x curves
filepath = paste0(dir_figures, "roc_curves.pdf")
g = ggplot(data=curves %>% filter(window %in% names(rocs)), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Window") +
  ggtitle("Sensitivity analysis: ICD10 and window")
ggsave(filepath, g, width=9, height=7)