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
model_path = "R_data/results/models/XGB_n7M_typeANY.rds"

# =========================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()

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

xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                     label=df$CaseControl)

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

rocs = list()
rocs[["All"]] = get_roc(df)

# =========================================================
# gender
male = df %>% filter(Gender==1)
rocs[["Male"]] = get_roc(male); gc(male)
female = df %>% filter(Gender==0)
rocs[["Female"]] = get_roc(female); gc(female)
df$Gender %>% mean

# =========================================================
# race
white_id = with(df, 1-pmax(Asian, Black, HawaiianPacific, IndianAlaskan))
White = df %>% filter(white_id==1)
rocs[["White"]] = get_roc(White); gc(White)
Black = df %>% filter(Black==1)
rocs[["Black"]] = get_roc(Black); gc(Black)
HawaiianPacific = df %>% filter(HawaiianPacific==1)
rocs[["HawaiianPacific"]] = get_roc(HawaiianPacific); gc(HawaiianPacific)
Asian = df %>% filter(Asian==1)
rocs[["Asian"]] = get_roc(Asian); gc(Asian)
IndianAlaskan = df %>% filter(IndianAlaskan==1)
rocs[["IndianAlaskan"]] = get_roc(IndianAlaskan); gc(IndianAlaskan)
NonWhite = df %>% filter(white_id==0)
rocs[["NonWhite"]] = get_roc(NonWhite); gc(NonWhite)

df %>% select(c(Asian, Black, HawaiianPacific, IndianAlaskan)) %>% summarise_all(mean)


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
# Gender
filepath = paste0(dir_figures, "roc_Gender.pdf")
g = ggplot(data=curves %>% filter(window %in% c("All", "Male", "Female")), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Sex") +
  ggtitle("Sensitivity analysis: Sex")
ggsave(filepath, g, width=6, height=4)


# =========================================================
# Gender
filepath = paste0(dir_figures, "roc_Race.pdf")
g = ggplot(data=curves %>% filter(window %in% 
              c("White", "Black", "HawaiianPacific", "IndianAlaskan", 
                "Asian", "All", "NonWhite")), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Race") +
  ggtitle("Sensitivity analysis: Race")
ggsave(filepath, g, width=6, height=4)


# ===
# tmp for comparison
df %<>% left_join(master, by="ID")
with(df, table(CaseControl, CancerType, Gender))


df$age_bin = cut(df$ageatindex, c(20, 30, 40, 50, 55, 60, 65, 70, 75, 80, 90))
counts = with(df, table(CaseControl, age_bin))
cc = counts %>% colSums() 
cond_prob = counts / matrix(cc, nrow=2, byrow=T, ncol=10)
round(cond_prob * 100000, 0)
