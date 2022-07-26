setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
model_path = "R_data/results/models/XGB_all_ANY.rds"

# =========================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df
# subset to test set
df %<>% filter(id %in% test_ids)
# imputation
set.seed(0)
df = impute_srs(df, quantiles)

xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                     label=df$casecontrol)

# =========================================================
# function to get ROC from a df
get_roc = function(df){
  # ensure correct column ordering for xgb model
  df %<>% select(c(id, casecontrol, xgb_fit$feature_names))
  y = df$casecontrol
  # convert to xgb format
  df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                   label=df$casecontrol)
  # get predicted risk and ROC curve
  proba = predict(xgb_fit, newdata=df)
  fg = proba[y==1]; bg = proba[y==0]
  
  
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  roc$curve = roc$curve[seq(1, nrow(roc$curve), by=ceiling(nrow(df)/1000)), ]
  roc$curve %<>% data.frame()
  colnames(roc$curve) = c("fpr", "recall", "tr")
  
  proc = pROC::roc(controls=bg, cases=fg)
  roc$ci = pROC::ci(proc, of="auc")
  roc$display.ci = paste0(
    round(roc$au, 3), " [",
    round(roc$ci[1], 3), ",",
    round(roc$ci[3], 3), "]"
  )
  roc$display = round(roc$au, 3)
  return(roc)
}

rocs = list()
rocs[["all"]] = get_roc(df)



# =========================================================
# age
df %<>% mutate(age_group=cut(age, c(0, 45, 60, 70, 80, 100)))
age_groups = df$age_group %>% levels
for(group in age_groups){
  rocs[[group]] = get_roc(df%>%filter(age_group==group))
}

# =========================================================
# gender
male = df %>% filter(gender==1)
rocs[["male"]] = get_roc(male); gc(male)
female = df %>% filter(gender==0)
rocs[["female"]] = get_roc(female); gc(female)
df$gender %>% mean

# =========================================================
# race
white_id = with(df, 1-pmax(asian, black, hawaiianpacific, indianalaskan))
white = df %>% filter(white_id==1)
rocs[["white"]] = get_roc(white); gc(white)
black = df %>% filter(black==1)
rocs[["black"]] = get_roc(black); gc(black)
hawaiianpacific = df %>% filter(hawaiianpacific==1)
rocs[["hawaiianpacific"]] = get_roc(hawaiianpacific); gc(hawaiianpacific)
asian = df %>% filter(asian==1)
rocs[["asian"]] = get_roc(asian); gc(asian)
indianalaskan = df %>% filter(indianalaskan==1)
rocs[["indianalaskan"]] = get_roc(indianalaskan); gc(indianalaskan)
nonwhite = df %>% filter(white_id==0)
rocs[["nonwhite"]] = get_roc(nonwhite); gc(nonwhite)

df %>% select(c(asian, black, hawaiianpacific, indianalaskan)) %>% summarise_all(mean)


# =========================================================
# post processing
aucs = sapply(rocs, function(roc) roc$auc)
curves = lapply(seq_along(rocs), function(i){
  curve = data.frame(rocs[[i]]$curve)
  colnames(curve) = c("fpr", "recall", "threshold")
  nm = names(rocs)[i]
  curve$window = nm
  curve$label = paste0(nm, " (AUC: ", rocs[[i]]$display.ci, ")")
  curve
})
curves %<>% bind_rows()

# =========================================================
# age
filepath = paste0(dir_figures, "roc_age.pdf")
g = ggplot(data=curves %>% filter(window %in% c(age_groups, "all")), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Age") +
  ggtitle("Test ROC stratified by age")
ggsave(filepath, g, width=8, height=4)

# =========================================================
# gender
filepath = paste0(dir_figures, "roc_gender.pdf")
g = ggplot(data=curves %>% filter(window %in% c("all", "male", "female")), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Sex") +
  ggtitle("Test ROC stratified by sex")
ggsave(filepath, g, width=8, height=4)


# =========================================================
# race
filepath = paste0(dir_figures, "roc_Race.pdf")
g = ggplot(data=curves %>% filter(window %in% 
              c("white", "black", "hawaiianpacific", "indianalaskan", 
                "asian", "all", "nonwhite")), 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Race") +
  ggtitle("Test ROC stratified by race")
ggsave(filepath, g, width=8, height=4)
