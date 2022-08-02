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
dir_figures = "R_code/hosea-project/figures/diagnosing/"
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

staging = read.csv("./R_data/staging.csv") %>% tibble::tibble()
colnames(staging) %<>% tolower()

patch_staging = function(df, staging){
  df %<>% left_join(staging, 
                    by=c("stagegroupclinical", "clinicalt", "clinicaln", "clinicalm"))
  df %<>% mutate(
    nccn_stage_2017=ifelse(
      (casecontrol==1) & (nccn_stage_2017==""), 
      "missing", 
      nccn_stage_2017
    ))
  return(df)
}

master %<>% patch_staging(staging)

# =========================================================
# wbc
dff = bind_rows(
  df %>% filter(casecontrol==0) %>% sample_n(20000),
  df %>% filter(casecontrol==1)
)

xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                      label=dff$casecontrol)
proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib=F)
proba %<>% data.frame()
proba %<>% select(starts_with("wbc"))
colnames(proba) = paste0(colnames(proba), "_shap")

df_wbc = bind_cols(
  dff %>% select(id, casecontrol, starts_with("wbc")) %>% 
    left_join(master %>% select(id, nccn_stage_2017), by="id"),
  proba
  )


# =========================================================
# plot
summaries = c("mean", "min", "max", "tv", "maxdiff", "mindiff")

summary = "max"

df_wbc %<>% mutate(staging=ifelse(nccn_stage_2017=="", "Control", 
                                  ifelse(nccn_stage_2017 %in% c("0", "1"), "0/I", 
                                  ifelse(nccn_stage_2017 %in% c("2", "3", "4", "3 or 4"), "II+", "other"))))

for(summary in summaries){
  g = GGally::ggpairs(
    data=df_wbc %>% filter(!(staging=="other")),
    mapping=aes(color=staging, alpha=0.2),
    columns=paste0("wbc_", summary, c("", "_shap"))
  )
  filepath = paste0(dir_figures, "wbc_", summary, ".pdf")
  ggsave(filepath, g, width=8, height=8)
}
