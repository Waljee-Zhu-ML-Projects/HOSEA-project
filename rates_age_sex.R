setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())
library(pdp)
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
file_path = paste0(dir_path, "5-1_test_merged.rds")
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

proba = predict(xgb_fit, newdata=xgb_df)



# =========================================================
# put in table

dff = data.frame(
  age=df$age,
  sex=ifelse(df$gender==1, "M", "F"),
  risk=proba
)

df_group = dff %>% group_by(age, sex) %>% summarise(risk=mean(risk))

# =========================================================
# plot

filepath = paste0(dir_figures, "risk_age_sex14.pdf")
g = ggplot() + 
  geom_line(
    data=df_group,
    mapping=aes(x=age, y=risk*100000/14, color=sex)
  ) + xlab("Age") + ylab("Pred. risk (/100,000/14)") +
  scale_y_log10(breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200))
ggsave(filepath, g, width=8, height=5)
