setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(ggridges)
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records_old/"
data_path = paste0(dir_path, "-5--1.rds")
dir_results = "R_data/hosea-project/results/sensitivity_missing/"
dir_figures = "R_code/hosea-project/figures/sensitivity_missing/"
model_path = "R_data/results/models/finalMP_resample_nall.rds"
n_imputations = 100

# =========================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()

# =========================================================
# get prediction and summarize
# read in data
df = readRDS(data_path)
# subset to test set
df %<>% filter(ID %in% test_ids)
# ensure correct column ordering for xgb model
df %<>% select(c(ID, CaseControl, xgb_fit$feature_names))

preds = sapply(seq(n_imputations), function(i){
  cat(i)
  # imputation
  set.seed(i)
  xgb_df = impute_srs(df, quantiles)
  # convert to xgb format
  xgb_df = xgb.DMatrix(as.matrix(xgb_df[-c(1,2)]),
                   label=xgb_df$CaseControl)
  # get predicted risk and ROC curve
  proba = predict(xgb_fit, newdata=xgb_df)
  return(proba)
}, simplify=T)

pred_df = data.frame(preds)
rownames(pred_df) = df$ID
summaries = data.frame(
  mean=apply(pred_df, 1, mean),
  sd=apply(pred_df, 1, sd),
  min=apply(pred_df, 1, min),
  median=apply(pred_df, 1, median),
  max=apply(pred_df, 1, max)
)
summaries$snr = summaries$mean / summaries$sd
summaries$na_count = df %>% select(-c(ID, CaseControl)) %>% is.na() %>% rowSums()
saveRDS(summaries, paste0(dir_path, n_imputations, "imputations_summaries.rds"))

# =========================================================
# analyse
summaries = readRDS(paste0(dir_path, n_imputations, "imputations_summaries.rds"))
df$na_count = df %>% select(-c(ID, CaseControl)) %>% is.na() %>% rowSums()

# drop complete data
plot_df = summaries %>% filter(na_count>0)
plot_df$ID = as.numeric(rownames(plot_df))

# percentage of missing values
quantile(plot_df$na_count, (0:10)/10)
bins = c(0, 6, 8, 11, 12, 14, 18, 23, 31, 42, 51, 122, 197, 223, 239)
plot_df$na_count_bin = cut(plot_df$na_count, bins, include.lowest=F)
table(plot_df$na_count_bin)
table(plot_df$na_count)

g = ggplot(plot_df, aes(x=snr, y=na_count_bin, fill=na_count_bin)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75) +
  scale_x_continuous(trans="log", breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)) +
  ylab("Nb. of missing values (/239)") +
  xlab("Predicted risk SNR") +
  theme(legend.position="none")
fig_path = paste0(dir_figures, "nb_missing_values.pdf")
pdf(fig_path, width=9, height=7)
g
dev.off()
# no idea why but this breaks R
# ggsave(g, fig_path, width=8, height=6)

# demographic
demo_vars = c(
  'ageatindex','Gender','bmi','weight',
  'Asian','Black','HawaiianPacific','IndianAlaskan',
  'smoke_current','smoke_former','agentorange'
)
na_demo = df %>% select(c(ID, demo_vars)) %>% mutate_at(demo_vars, is.na)
plot_df %<>% left_join(na_demo, by="ID")
plot_df %<>% mutate(demo_any=max(across(demo_vars))==1)

dfs = list()
dfs[["any"]] = plot_df %>% filter(demo_any) %>% 
  select(c(ID, snr)) %>% mutate(variable="demo_any")
for(var in demo_vars){
  dfs[[var]] = plot_df %>% filter(plot_df %>% select(var)) %>% 
    select(c(ID, snr)) %>% mutate(variable=var)
}
demo_df = bind_rows(dfs)
demo_df$variable = factor(demo_df$variable, levels = rev(c("demo_any", demo_vars)))
g = ggplot(demo_df, aes(x=snr, y=variable, fill=variable)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75) +
  scale_x_continuous(trans="log", breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)) +
  ylab("Missing demographic variable") +
  xlab("Predicted risk SNR") +
  theme(legend.position="none")
fig_path = paste0(dir_figures, "demo_missing_values.pdf")
pdf(fig_path, width=9, height=7)
g
dev.off()


