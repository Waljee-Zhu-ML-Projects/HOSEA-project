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
dir_path = "R_data/processed_records/"
data_path = paste0(dir_path, "5-1.rds")
dir_results = "R_data/hosea-project/results/"
dir_figures = "R_code/hosea-project/figures/"
model_path = "R_data/results/models/XGB_nALL_typeANY.rds"
n_imputations = 210

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
master = df$master
df = df$df
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
  xgb_df = xgb.DMatrix(xgb_df %>% select(-c(ID, CaseControl)) %>% as.matrix(),
                   label=xgb_df$CaseControl)
  # get predicted risk
  proba = predict(xgb_fit, newdata=xgb_df)
  return(proba)
}, simplify=T)

pred_df = data.frame(preds)
rownames(pred_df) = df$ID
saveRDS(pred_df, paste0(dir_path, n_imputations, "imputations.rds"))

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
# Test ROC per number of imputation
n_imputations = c(1, 2, 5, 10, 20, 50, 100)
u = cumsum(n_imputations)
l = c(1, u[-length(u)]+1)

agg_preds = sapply(seq(n_imputations), function(i){
  cat(i)
  # get predicted risk
  proba = preds[, l[i]:u[i]] %>% data.frame() %>% rowMeans()
  return(proba)
}, simplify=T)

agg_preds %<>% data.frame()
rownames(agg_preds) = df$ID


y = df$CaseControl

get_roc = function(proba){
  print(length(proba))
  fg = proba[y==1]; bg = proba[y==0]
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  roc$curve = roc$curve[seq(1, nrow(roc$curve), 
                            by=ceiling(length(y)/1000)), ]
  return(roc)
}

rocs = apply(agg_preds, 2, get_roc)
names(rocs) = n_imputations

# post processing
aucs = sapply(rocs, function(roc) roc$auc)
curves = lapply(seq_along(rocs), function(i){
  curve = data.frame(rocs[[i]]$curve)
  colnames(curve) = c("fpr", "recall", "threshold")
  nm = names(rocs)[i]
  curve$window = nm
  curve$label = paste0(sprintf("%03s",nm), " (AUC=", round(aucs[nm], 3), ")")
  curve
})
curves %<>% bind_rows()

# plot
filepath = paste0(dir_figures, "roc_n_imputations.pdf")
g = ggplot(data=curves, 
           aes(x=fpr, y=recall, color=label)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Nb. imputations") +
  ggtitle("Sensitivity analysis: Multiple imputations") + 
  scale_fill_hue(breaks=n_imputations)
ggsave(filepath, g, width=6, height=5)

# =========================================================
# same but aggregated every 10
n_imputations = rep(10, 21)
u = cumsum(n_imputations)
l = c(1, u[-length(u)]+1)

agg_preds = sapply(seq(n_imputations), function(i){
  cat(i)
  # get predicted risk
  proba = preds[, l[i]:u[i]] %>% data.frame() %>% rowMeans()
  return(proba)
}, simplify=T)

agg_preds %<>% data.frame()
rownames(agg_preds) = df$ID

# Number of flips
trs = c(100, 150, 200, 250, 300)/100000

n_pred_df = sapply(trs, function(tr){
  pred = agg_preds > tr
  n_pred = pred %>% rowMeans()
  return(n_pred)
}) %>% data.frame()

colnames(n_pred_df) = trs*100000

n_pred_df %<>% tidyr::pivot_longer(everything())

filepath = paste0(dir_figures, "prediction_switch_21x10.pdf")
g = ggplot(n_pred_df, aes(x=value, fill=name)) + 
  geom_histogram(breaks=(0:21)/21, position="dodge", alpha=1.0) + 
  ggtitle("Sensitivity analysis: 10 imputations") +
  scale_y_continuous(trans="log10") +
  labs(fill="Threshold\n(/100,000)") +
  xlab("Proportion of predicted case")
ggsave(filepath, g, width=6, height=4)




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


