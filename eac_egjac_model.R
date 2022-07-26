setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/utils_xgb.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_model = "R_data/results/models/"
dir_figures = "R_code/hosea-project/figures/eac_v_egjac/"
dir_results = "R_data/results/analyses/"

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df

# =========================================================
# get outcomes
outcomes = master %>% select(id, cancertype)
# subset to cases
cases_df = df %>% left_join(outcomes, by="id") %>% 
  filter(!cancertype=="") %>% 
  mutate(eac=(cancertype=="EAC")) %>% 
  select(-casecontrol, -cancertype)
cases_df %<>% mutate(casecontrol=ifelse(eac, 1, 0)) #1=EAC, 0=EGJAC
cases_df %<>% select(-eac)

# =========================================================
# fit
cases_df %<>% select(-any_of(c(
  "colonoscopy_n", "colonoscopy_maxdiff",
  "labs_fobt_n", "labs_fobt_maxdiff",
  "hgb_mean", "hgb_min", "hgb_max", "hgb_mindiff", "hgb_maxdiff", "hgb_tv",
  "rbc_mean", "rbc_min", "rbc_max", "rbc_mindiff", "rbc_maxdiff", "rbc_tv",
  "chol_mean", "chol_min", "chol_max", "chol_mindiff", "chol_maxdiff", "chol_tv"
)))

# train valid split 75/25
set.seed(0)
ids_train = cases_df %>% sample_n(round(nrow(.) * 0.75)) %>% pull(id)
train = cases_df %>% filter(id %in% ids_train)
valid = cases_df %>% filter(!(id %in% ids_train))

# replciability
n_quantiles = 10000
quantiles = compute_quantiles(train, n_quantiles)
set.seed(0)
train_ = impute_srs(train, quantiles)
valid_ = impute_srs(valid, quantiles)
dwatchlist = xgb_prep(train=train_, valid=valid_, test=valid_)
dwatchlist$test = NULL

param_xg = list(
  max_depth = 2,
  subsample = 0.5,
  eta = .05,
  objective = 'binary:logistic',
  eval_metric = 'auc',
  nthread=-1
)

xgb_fit = xgb.train(param_xg,
                    dwatchlist$train,
                    nrounds=20000,
                    dwatchlist,
                    verbose=1,print_every_n=10,
                    early_stopping_rounds=100)


xgb_df = dwatchlist$valid
shap = predict(xgb_fit, newdata=xgb_df, predcontrib=TRUE, approxcontrib=F)

# group features
features = data.frame(name=xgb_fit$feature_names)
print(paste(xgb_fit$feature_names, collapse="   "))
features$group = c(
   'gender', 'bmi_weight', 'bmi_weight', 
  rep("race", 4), 'agentorange', 'age', rep("smoke", 2), 
  'gerd', 'chf', 'ctd', 'dem', 'diab_c', 'hiv', 'mld', 'msld', 
  'para', 'rd', 'cd', 'copd', 'diab_nc', 'mi', 'pud', 'pvd', 
  rep("h2r", 5), rep("ppi", 5), 
  rep("a1c", 6), rep("bun", 6),  rep("calc", 6),  
  rep("chlor", 6),  rep("co2", 6),  rep("creat", 6), 
  rep("gluc", 6), rep("k", 6),  rep("na", 6), 
  rep("baso", 6),  rep("eos", 6),  rep("hct", 6), 
  rep("lymph", 6),  rep("mch", 6), 
  rep("mchc", 6),  rep("mcv", 6),  rep("mono", 6), 
  rep("mpv", 6),  rep("neut", 6),  rep("platelet", 6), 
  # rep("rbw", 6), 
  rep("wbc", 6),  rep("crp", 6), 
  rep("alkphos", 6),  rep("alt", 6),  rep("ast", 6), 
  rep("totprot", 6),  rep("hdl", 6), rep("ldl", 6),
  rep("trig", 6)
)
features$category = c( 
  rep("Demographic", 11),
  rep("Comorbidities", 16), 
  # rep("Clinical", 4), 
  rep("Medication", 2*5),
  rep("Lab", 29*6)
  # rep("Lab", 198)
)


groups = lapply(unique(features$group), 
                function(gr) unique(features$name[features$group==gr]))
names(groups) = unique(features$group)
shap_groups = lapply(names(groups), function(gr){
  pr = shap[, groups[[gr]]]
  if(!is.vector(pr)) pr %<>% rowSums()
  pr
})
names(shap_groups) = names(groups)
shap_groups = bind_cols(shap_groups)
shap_group_agg = abs(shap_groups) %>% colMeans()
df_shap = data.frame(sort(shap_group_agg, decreasing=F))
colnames(df_shap) = c("SHAP")
df_shap$feature = factor(rownames(df_shap), levels=rownames(df_shap))

# plot by group
g = ggplot(df_shap, aes(x=SHAP, y=feature)) + geom_bar(stat="identity") +
  xlab("mean|SHAP|") + ylab("") 
filepath = paste0(dir_figures, "vs_shap_groups.pdf")
ggsave(filepath, g, width=6, height=6)


# SHAP marginal plots
shap_vals = data.frame(shap)
for(var in xgb_fit$feature_names){
  deciles = quantile(valid[[var]], (1:9)/10, na.rm=T)
  deciles = data.frame(x=deciles, xend=deciles, y=-0.1, yend=0)
  xlims = quantile(valid[[var]], c(0.05, 0.95), na.rm=T)
  cat(paste0("Feature: ", var, "...\n"))
  plotdf = data.frame(
    Case=valid$casecontrol %>% factor(levels=c(0,1), labels=c("EGJAC", "EAC")),
    var=valid[[var]],
    shap=exp(shap_vals[[var]])
  )
  g = ggplot() + 
    geom_point(data=plotdf, mapping=aes(x=var, y=shap, color=Case), alpha=0.5) +
    geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
                 inherit.aes=F) + ylab("exp(SHAP)") + xlab(var) + 
    geom_smooth(
      data=plotdf, mapping=aes(x=var, y=shap),
      method=ifelse(length(unique(valid[[var]]))>3, "gam", "lm")) +
    xlim(xlims)
  g
  filepath = paste0(dir_figures, "vs_shap/", var, ".pdf")
  ggsave(filepath, g, width=6, height=4)
  cat("...done!\n")
}



