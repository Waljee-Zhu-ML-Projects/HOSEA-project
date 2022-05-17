setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
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


# =========================================================
# group features
features = data.frame(name=xgb_fit$feature_names)
print(paste(xgb_fit$feature_names, collapse="   "))
features$group = c(
  'gender', 'bmi', 'weight', 
  rep("race", 4), 'agentorange', 'age', rep("smoke", 2), 
  'gerd', 'chf', 'ctd', 'dem', 'diab_c', 'hiv', 'mld', 'msld', 
  'para', 'rd', 'cd', 'copd', 'diab_nc', 'mi', 'pud', 'pvd', 
  # rep("colonoscopy", 2), rep("labs_fobt", 2), 
  rep("h2r", 5), rep("ppi", 5), 
  rep("a1c", 6), rep("bun", 6),  rep("calc", 6),  
  rep("chlor", 6),  rep("co2", 6),  rep("creat", 6), 
  rep("gluc", 6), rep("k", 6),  rep("na", 6), 
  rep("baso", 6),  rep("eos", 6),  rep("hct", 6), 
   rep("lymph", 6),  rep("mch", 6), 
  # rep("hgb", 6),  rep("lymph", 6),  rep("mch", 6), 
  rep("mchc", 6),  rep("mcv", 6),  rep("mono", 6), 
  rep("mpv", 6),  rep("neut", 6),  rep("platelet", 6), 
  # rep("rbw", 6),  
  rep("wbc", 6),  rep("crp", 6), 
  # rep("rbc", 6),  rep("rbw", 6),  rep("wbc", 6),  rep("CRP", 6), 
  rep("alkphos", 6),  rep("alt", 6),  rep("ast", 6), 
  rep("totprot", 6),  rep("hdl", 6), 
  # rep("totprot", 6),  rep("chol", 6),  rep("hdl", 6), rep("ldl", 6),  
  rep("trig", 6)
)
features$category = c( 
  rep("Demographic", 11),
  rep("Comorbidities", 16), 
  # rep("Clinical", 4), 
  rep("Medication", 2*5),
  rep("Lab", 28*6)
  # rep("Lab", 198)
)

# =========================================================
# variable importance
vi_raw = xgboost::xgb.importance(xgb_fit$feature_names, xgb_fit)
vi = vi_raw %>% left_join(features, by=c("Feature"="name")) %>%
  select(-c(Cover, Frequency))
vi_cat = vi %>% group_by(category) %>% summarize(Gain=sum(Gain))
vi_group = vi %>% group_by(group) %>% summarize(Gain=sum(Gain))
vi_group %<>% left_join(
  features %>% group_by(group, category) %>%
    summarize(group=first(group), category=first(category)), 
  by="group")
vi_group %<>% arrange(category)

vi_all = vi %>% bind_rows(vi_group) %>% bind_rows(vi_cat)
vi_all %<>% arrange(category, group, Feature)
vi_all = vi_all[, c("category", "group", "Feature", "Gain")]
vi_group = vi_group[, c("category", "group", "Gain")]
vi_cat = vi_cat[, c("category", "Gain")]

write.csv(vi_all, paste0(dir_figures, "vi.csv"))
write.csv(vi_group, paste0(dir_figures, "vi_group.
                           csv"))
write.csv(vi_cat, paste0(dir_figures, "vi_cat.csv"))



# =========================================================
# variable importance plots

g = ggplot(data=vi_cat, aes(y=reorder(category, Gain), x=Gain)) +
  geom_bar(stat="identity") +
  ggtitle("Variable importance by category") +
  xlim(0, 1) + ylab("Category")
filename = paste0(dir_figures, "pdp_new/vi_cat.pdf")
ggsave(filename, g, width=5, height=5)

g = ggplot(data=vi_group, aes(y=reorder(group, Gain), x=Gain)) +
  geom_bar(stat="identity") +
  ggtitle("Variable importance by group") +
  ylab("Feature group")
filename = paste0(dir_figures, "pdp_new/vi_group.pdf")
ggsave(filename, g, width=5, height=5)

for(cat in features$category %>% unique()){
  g = ggplot(data=vi_group %>% filter(category==!!cat), 
             aes(y=reorder(group, Gain), x=Gain)) +
    geom_bar(stat="identity") +
    ggtitle(paste0("Variable importance: ", cat))  +
    xlim(0, 1) + ylab("Group")
  filename = paste0(dir_figures, "pdp_new/vi_group_", cat, ".pdf")
  ggsave(filename, g, width=5, height=5)
}



# =========================================================
# PDPs

train = df %>% select(xgb_fit$feature_names)
set.seed(0) # sample some rows, otherwise waaaay too long
train %<>% sample_n(10000)

for(var in xgb_fit$feature_names){
  deciles = quantile(df[[var]], (1:9)/10)
  deciles = data.frame(x=deciles, xend=deciles, y=0, yend=10)
  cat(paste0("Feature: ", var, "...\n"))
  cat = features %>% filter(name==!!var) %>% use_series(category)
  group = features %>% filter(name==!!var) %>% use_series(group)
  vi_var = vi_raw %>% filter(Feature==!!var) %>% use_series(Gain)
  out = pdp::partial(xgb_fit, pred.var=var, train=train, 
                     type="classification", prob=T, which.class=1,
                     plot=F, progress="text")
  out$yhat = out$yhat * 100000
  g = ggplot(out, aes(x=get(var), y=yhat)) + 
    geom_line() + ylim(0, max(out$yhat))+ 
    ylab("Predicted risk (/100,000)") + xlab(var) +
    ggtitle(paste0("PDP: ", var, " (VI=", round(vi_var, 3), ")")) +
    geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
                 inherit.aes=F)
  if(length(unique(df[[var]]))>100) g = g + xlim(quantile(df[[var]], c(0.05, 0.95)))
  g
  filepath = paste0(dir_figures, "pdp_new/", var, ".pdf")
  ggsave(filepath, g, width=5, height=5)
  cat("...done!\n")
}



# =========================================================
# Missing values
missing_prop = df %>% is.na() %>% colMeans()
missing_prop = data.frame(missing_prop=missing_prop*100)

g = ggplot(missing_prop, aes(x=missing_prop)) + 
  geom_histogram(breaks=(0:10)/0.1) + 
  xlab("% missing by feature") + ylab("Frequency")
ggsave(paste0(dir_figures, "missing_prop_feature.png"), g, height=4, width=5)


missing_prop = df %>% is.na() %>% rowMeans()
missing_prop = data.frame(missing_prop=missing_prop*100)

g = ggplot(missing_prop, aes(x=missing_prop)) + 
  geom_histogram(breaks=c(0, 100/238, (1:30)/0.3)) + 
  xlab("% missing by patient") + ylab("Frequency")
ggsave(paste0(dir_figures, "missing_prop_patient.png"), g, height=4, width=5)

missing_prop_0 = df %>% filter(CaseControl==0) %>% is.na() %>% 
  data.frame() %>% summarize_all(mean)
missing_prop_1 = df %>% filter(CaseControl==1) %>% is.na() %>% 
  data.frame() %>% summarize_all(mean)
round(rbind(missing_prop_y, missing_prop_1) %>% t() * 100, 0)

# =========================================================
# SHAP values

dff = bind_rows(
  df %>% filter(casecontrol==1),
  df %>% filter(casecontrol==0) %>% sample_n(5000)
)
xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                      label=dff$casecontrol)
proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib=F)

#by group
groups = lapply(unique(features$group), 
                function(gr) unique(features$name[features$group==gr]))
names(groups) = unique(features$group)
shap_groups = lapply(names(groups), function(gr){
  pr = proba[, groups[[gr]]]
  if(!is.vector(pr)) pr %<>% rowSums()
  pr
})
names(shap_groups) = names(groups)
shap_groups = bind_cols(shap_groups)
shap_group_agg = abs(shap_groups) %>% colMeans()
df_shap = data.frame(sort(shap_group_agg, decreasing=F))
colnames(df_shap) = c("SHAP")
df_shap$feature = factor(rownames(df_shap), levels=rownames(df_shap))

g = ggplot(df_shap, aes(x=SHAP, y=feature)) + geom_bar(stat="identity") +
  xlab("mean|SHAP|") + ylab("") 
filepath = paste0(dir_figures, "shap_new/shap_groups.pdf")
ggsave(filepath, g, width=5, height=5)

corr_shap_groups = cor(shap_groups)
corr_shap_groups_long = 
  corr_shap_groups %>% 
  data.frame() %>% 
  mutate(feature=rownames(.)) %>% 
  tidyr::pivot_longer(rownames(.))
corr_shap_groups_long %>% 
  filter(feature!=name) %>% 
  arrange(desc(value))
filepath = paste0(dir_figures, "shap_new/shap_corr.pdf")
pdf(filepath, width=8, height=8)
corrplot::corrplot(corr_shap_groups)
dev.off()




#by cat
groups = lapply(unique(features$category), 
                function(gr) unique(features$name[features$category==gr]))
names(groups) = unique(features$category)
shap_groups = lapply(names(groups), function(gr){
  pr = proba[, groups[[gr]]]
  if(!is.vector(pr)) pr %<>% rowSums()
  pr
})
names(shap_groups) = names(groups)
shap_groups = bind_cols(shap_groups)
shap_group_agg = abs(shap_groups) %>% colMeans()
df_shap = data.frame(sort(shap_group_agg, decreasing=F))
colnames(df_shap) = c("SHAP")
df_shap$feature = factor(rownames(df_shap), levels=rownames(df_shap))

g = ggplot(df_shap, aes(x=SHAP, y=feature)) + geom_bar(stat="identity") +
  xlab("mean|SHAP|") + ylab("") 
filepath = paste0(dir_figures, "shap_new/shap_category.pdf")
ggsave(filepath, g, width=6, height=5)

shap_vals = data.frame(proba)
for(var in xgb_fit$feature_names){
  deciles = quantile(df[[var]], (1:9)/10)
  deciles = data.frame(x=deciles, xend=deciles, y=-0.1, yend=0)
  cat(paste0("Feature: ", var, "...\n"))
  cat = features %>% filter(name==!!var) %>% use_series(category)
  group = features %>% filter(name==!!var) %>% use_series(group)
  vi_var = vi_raw %>% filter(Feature==!!var) %>% use_series(Gain)
  plotdf = data.frame(
    var=dff[[var]],
    shap=exp(shap_vals[[var]])
  )
  g = ggplot(plotdf, aes(x=var, y=shap)) + geom_point(alpha=0.05) +
    geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
                 inherit.aes=F) + ylab("exp(SHAP)") + xlab(var) + 
    geom_smooth(method=ifelse(length(unique(dff[[var]]))>3, "gam", "lm"))
  g
  filepath = paste0(dir_figures, "shap_new/", var, ".pdf")
  ggsave(filepath, g, width=5, height=5)
  cat("...done!\n")
}

