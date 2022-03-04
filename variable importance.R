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
dir_path = "R_data/processed_records_old/"
dir_figures = "R_code/hosea-project/figures_old/"
dir_results = "R_data/results/analyses/"
model_path = "R_data/results/models_old/finalMP_resample_nall_d.rds"

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
# master = df$master
# df = df$df
repeated = table(master$ID)
repeated = names(repeated)[repeated>1]
df %<>% filter(!((ID %in% repeated) & (CaseControl==0)))
# subset to test set
df %<>% filter(ID %in% test_ids)
# imputation
set.seed(0)
df = impute_srs(df, quantiles)

xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                     label=df$CaseControl)


# =========================================================
# group features
features = data.frame(name=xgb_fit$feature_names)
print(paste(xgb_fit$feature_names, collapse="', '"))
features$group = c(
  'ageatindex', 'Gender', 'bmi', 'weight', 
  rep("Race", 4), 'agentorange', rep("Smoke", 2), 
  'GerdAtIndex', 'CHF', 'CTD', 'DEM', 'DIAB_C', 'HIV', 'MLD', 'MSLD', 
  'PARA', 'RD', 'cd', 'copd', 'diab_nc', 'mi', 'pud', 'pvd', 
  rep("colonoscopy", 2), rep("labs_fobt", 2), 
  rep("h2r", 5), rep("ppi", 5), 
  rep("A1c", 6), rep("bun", 6),  rep("calc", 6),  
  rep("chlor", 6),  rep("co2", 6),  rep("creat", 6), 
  rep("gluc", 6), rep("k", 6),  rep("na", 6), 
  rep("baso", 6),  rep("eos", 6),  rep("hct", 6), 
  rep("hgb", 6),  rep("lymph", 6),  rep("mch", 6), 
  rep("mchc", 6),  rep("mcv", 6),  rep("mono", 6), 
  rep("mpv", 6),  rep("neut", 6),  rep("platelet", 6), 
  rep("rbc", 6),  rep("rbw", 6),  rep("wbc", 6),  rep("CRP", 6), 
  rep("alkphos", 6),  rep("alt", 6),  rep("ast", 6), 
  rep("totprot", 6),  rep("chol", 6),  rep("hdl", 6), 
  rep("ldl", 6),  rep("trig", 6)
)
features$category = c( 
  rep("Demographic", 11),
  rep("Comorbidities", 16), 
  rep("Clinical", 4), 
  rep("Medication", 10),
  rep("Lab", 198)
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
write.csv(vi_group, paste0(dir_figures, "vi_group.csv"))
write.csv(vi_cat, paste0(dir_figures, "vi_cat.csv"))



# =========================================================
# variable importance plots

g = ggplot(data=vi_cat, aes(y=reorder(category, Gain), x=Gain)) +
  geom_bar(stat="identity") +
  ggtitle("Variable importance by category") +
  xlim(0, 1) + ylab("Category")
filename = paste0(dir_figures, "vi_cat.pdf")
ggsave(filename, g, width=5, height=5)

for(cat in features$category %>% unique()){
  g = ggplot(data=vi_group %>% filter(category==!!cat), 
             aes(y=reorder(group, Gain), x=Gain)) +
    geom_bar(stat="identity") +
    ggtitle(paste0("Variable importance: ", cat))  +
    xlim(0, 1) + ylab("Group")
  filename = paste0(dir_figures, "vi_group_", cat, "_.pdf")
  ggsave(filename, g, width=5, height=5)
}



# =========================================================
# PDP

vars_to_plot = unique(c(
  features %>% 
    filter(category %in% c("Demegraphic", "Comorbidities")) %>% 
    use_series(name),
  vi_raw %>% arrange(desc(Gain)) %>%
    use_series(Feature)%>% head(10)
))


train = df %>% select(xgb_fit$feature_names)
set.seed(0) # sample some rows, otherwise waaaay too long
train = train[sample.int(nrow(train), 10000), ]

for(var in vars_to_plot){
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
    geom_line() + ylim(0, 0.003*100000)+ 
    ylab("Predicted risk (/100,000)") + xlab(var) +
    ggtitle(paste0("PDP: ", var, " (VI=", round(vi_var, 3), ")")) +
    geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
                 inherit.aes=F)
  filepath = paste0(dir_figures, "pdp/", var, ".pdf")
  ggsave(filepath, g, width=3, height=4)
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
