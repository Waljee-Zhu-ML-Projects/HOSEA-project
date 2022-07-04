setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(patchwork)
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


# =========================================================
# group features
features = data.frame(name=xgb_fit$feature_names)
print(paste(xgb_fit$feature_names, collapse="   "))
features$group = c(
  'gender', 'bmi', 'weight', 
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


# =========================================================
# SHAP values
set.seed(0)
dff = bind_rows(
  df %>% filter(casecontrol==1),
  df %>% filter(casecontrol==0) %>% sample_n(10000)
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


shap_groups$id = dff$id

# =========================================================
# Get time series
which = lab_vars
master %<>% filter(id %in% dff$id)
dfs = list()

lab_vars = unique(features %>% filter(category=="Lab") %>% pull(group))
files_labs=c('labs_a1c.sas7bdat', 'labs_bmp.sas7bdat', 'labs_cbc.sas7bdat', 
             'labs_crp.sas7bdat', 'labs_lft.sas7bdat', 'labs_lipid.sas7bdat')

for(dir in c("unzipped_data/", "unzipped_data/cases/")){
  for(file in files_labs){
    cat(paste0("- ", dir, file, " ...\n"))
    src_df = HOSEA:::load_sas(paste0(dir, file), "labs", verbose=1)
    colnames(src_df) %<>% tolower()
    src_df %<>% left_join(master %>% select(id, start, end, casecontrol), by="id")
    src_df %<>% filter((labdate>=start)&(labdate<=end))
    src_df %<>% arrange(id, labdate)
    src_df %<>% mutate(date_to_end = end-labdate)
    df_out = src_df %>% select(-start, -end, -labdate) %>% 
      tidyr::pivot_longer(any_of(lab_vars), names_to="labtype") 
    dfs[[paste0(dir, file)]] = df_out
    cat("\n  ...done.\n")
  }
}
df_lab = bind_rows(dfs)

stat = "maxdiff"


for(labtype in lab_vars){
  lab_stat = paste0(labtype, "_", stat)
  yrange = range(df_lab %>% filter(labtype==!!labtype) %>% pull(value), na.rm=T)
  non_na = df_lab %>% filter(labtype%in%!!labtype) %>% group_by(id) %>%
    summarise(n_notna=sum(!is.na(value))) %>% filter(n_notna>1) %>% pull(id)
  df_tmp = df_lab %>% filter(labtype%in%!!labtype) %>% filter(id %in% non_na) %>% 
    select(id, date_to_end, value, casecontrol) %>% filter(date_to_end > 365)
  df_tmp %<>% left_join(dff %>% select(id, starts_with(labtype)), by="id")
  shap_tmp = shap_groups %>% filter(id %in% non_na)
  # get ids
  qs = quantile(shap_tmp%>%pull(!!labtype), c(0.50, 0.9))
  top = shap_tmp %>% filter(shap_tmp%>%pull(!!labtype)>qs[2]) %>% pull(id)
  bottom = shap_tmp %>% filter(shap_tmp%>%pull(!!labtype)<qs[1]) %>% sample_n(min(1000, nrow(.))) %>% pull(id)
  df_top = df_tmp %>% filter(id %in% top) %>% tidyr::drop_na(value)
  df_bottom = df_tmp %>% filter(id %in% bottom) %>% tidyr::drop_na(value)
  top_shap = range(shap_tmp%>%filter(id%in%top)%>%pull(!!labtype))
  bottom_shap = range(shap_tmp%>%filter(id%in%bottom)%>%pull(!!labtype))
  cases = df_tmp %>% filter(casecontrol==1) %>% pull(id)
  controls = df_tmp %>% filter(casecontrol==0) %>% sample_n(1000) %>% pull(id)
  df_cases = df_tmp %>% filter(id%in%cases) %>% tidyr::drop_na(value)
  df_controls = df_tmp %>% filter(id%in%controls) %>% tidyr::drop_na(value)
  # plot
  g_top = ggplot() + 
    geom_line(
      data=df_top, 
      mapping=aes(x=-date_to_end, y=value, group=id, color=.data[[lab_stat]]),
      alpha=0.5
    ) + ylim(yrange[1], yrange[2]) + 
    ggtitle(paste0(labtype, ": Top 10% SHAP values [", 
                   round(top_shap[1], 2), ",", round(top_shap[2], 2), "]")) +
    xlab("Days to index date") + ylab(labtype) + 
    scale_color_gradient(low="orange", high="purple")
  g_bottom = ggplot() + 
    geom_line(
      data=df_bottom, 
      mapping=aes(x=-date_to_end, y=value, group=id, color=.data[[lab_stat]]),
      alpha=0.5
    ) + ylim(yrange[1], yrange[2]) + 
    ggtitle(paste0(labtype, ": Bottom 50% SHAP values [", 
                   round(bottom_shap[1], 2), ",", round(bottom_shap[2], 2), "]")) +
    xlab("Days to index date") + ylab(labtype) + 
    scale_color_gradient(low="orange", high="purple")
  g_cases = ggplot() + 
    geom_line(
      data=df_cases, 
      mapping=aes(x=-date_to_end, y=value, group=id),
      alpha=0.1
    ) + ylim(yrange[1], yrange[2]) + 
    ggtitle(paste0(labtype, ": Cases")) +
    xlab("Days to index date") + ylab(labtype)
  g_controls = ggplot() + 
    geom_line(
      data=df_controls, 
      mapping=aes(x=-date_to_end, y=value, group=id),
      alpha=0.1
    ) + ylim(yrange[1], yrange[2]) + 
    ggtitle(paste0(labtype, ": Controls")) +
    xlab("Days to index date") + ylab(labtype)
  g = g_top + g_cases + g_bottom + g_controls + plot_layout(ncol=2)
  ggsave(paste0(dir_figures, "lab/", labtype, ".pdf"), width=8, height=8)
}

# compare with n labs
n_labs = df_lab %>% group_by(id) %>% summarize(n=sum(!is.na(value)))

df = shap_groups %>% select(id, !!labtype) %>% left_join(n_labs, by="id")

ggplot(df, aes(x=n+1, y=mch)) + 
  geom_point() + 
  geom_smooth() + 
  scale_x_log10() + 
  xlab("Nb. lab results") + ylab("SHAP") +
  ggtitle(labtype)


# test
dir = "unzipped_data/cases/"
file = 'labs_cbc.sas7bdat'
type="hct"
verbose=2
i = 20000699
src_df = HOSEA:::load_sas(paste0(dir, file), "labs", verbose=verbose-1)
colnames(src_df) %<>% tolower()
# restrict to prediction window
src_df %<>% left_join(master %>% select(id, start, end), by="id")
src_df %<>% filter((labdate>=start)&(labdate<=end))
# ensure ordered
src_df %<>% arrange(id, labdate)
src_df %<>% filter(id==i)
# compute vars
tmp = src_df %>% select(id, labdate, !!type)
tmp %<>% tidyr::drop_na(!!type)
colnames(tmp) = c("id", "labdate", "var")
# compute lag variables
tmp %<>% mutate(
  labdate_lag = lag(labdate),
  var_lag = lag(var)
)
# compute diff and slope
tmp %<>% mutate(
  dlabdate = pmax(1, labdate - labdate_lag),
  dvar = var - var_lag
)
tmp %<>% mutate(
  svar = dvar / dlabdate
)
tmp %<>% select(id, labdate, labdate_lag, dlabdate, var, var_lag, dvar, svar)
# compute summaries
out = tmp %>% group_by(id) %>%
  summarize(
    mean = HOSEA:::safe_mean(var),
    max = HOSEA:::safe_max(var),
    min = HOSEA:::safe_min(var),
    maxdiff = HOSEA:::safe_max(svar),
    mindiff = HOSEA:::safe_min(svar),
    tv = HOSEA:::safe_mean(abs(svar)),
  )

dff %>% select(id, casecontrol, hct_maxdiff) %>% 
  filter(casecontrol==1) %>%
  arrange(hct_maxdiff, decreasing=T)
