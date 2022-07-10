setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(tidyr)
library(xgboost)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())
library(pdp)
library(HOSEA)
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/eac_v_egjac/"
outcome_path = "R_data/processed_records/outcome.rds"
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
# df = impute_srs(df, quantiles)

# =========================================================
# get outcomes
outcomes = master %>% select(id, cancertype)
outcomes %<>% mutate(
  ANY=ifelse(cancertype=="", 0, 1),
  EAC=ifelse(cancertype=="EAC", 1, 0),
  EGJAC=ifelse(cancertype=="EGJAC", 1, 0)
)
# merge into df
df %<>% left_join(outcomes, by="id")
outcome_names = c("ANY", "EAC", "EGJAC")


# =========================================================
# predicted risk per outcome
xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                 label=df$casecontrol)
proba = predict(xgb_fit, newdata=xgb_df)

risks = data.frame(cancertype=df$cancertype, risk=proba*100000)

filepath = paste0(dir_figures, "risk_distribution.pdf")
g = ggplot(data=risks, aes(x=risk, colour=cancertype, fill=cancertype)) + 
  geom_density(alpha=0.2) + xlab("Predicted risk (/100,000)") +
  ylab("Density") + xlim(0, 1000) + scale_x_continuous(trans="log10") + 
  ggtitle(paste0("Risk distribution"))
g
ggsave(filepath, g, width=6, height=3)

# =========================================================
# difference in features

cases = df %>% filter(casecontrol == 1)

cases = df

means = cases %>% group_by(cancertype) %>% select(xgb_fit$feature_names) %>%
  summarise_all(mean, na.rm=T)
means = data.frame(means)
rownames(means) = c("Control", "EAC", "EGJAC")
means$cancertype = NULL
means = data.frame(t(means))

mean_diff = means$EAC - means$EGJAC

sds = cases %>% group_by(cancertype) %>% select(xgb_fit$feature_names) %>%
  summarise_all(sd, na.rm=T)
sds = data.frame(sds)
rownames(sds) = c("Control", "EAC", "EGJAC")
sds$cancertype = NULL
sds = data.frame(t(sds))

ns = cases %>% group_by(cancertype) %>% select(xgb_fit$feature_names) %>%
  summarise(across(everything(), ~ sum(!is.na(.))))

test_df = data.frame(
  mean_control=means$Control,
  mean_EAC=means$EAC,
  mean_EGJAC=means$EGJAC,
  mean_diff=mean_diff,
  sd_EAC=sds$EAC,
  sd_EGJAC=sds$EGJAC,
  n_EAC=ns[2, 2:212] %>% unlist(., use.names=FALSE),
  n_EGJAC=ns[3, 2:212] %>% unlist(., use.names=FALSE)
)
test_df$prop_nonNA_EAC = test_df$n_EAC / max(test_df$n_EAC)
test_df$prop_nonNA_EGJAC = test_df$n_EGJAC / max(test_df$n_EGJAC)

test_df$sd_pooled = with(
  test_df, 
  sqrt(((n_EAC-1)*sd_EAC^2 +(n_EGJAC-1)*sd_EGJAC^2) / (n_EAC+n_EGJAC -2))
)
test_df$denum = with(
  test_df, 
  sd_pooled * sqrt(1/n_EAC + 1/n_EGJAC)
)
test_df$t_stat = with(
  test_df, 
  abs(mean_diff) / denum
)
test_df$pvalue = with(
  test_df, 
  pt(t_stat, n_EAC+n_EGJAC-2, lower.tail=F)*2
)
rownames(test_df) = xgb_fit$feature_names
test_df[, c("t_stat", "pvalue")]
test_df$pvalue.adj = p.adjust(test_df$pvalue, "fdr")

test_df %>% filter(pvalue.adj < 0.05) %>% 
  select(mean_control, mean_EAC, mean_EGJAC, pvalue.adj) %>%
  xtable::xtable(digits=3)

prop_df = test_df %>% 
  select(n_EAC, n_EGJAC) %>% rename(n_nonNA_EAC=n_EAC, n_nonNA_EGJAC=n_EGJAC) %>%
  mutate(n_EAC=max(n_nonNA_EAC), n_EGJAC=max(n_nonNA_EGJAC))

prop_df %<>% mutate(
  pvalue=prop.test(c(n_nonNA_EAC, n_nonNA_EGJAC), c(n_EAC, n_EGJAC))
)
prop_df$pvalue = sapply(seq(nrow(prop_df)), function(i) prop.test(
  c(prop_df$n_nonNA_EAC[i], prop_df$n_nonNA_EGJAC[i]),
  c(prop_df$n_EAC[i], prop_df$n_EGJAC[i])
  )$p.value)
prop_df$pvalue.adj = p.adjust(prop_df$pvalue, "fdr")

prop_df %<>% mutate(
  prop_nonNA_EAC=n_nonNA_EAC / n_EAC,
  prop_nonNA_EGJAC=n_nonNA_EGJAC / n_EGJAC
)

prop_df %>% filter(pvalue.adj < 0.05) %>% 
  select(prop_nonNA_EAC, prop_nonNA_EGJAC, pvalue.adj) %>%
  xtable::xtable(digits=3)

g = ggplot() + 
  geom_point(
    data=prop_df,
    mapping=aes(x=prop_nonNA_EAC, y=prop_nonNA_EGJAC, color=pvalue.adj)
  ) + xlim(0, 1) + ylim(0, 1) + 
  geom_abline(slope=1, intercept=0) + 
  scale_color_gradientn(name="adj. p value", colors=rainbow(3), 
                        values=c(0, 0.05, 1)) + 
  xlab("Prop. missing EAC")+ 
  ylab("Prop. missing EGJAC")

filepath = paste0(dir_figures, "vs_propNA.pdf")
ggsave(filepath, g, width=6, height=5)
# =========================================================
# SHAP values

dff = bind_rows(
  df %>% filter(casecontrol==1),
  df %>% filter(casecontrol==0) %>% sample_n(5000)
)
xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                     label=dff$casecontrol)
proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib = F)

ids = dff$EAC == 1
ids = dff$EGJAC == 1
ids = dff$casecontrol == 0

shapplot = xgboost:xgb.plot.shap(
  data=dff[ids, ] %>% select(xgb_fit$feature_names) %>% as.matrix(),
  shap_contrib=proba[ids, ], 
  top_n=9,
  # features=xgb_fit$feature_names[1:9], 
  n_col=3,
  model=xgb_fit,
  plot=T
)



# =========================================================
# age
filepath = paste0(dir_figures, "shap_age.pdf")

shapplot = xgboost::xgb.plot.shap(
  data=dff %>% select(xgb_fit$feature_names) %>% as.matrix(),
  shap_contrib=proba, 
  features="ageatindex", 
  model=xgb_fit,
  plot=F
)
shap = exp(shapplot$shap_contrib)
df_shap = data.frame(
  casecontrol=as.factor(dff$casecontrol),
  SHAP=shap,
  age=shapplot$data
)
colnames(df_shap) = c("casecontrol", "SHAP", "age")

g = ggplot(df_shap, aes(x=age, y=SHAP, group=casecontrol, color=casecontrol)) + 
  geom_point(alpha=0.1) + geom_smooth(method="gam") + ylab("exp(SHAP)")
filepath = paste0(dir_figures, "shap_age.pdf")
ggsave(filepath, g, width=8, height=5)  

tab = with(df, table(ageatindex, casecontrol))
tab = cbind(tab, rowSums(tab))
tab = cbind(tab, 100000*tab[, 2]/tab[, 3])
tab = data.frame(tab)
colnames(tab) = c("controls", "cases", "n", "case_prop")
tab$age = as.numeric(rownames(tab))
filepath = paste0(dir_figures, "cond_prop_age.pdf")
g = ggplot(tab, aes(x=age, y=case_prop)) + 
  geom_point() + ylab("Case proportion (/100,000)")
ggsave(filepath, g, width=6, height=5)  



proba = predict(xgb_fit, newdata=xgb_df)
df_risk = data.frame(
  casecontrol=as.factor(df$casecontrol),
  age=df$ageatindex,
  risk=proba*100000
)

df_risk_agg = df_risk %>% group_by(casecontrol, age) %>% summarise(
  med_risk=median(risk),
  q25_risk=quantile(risk, 0.25),
  q75_risk=quantile(risk, 0.75)
)

g = ggplot(df_risk_agg, aes(x=age, y=med_risk, 
                        group=casecontrol, color=casecontrol, fill=casecontrol,
                        ymin=q25_risk, ymax=q75_risk)) +
  geom_line() + ylab("Predicted risk") +
  geom_ribbon(alpha=0.2) + scale_y_continuous(trans="log10")

filepath = paste0(dir_figures, "risk_age.pdf")
ggsave(filepath, g, width=6, height=5)  



# ==============================================================================
# Comparing SHAP values

# read in data
file_path = paste0(dir_path, "5-1_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df
# subset to test set
df %<>% filter(id %in% test_ids)

# sample ids
ids = df %>% sample_n(50000) %>% pull(id)
dfn = df %>% filter(id %in% ids)
dfn %<>% left_join(outcomes, by="id")

# read in models
dir_model = "R_data/results/models/"
model_names = c(
  ANY="XGB_all_ANY.rds",
  EAC="XGB_all_EAC.rds",
  EGJAC="XGB_all_EGJAC.rds"
)
models = lapply(model_names, function(file){
  results = readRDS(paste0(dir_model, file))
  return(results)
})

# get SHAP values for all three models
shap_model = lapply(outcome_names, function(outcome){
  xgb_fit = models[[outcome]]$xgb_fit
  set.seed(0)
  dfni = impute_srs(dfn, models[[outcome]]$quantiles)
  dfni %<>% mutate(casecontrol = get(outcome))
  xgb_df = xgb.DMatrix(as.matrix(dfni %>% select(xgb_fit$feature_names)),
                        label=dfni$casecontrol)
  proba = predict(xgb_fit, newdata=xgb_df, predcontrib=TRUE, approxcontrib=F)
  return(proba)
})
names(shap_model) = outcome_names
# variable importance

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


shap_by_group = lapply(shap_model, function(proba){
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
  return(df_shap)
})

shap_by_group %<>% bind_rows(.id="model")

feature_order = shap_by_group %>% filter(model=="ANY") %>% arrange(SHAP) %>% pull(feature)

shap_by_group$feature %<>% factor(levels=feature_order, ordered=T)

# plot by group
g = ggplot(shap_by_group, aes(x=SHAP, y=feature, fill=model)) + 
  geom_bar(stat="identity", position="dodge") +
  xlab("mean|SHAP|") + ylab("") 
g
filepath = paste0(dir_figures, "shap_groups_cancertype.pdf")
ggsave(filepath, g, width=6, height=8)

# shap pdp
shap_model_ = lapply(shap_model, function(df) df %>% data.frame() %>% mutate(id=dfn$id))
shap_model_df = shap_model_ %>% bind_rows(.id="model")
for(var in xgb_fit$feature_names){
  deciles = quantile(dfn[[var]], (1:9)/10, na.rm=T)
  xlims = quantile(dfn[[var]], c(0.05, 0.95), na.rm=T)
  deciles = data.frame(x=deciles, xend=deciles, y=-0.1, yend=0)
  cat(paste0("Feature: ", var, "...\n"))
  plotdf = shap_model_df %>% 
    select(id, model, !!var) %>% 
    rename(shap=!!var) %>%
    left_join(dfn %>% select(id, casecontrol, !!var), by="id")
  plotdf$Case = plotdf$casecontrol %>% factor(levels=c(0,1), labels=c("Control", "Case"))
  plotdf %<>% mutate(shap=exp(shap))
  g = ggplot() + 
    geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
                 inherit.aes=F) + ylab("exp(SHAP)") + xlab(var) + 
    geom_smooth(
      data=plotdf, mapping=aes_string(x=var, y="shap", color="model"),
      method=ifelse(length(unique(dff[[var]]))>3, "gam", "lm")) +
    xlim(xlims)
  g
  filepath = paste0(dir_figures, "shap/", var, ".pdf")
  ggsave(filepath, g, width=6, height=4)
  cat("...done!\n")
}

# ==============================================================================
# comparison


# read in data
file_path = paste0(dir_path, "5-1_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df
# subset to test set
df %<>% filter(id %in% test_ids)

# sample ids
ids = df %>% sample_n(50000) %>% pull(id)
dfn = df %>% filter(id %in% ids)
dfn %<>% left_join(outcomes, by="id")

# read in models
dir_model = "R_data/results/models/"
model_names = c(
  ANY="XGB_all_ANY.rds",
  EAC="XGB_all_EAC.rds",
  EGJAC="XGB_all_EGJAC.rds"
)
models = lapply(model_names, function(file){
  results = readRDS(paste0(dir_model, file))
  return(results)
})

# get predictions for all three models
risk_model = lapply(outcome_names, function(outcome){
  xgb_fit = models[[outcome]]$xgb_fit
  set.seed(0)
  dfni = impute_srs(dfn, models[[outcome]]$quantiles)
  dfni %<>% mutate(casecontrol = get(outcome))
  xgb_df = xgb.DMatrix(as.matrix(dfni %>% select(xgb_fit$feature_names)),
                       label=dfni$casecontrol)
  proba = predict(xgb_fit, newdata=xgb_df)
  out = data.frame(
    id=dfni$id, 
    risk=log10(proba*100000)
  )
  return(out)
})
names(risk_model) = outcome_names

risk_model %<>% bind_rows(.id="model")

risk_model %<>% pivot_wider(id_cols=id, names_from=model, values_from=risk)

risk_model %<>% left_join(outcomes%>%select(id, cancertype), by="id")

g = GGally::ggpairs(
  data=risk_model %>% select(-id),
  mapping=aes(color=cancertype, alpha=0.2),
  columns=outcome_names
)
g
filepath = paste0(dir_figures, "risk_scatter.pdf")
ggsave(filepath, g, width=6, height=6)
