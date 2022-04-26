setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
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
model_path = "R_data/results/models/XGB_n7M_typeANY.rds"

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
master = df$master
df = df$df
# subset to test set
df %<>% filter(ID %in% test_ids)
# imputation
set.seed(0)
df = impute_srs(df, quantiles)

# =========================================================
# get outcomes
outcomes = master %>% select(ID, CancerType)
outcomes %<>% mutate(
  ANY=ifelse(CancerType=="", 0, 1),
  EAC=ifelse(CancerType=="EAC", 1, 0),
  EGJAC=ifelse(CancerType=="EGJAC", 1, 0)
)
# merge into df
df %<>% left_join(outcomes, by="ID")
outcome_names = c("ANY", "EAC", "EGJAC")


# =========================================================
# predicted risk per outcome
xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                 label=df$CaseControl)
proba = predict(xgb_fit, newdata=xgb_df)

risks = data.frame(CancerType=df$CancerType, risk=proba*100000)

filepath = paste0(dir_figures, "risk_distribution.pdf")
g = ggplot(data=risks, aes(x=risk, colour=CancerType, fill=CancerType)) + 
  geom_density(alpha=0.2) + xlab("Predicted risk (/100,000)") +
  ylab("Density") + xlim(0, 1000) +# scale_x_continuous(trans="log10") + 
  ggtitle(paste0("Risk distribution"))
g
ggsave(filepath, g, width=6, height=3)

# =========================================================
# difference in features

cases = df %>% filter(CaseControl == 1)

cases = df

means = cases %>% group_by(CancerType) %>% select(xgb_fit$feature_names) %>%
  summarise_all(mean)
means = data.frame(means)
rownames(means) = c("Control", "EAC", "EGJAC")
means$CancerType = NULL
means = data.frame(t(means))

mean_diff = means$EAC - means$EGJAC

sds = cases %>% group_by(CancerType) %>% select(xgb_fit$feature_names) %>%
  summarise_all(sd)
sds = data.frame(sds)
rownames(sds) = c("Control", "EAC", "EGJAC")
sds$CancerType = NULL
sds = data.frame(t(sds))

ns = (cases %>% group_by(CancerType) %>% summarize(n=n()))$n

test_df = data.frame(
  mean_control=means$Control,
  mean_EAC=means$EAC,
  mean_EGJAC=means$EGJAC,
  mean_diff=mean_diff,
  sd_EAC=sds$EAC,
  sd_EGJAC=sds$EGJAC,
  n_EAC=rep(ns[1], 2088),
  n_EGJAC=rep(ns[2], 761)
)

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

# black and female!?
tab = test_df[1:40, 1:4]
xtable::xtable(tab, digits=3)

xtable::xtable(test_df[test_df$pvalue < 0.01 / nrow(test_df), c(1:4)], digits=3)

t.test(
  df %>% filter(EAC==1) %>% pull(baso_max),
  df %>% filter(EGJAC==1) %>% pull(baso_max)
)

# =========================================================
# SHAP values

dff = bind_rows(
  df %>% filter(CaseControl==1),
  df %>% filter(CaseControl==0) %>% sample_n(5000)
)
xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                     label=dff$CaseControl)
proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib = F)

ids = dff$EAC == 1
ids = dff$EGJAC == 1
ids = dff$CaseControl == 0

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
  CaseControl=as.factor(dff$CaseControl),
  SHAP=shap,
  age=shapplot$data
)
colnames(df_shap) = c("CaseControl", "SHAP", "age")

g = ggplot(df_shap, aes(x=age, y=SHAP, group=CaseControl, color=CaseControl)) + 
  geom_point(alpha=0.1) + geom_smooth(method="gam") + ylab("exp(SHAP)")
filepath = paste0(dir_figures, "shap_age.pdf")
ggsave(filepath, g, width=8, height=5)  

tab = with(df, table(ageatindex, CaseControl))
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
  CaseControl=as.factor(df$CaseControl),
  age=df$ageatindex,
  risk=proba*100000
)

df_risk_agg = df_risk %>% group_by(CaseControl, age) %>% summarise(
  med_risk=median(risk),
  q25_risk=quantile(risk, 0.25),
  q75_risk=quantile(risk, 0.75)
)

g = ggplot(df_risk_agg, aes(x=age, y=med_risk, 
                        group=CaseControl, color=CaseControl, fill=CaseControl,
                        ymin=q25_risk, ymax=q75_risk)) +
  geom_line() + ylab("Predicted risk") +
  geom_ribbon(alpha=0.2) + scale_y_continuous(trans="log10")

filepath = paste0(dir_figures, "risk_age.pdf")
ggsave(filepath, g, width=6, height=5)  
