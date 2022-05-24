setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/birthyear/"
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
master$indexdate = master$end + 365
master %<>% left_join(df%>%select(id, age), by="id")
year_end = 18690 + 365
master %<>% mutate(indexyear = (indexdate - year_end)/365 + 2019)
hist(master$indexyear, breaks=2004:2020)
master %<>% mutate(birthyear = indexyear - age)
hist(master$birthyear)
master %<>% mutate(birth35 = birthyear < 1935)
master %<>% mutate(agebin = cut(age, c(70, 75, 80, 90)))

# prevalance
ctab = with(master, table(casecontrol, birth35, agebin))
prop_cases = 100000 * ctab[2,,] / (ctab[1,,] + ctab[2,,])
# proportion
ctab[1,,] + ctab[2,,]


dff = df %>% left_join(master%>%select(id, birthyear), by="id")
dff %<>% mutate(birth_lt_1935=birthyear<1935)
dff %<>% mutate(birth_cohort=cut(birthyear, c(1900, 1935, 1945, 1955, 1965, 2000)))
ctab2 = with(dff, table(age, birth_lt_1935, casecontrol))
ctab2 = with(dff, table(age, birth_cohort, casecontrol))

prop_cases = ctab2[,,2] *100000 / (ctab2[,,1]+ctab2[,,2])
prop_cases %<>% data.frame()
prop_cases$age = as.numeric(prop_cases$age)+17

g = ggplot(prop_cases, aes(x=age, y=Freq, color=birth_cohort)) +
  geom_line() + ylab("Prop. cases (/100,000)")
ggsave(paste0(dir_figures, "prop_cases_age_birth.pdf"), g, width=8, height=4)



g = ggplot(dff, aes(x=age, fill=birth_cohort)) +
  geom_histogram(binwidth=1, position="fill", alpha=0.8)
ggsave(paste0(dir_figures, "age_birth.pdf"), g, width=8, height=4)

# imputation
dff = df %>% filter(!is.na(age))
set.seed(0)
dff = impute_srs(dff, quantiles)
dff %<>% sample_n(50000)

xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                      label=dff$casecontrol)
proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib=F)
risk = predict(xgb_fit, newdata=xgb_dff)

plot_df = data.frame(
  id = dff$id,
  age = dff$age,
  SHAP_age = proba[, 1],
  risk = risk
)
plot_df %<>% left_join(master%>%select(id, birthyear), by="id")
plot_df %<>% mutate(birth_lt_1935=birthyear<1935)
plot_df %<>% mutate(birth_cohort=cut(birthyear, c(1900, 1935, 1945, 1955, 1965, 2000)))

g = ggplot(plot_df, aes(x=age, y=risk, color=birth_cohort)) +
  geom_point(alpha=0.05) + scale_y_continuous(trans="log10") +
  geom_smooth(method="gam")
ggsave(paste0(dir_figures, "risk_age_birth.pdf"), g, width=6, height=4)

g = ggplot(plot_df, aes(x=age, y=SHAP_age, color=birth_cohort)) +
  geom_point(alpha=0.05) +
  geom_smooth(method="gam") + ylim(0, 0.1)
ggsave(paste0(dir_figures, "shap_age_birth.pdf"), g, width=6, height=4)

# pdps

library(pdp)
dff = dff %>% left_join(master%>%select(id, birthyear), by="id")

dff_lt_1935 = dff %>% filter(birthyear<1935) %>% select(xgb_fit$feature_names)
dff_gt_1935 = dff %>% filter(birthyear>=1935) %>% select(xgb_fit$feature_names)

out_lt_1935 = pdp::partial(xgb_fit, pred.var="age", 
                           train=dff_lt_1935, 
                           type="classification", prob=T, which.class=1,
                           plot=F, progress="text")

out_gt_1935 = pdp::partial(xgb_fit, pred.var="age", 
                           train=dff_gt_1935, 
                           type="classification", prob=T, which.class=1,
                           plot=F, progress="text")

out_lt_1935$yhat = out_lt_1935$yhat * 100000
out_gt_1935$yhat = out_gt_1935$yhat * 100000
out_lt_1935$birth_lt_1935 = T
out_gt_1935$birth_lt_1935 = F
out = rbind(out_lt_1935, out_gt_1935)

g = ggplot(out, aes(x=age, y=yhat, color=birth_lt_1935)) + 
  geom_line() + ylim(0, max(out$yhat))+ 
  ylab("PDP (/100,000)")
filepath = paste0(dir_figures, "age_pdp.pdf")
ggsave(filepath, g, width=8, height=4)

# ==============================================================================
# = Fitted models ==============================================================
# ==============================================================================
MODELS_PATHS = list(
  NONE = "R_data/results/models/test/XGB_n1M_typeANY_yearNONE.rds",
  IND = "R_data/results/models/test/XGB_n1M_typeANY_yearIND.rds",
  NUM = "R_data/results/models/test/XGB_n1M_typeANY_yearNUM.rds",
  DROP = "R_data/results/models/test/XGB_n1M_typeANY_yearDROP.rds"
) 
models = lapply(names(MODELS_PATHS), function(name){
  file_path = MODELS_PATHS[[name]]
  return(readRDS(file_path))
})
names(models) = names(MODELS_PATHS)
lapply(models, function(model) model$xgb_fit)

# import data
complete_data = readRDS('R_data/processed_records/5-1_test_merged.rds')
master = complete_data$master
complete_data = complete_data$df

# subset to test ids
complete_data %<>% filter(id %in% models[["NONE"]]$test_ids)

# add in year
master$indexdate = master$end + 365
master %<>% left_join(complete_data %>% select(id, age), by="id")
year_end = 19055
master %<>% mutate(indexyear = (indexdate - year_end)/365 + 2019)
hist(master$indexyear)
master %<>% mutate(birthyear = indexyear - age)
hist(master$birthyear)
master %<>% mutate(birth35 = birthyear >= 1935)
complete_data %<>% left_join(master%>%select(id,birthyear,birth35), by="id")

# merge into single df
pred = HOSEA::predict.HOSEA(complete_data, 1, models)
df = complete_data %>% select(id, casecontrol, age, birthyear)
df %<>% left_join(pred, by="id")

df %<>% mutate(birth_cohort=cut(birthyear, c(1900, 1925, 1935, 1945, 1955, 1965, 1975, 2000)))

subsets = list(
  all=df%>%pull(id),
  cohort_1900_1925=df%>%filter(birth_cohort=="(1900,1925]")%>%pull(id),
  cohort_1925_1935=df%>%filter(birth_cohort=="(1925,1935]")%>%pull(id),
  cohort_1935_1945=df%>%filter(birth_cohort=="(1935,1945]")%>%pull(id),
  cohort_1945_1955=df%>%filter(birth_cohort=="(1945,1955]")%>%pull(id),
  cohort_1955_1965=df%>%filter(birth_cohort=="(1955,1965]")%>%pull(id),
  cohort_1965_1975=df%>%filter(birth_cohort=="(1965,1975]")%>%pull(id),
  cohort_1975_2000=df%>%filter(birth_cohort=="(1975,2000]")%>%pull(id)
)

auc_mat = sapply(names(subsets), function(sname){
  ids = subsets[[sname]]
  y = df%>%filter(id %in% ids) %>% pull(casecontrol)
  dff = complete_data %>% filter(id %in% ids)
  aucs = sapply(names(models), function(mname){
    xgb_fit = models[[mname]]$xgb_fit
    # convert to xgb format
    xgb_df = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                     label=dff$casecontrol)
    # get predicted risk and ROC curve
    proba = predict(xgb_fit, newdata=xgb_df)
    fg = proba[y==1]; bg = proba[y==0]
    roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
    return(roc$auc)
  })
}, simplify=T)

xtable::xtable(t(auc_mat), digits=3)

# shap and pdp
test_df = complete_data %>% sample_n(25000)
shap_dfs = list()
pdp_dfs = list()

for(mname in names(models)){
  xgb_fit = models[[mname]]$xgb_fit
  quantiles = models[[mname]]$quantiles
  set.seed(0)
  dff = impute_srs(test_df, quantiles)
  xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                        label=dff$casecontrol)
  proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib=F)
  shap_df = proba %>% data.frame() %>% mutate(id=dff$id) %>% select(id, age)
  pdp_df = pdp::partial(xgb_fit, pred.var="age", train=dff %>% select(xgb_fit$feature_names), 
               type="classification", prob=T, which.class=1,
               plot=F, progress="text")
  shap_dfs[[mname]] = shap_df
  pdp_dfs[[mname]] = data.frame(pdp_df)
}

shap_df = bind_rows(shap_dfs, .id="model")
pdp_df = bind_rows(pdp_dfs, .id="model")

# plot
shap_df %<>% left_join(test_df%>%select(id, age), by="id")
g = ggplot(shap_df, aes(x=age.y, y=age.x, color=model)) + 
  geom_point(alpha=0.05) +
  geom_smooth(method="gam") + xlab("age") + ylab("SHAP")
ggsave(paste0(dir_figures, "shap_models.pdf"), width=8, height=4)

pdp_df$yhat = 100000*pdp_df$yhat
g = ggplot(pdp_df, aes(x=age, y=yhat, color=model)) + 
  geom_line() + xlab("age") + ylab("PDP (/100,000)")
ggsave(paste0(dir_figures, "pdp_models.pdf"), width=8, height=4)

# pdp and shap for birth cohort
mname = "NUM"
xgb_fit = models[[mname]]$xgb_fit
quantiles = models[[mname]]$quantiles
set.seed(0)
dff = impute_srs(test_df, quantiles)
xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                      label=dff$casecontrol)
proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib=F)
shap_df_ = proba %>% data.frame() %>% mutate(id=dff$id) %>% select(id, age, birthyear)
pdp_df_ = pdp::partial(xgb_fit, pred.var="birthyear", train=dff %>% select(xgb_fit$feature_names), 
                      type="classification", prob=T, which.class=1,
                      plot=F, progress="text")

shap_df_ %<>% left_join(test_df%>%select(id, age, birthyear), by="id")

g=ggplot(shap_df_, aes(x=age.y, y=age.x)) + 
  geom_point(alpha=0.05) +
  geom_smooth(method="gam") + xlab("age") + ylab("SHAP")
ggsave(paste0(dir_figures, "shap_num_age.pdf"), width=8, height=4)
g=ggplot(shap_df_, aes(x=birthyear.y, y=birthyear.x)) + 
  geom_point(alpha=0.05) +
  geom_smooth(method="gam") + xlab("birth year") + ylab("SHAP")
ggsave(paste0(dir_figures, "shap_num_by.pdf"), width=8, height=4)
g=ggplot(shap_df_, aes(x=age.y, y=birthyear.x+age.x, color=birthyear.y)) + 
  geom_point(alpha=0.05) +
  geom_smooth(method="gam") + xlab("age") + ylab("SHAP(age+by)") +
  scale_color_gradientn(colors=rainbow(6))
ggsave(paste0(dir_figures, "shap_num_age_by.pdf"), width=8, height=4)


pdp_df_ = pdp::partial(xgb_fit, pred.var="birthyear", train=dff %>% select(xgb_fit$feature_names), 
                      type="classification", prob=T, which.class=1,
                      plot=F, progress="text")
pdp_df_$yhat = 100000*pdp_df_$yhat
g = ggplot(pdp_df_, aes(x=birthyear, y=yhat)) + 
  geom_line() + xlab("birth year") + ylab("PDP (/100,000)")
ggsave(paste0(dir_figures, "pdp_by.pdf"), width=8, height=4)

plot_df = with(complete_data, table(casecontrol, 
                                    cut(age, breaks=c(20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100)),
                                    cut(birthyear, seq(1915, 2000, 5))))
plot_df = plot_df[2, , ]/ (plot_df[1,,]+plot_df[2,,])
plot_df %<>% data.frame()
colnames(plot_df) = c("age", "birthyear", "propcases")
plot_df$birthyear = as.numeric(plot_df$birthyear)*5 + 1915+2.5

ggplot(plot_df, aes(x=birthyear, y=propcases, color=age)) + geom_line()


# ==============================================================================
# = Calibration
# ==============================================================================
dff = df %>% filter(!is.na(age))
set.seed(0)
dff = impute_srs(dff, quantiles)
dff %<>% left_join(master%>%select(id, birthyear), by="id")

dff %<>% mutate(
  birthcohort = cut(birthyear, seq(1910, 2000, 10),
                    labels=paste0(
                      "(",
                      seq(1910, 1990, 10),",",
                      seq(1920, 2000, 10), "]"
                    )),
  age_bin = cut(age, seq(20, 90, 10))
)

# birth cohort
calib_bc = lapply(dff%>%pull(birthcohort)%>%unique(), function(bc) {
  if(is.na(bc)) return(NULL)
  dfff = dff %>% filter(birthcohort==bc)
  xgb_df = xgb.DMatrix(
    as.matrix(dfff %>% select(xgb_fit$feature_names)), 
    label=dfff$casecontrol
  )
  calib  = calibration_curve(xgb_fit, xgb_df, nbins=50)
  calib$birthcohort = bc
  print(calib)
  return(calib)
})

calib_bc_merged = calib_bc %>% bind_rows()

log = F

g = ggplot(data=calib_bc_merged, aes(x=mid*100000, y=propcase*100000, color=birthcohort)) + 
  theme(aspect.ratio=1) + 
  geom_line(alpha=0.5)  +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ylab("Observed (/100,000)") + xlab("Predicted (/100,000)")+ 
  ggtitle(paste0("Calibration by Birth Cohort", ifelse(log, " (log-log)", "")))
if(log){
  g = g + scale_x_log10(limits=c(1, 15000)) + scale_y_log10(limits=c(1, 15000))
}else{
  g = g + xlim(50, 250) + ylim(50, 250)
}
g
filename = paste0(dir_figures, "calibration_birthcohort",  ifelse(log, "_log", "_zoom"), ".pdf")
ggsave(filename, g, width=4, height=4)


# age
calib_age = lapply(dff%>%pull(age_bin)%>%unique(), function(bc) {
  if(is.na(bc)) return(NULL)
  dfff = dff %>% filter(age_bin==bc)
  xgb_df = xgb.DMatrix(
    as.matrix(dfff %>% select(xgb_fit$feature_names)), 
    label=dfff$casecontrol
  )
  calib  = calibration_curve(xgb_fit, xgb_df, nbins=50)
  calib$age_bin = bc
  print(calib)
  return(calib)
})

calib_age_merged = calib_age %>% bind_rows()

log = T

g = ggplot(data=calib_age_merged, aes(x=mid*100000, y=propcase*100000, color=age_bin)) + 
  theme(aspect.ratio=1) + 
  geom_line(alpha=1.0)  +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ylab("Observed (/100,000)") + xlab("Predicted (/100,000)")+ 
  ggtitle(paste0("Calibration by Age group", ifelse(log, " (log-log)", "")))
if(log){
  g = g + scale_x_log10(limits=c(1, 1000)) + scale_y_log10(limits=c(1, 1000))
}else{
  g = g + xlim(0, 250) + ylim(0, 250)
}
g
filename = paste0(dir_figures, "calibration_age",  ifelse(log, "_log", "_zoom"), ".pdf")
ggsave(filename, g, width=4, height=4)