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
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
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
master$indexdate = master$end + 365
master %<>% left_join(df%>%select(ID, ageatindex), by="ID")
year_end = 18690
master %<>% mutate(indexyear = (indexdate - year_end)/365 + 2019)
hist(master$indexyear)
master %<>% mutate(birthyear = indexyear - ageatindex)
hist(master$birthyear)
master %<>% mutate(birth35 = birthyear < 1935)
master %<>% mutate(agebin = cut(ageatindex, c(70, 75, 80, 90)))

# prevalance
ctab = with(master, table(CaseControl, birth35, agebin))
prop_cases = 100000 * ctab[2,,] / (ctab[1,,] + ctab[2,,])
# proportion
ctab[1,,] + ctab[2,,]


dff = df %>% left_join(master%>%select(ID, birthyear), by="ID")
dff %<>% mutate(birth_lt_1935=birthyear<1935)
ctab2 = with(dff, table(ageatindex, birth_lt_1935, CaseControl))

prop_cases = ctab2[,,2] *100000 / (ctab2[,,1]+ctab2[,,2])
prop_cases %<>% data.frame()
prop_cases$ageatindex = as.numeric(prop_cases$ageatindex)+17

g = ggplot(prop_cases, aes(x=ageatindex, y=Freq, color=birth_lt_1935)) +
  geom_line() + ylab("Prop. cases (/100,000)")
ggsave(paste0(dir_figures, "prop_cases_age_birth.pdf"), g, width=8, height=4)



g = ggplot(dff, aes(x=ageatindex, fill=birth_lt_1935)) +
  geom_histogram(binwidth=1, position="fill", alpha=0.8)
ggsave(paste0(dir_figures, "age_birth.pdf"), g, width=8, height=4)

# imputation
dff = df %>% filter(ageatindex >= 60)
set.seed(0)
dff = impute_srs(dff, quantiles)
dff %<>% sample_n(20000)

xgb_dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                      label=dff$CaseControl)
proba = predict(xgb_fit, newdata=xgb_dff, predcontrib=TRUE, approxcontrib=F)
risk = predict(xgb_fit, newdata=xgb_dff)

plot_df = data.frame(
  ID = dff$ID,
  ageatindex = dff$ageatindex,
  SHAP_ageatindex = proba[, 1],
  risk = risk
)
plot_df %<>% left_join(master%>%select(ID, birthyear), by="ID")
plot_df %<>% mutate(birth_lt_1935=birthyear<1935)

g = ggplot(plot_df, aes(x=ageatindex, y=risk, color=birth_lt_1935)) +
  geom_point(alpha=0.05) + scale_y_continuous(trans="log10") +
  geom_smooth()
ggsave(paste0(dir_figures, "risk_age_birth.pdf"), g, width=6, height=4)

g = ggplot(plot_df, aes(x=ageatindex, y=SHAP_ageatindex, color=birth_lt_1935)) +
  geom_point(alpha=0.05) +
  geom_smooth()
ggsave(paste0(dir_figures, "shape_age_birth.pdf"), g, width=6, height=4)

# pdps

library(pdp)
dff = dff %>% left_join(master%>%select(ID, birthyear), by="ID")

dff_lt_1935 = dff %>% filter(birthyear<1935) %>% select(xgb_fit$feature_names)
dff_gt_1935 = dff %>% filter(birthyear>=1935) %>% select(xgb_fit$feature_names)

out_lt_1935 = pdp::partial(xgb_fit, pred.var="ageatindex", 
                           train=dff_lt_1935, 
                           type="classification", prob=T, which.class=1,
                           plot=F, progress="text")

out_gt_1935 = pdp::partial(xgb_fit, pred.var="ageatindex", 
                           train=dff_gt_1935, 
                           type="classification", prob=T, which.class=1,
                           plot=F, progress="text")

out_lt_1935$yhat = out_lt_1935$yhat * 100000
out_gt_1935$yhat = out_gt_1935$yhat * 100000
out_lt_1935$birth_lt_1935 = T
out_gt_1935$birth_lt_1935 = F
out = rbind(out_lt_1935, out_gt_1935)

g = ggplot(out, aes(x=ageatindex, y=yhat, color=birth_lt_1935)) + 
  geom_line() + ylim(0, max(out$yhat))+ 
  ylab("PDP (/100,000)")
if(length(unique(df[[var]]))>100) g = g + xlim(quantile(df[[var]], c(0.05, 0.95)))
g
filepath = paste0(dir_figures, "pdp/", var, ".pdf")
ggsave(filepath, g, width=3, height=4)

