setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/evaluation_split.R')

# -----------------------------------------------
# original data (2--5)
# import data
complete_data = readRDS('R_data/complete_data_raw.rds')
complete_data$n_visits = NULL

# train-test split
set.seed(0)
out = train_test_split(df=complete_data, weights=c(3, 1))
train = out[[1]]
test = out[[2]]
rm(complete_data)

set.seed(0)
out = train_test_split(df=train, weights=c(2, 1))
rm(train);gc()
train_n = out[[1]]
# imputation
n_quantiles = 10000
quantiles = compute_quantiles(train_n, n_quantiles)
set.seed(0)
test_25 = impute_srs(test, quantiles)
rm(out, train_n, test);gc()
test_25 = xgb.DMatrix(as.matrix(test_25[-c(1,2)]),
                      label=test_25$CaseControl)
# -----------------------------------------------
# new data (4--5)
complete_data = readRDS('R_data/y45/complete_data_raw.rds')
complete_data$n_visits = NULL

# train-test split
set.seed(0)
out = train_test_split(df=complete_data, weights=c(3, 1))
test = out[[2]]
rm(complete_data)
# imputation
set.seed(0)
test_45 = impute_srs(test, quantiles)
rm(out, test);gc()
test_45 = xgb.DMatrix(as.matrix(test_45[-c(1,2)]),
                      label=test_45$CaseControl)

# -----------------------------------------------
# xgboost
xgb_fit = readRDS("R_data/results/models/resample_nall.rds")$xgb_fit


proba = predict(xgb_fit, newdata=test_25)
y = xgboost::getinfo(test_25, "label")
fg = proba[y==1]; bg = proba[y==0]
roc_25 = PRROC::roc.curve(fg, bg ,curve=TRUE)

proba = predict(xgb_fit, newdata=test_45)
y = xgboost::getinfo(test_45, "label")
fg = proba[y==1]; bg = proba[y==0]
roc_45 = PRROC::roc.curve(fg, bg ,curve=TRUE)

curves = data.frame(roc_45$curve[, 1:2], "4-5")
colnames(curves) = c("fpr", "recall", "years")
cdf = curves
curves = data.frame(roc_25$curve[, 1:2], "2-5")
colnames(curves) = c("fpr", "recall", "years")
cdf = rbind(cdf, curves)

library(ggplot2)
filepath = paste0("R_code/hosea-project/figures/window_roc.pdf")
g = ggplot(data=cdf, aes(x=fpr, y=recall, color=years)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted")
ggsave(filepath, g, width=8, height=7)
