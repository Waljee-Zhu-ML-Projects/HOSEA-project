setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(pdp)
library(cowplot)
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
model_path = "R_data/results/models/XGB_nALL_typeANY.rds"

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

xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                     label=df$CaseControl)


# =========================================================
# 2d pdp


train = df %>% select(xgb_fit$feature_names)
set.seed(0) # sample some rows, otherwise waaaay too long
train = train[sample.int(nrow(train), 10000), ]


# univariate pdps

vars = c("hgb_max", "hct_max")

var = vars[1]
deciles = quantile(df[[var]], (1:9)/10)
deciles = data.frame(x=deciles, xend=deciles, y=0, yend=10)
out = pdp::partial(xgb_fit, pred.var=var, train=train, 
                   type="classification", prob=T, which.class=1,
                   plot=F, progress="text")
out$yhat = out$yhat * 100000
g1 = ggplot(out, aes(x=get(var), y=yhat)) + 
  geom_line() + ylim(0, 300)+ 
  ylab("Predicted risk (/100,000)") + xlab(var) +
  ggtitle(paste0("PDP: ", var)) +
  geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
               inherit.aes=F)
ggsave(paste0(dir_figures, "pdp_hgb.pdf"), width=4, height=4)

var = vars[2]
deciles = quantile(df[[var]], (1:9)/10)
deciles = data.frame(x=deciles, xend=deciles, y=0, yend=10)
out = pdp::partial(xgb_fit, pred.var=var, train=train, 
                   type="classification", prob=T, which.class=1,
                   plot=F, progress="text")
out$yhat = out$yhat * 100000
g2 = ggplot(out, aes(x=get(var), y=yhat)) + 
  geom_line() + ylim(0, 300)+ 
  ylab("Predicted risk (/100,000)") + xlab(var) +
  ggtitle(paste0("PDP: ", var)) +
  geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
               inherit.aes=F)
ggsave(paste0(dir_figures, "pdp_hct.pdf"), width=4, height=4)

# bivariate pdp
out2 = pdp::partial(xgb_fit, pred.var=vars, train=train, 
                    type="classification", prob=T, which.class=1,
                    plot=F, progress="text")

out2$yhat = out2$yhat * 100000
g3 = ggplot(out2, aes(x=hgb_max, y=hct_max, fill=yhat)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") + 
  xlim(11,16) + ylim(33, 46)
ggsave(paste0(dir_figures, "pdp_hct_hgb.pdf"), width=4, height=4)







vars = c("A1c_mean", "gluc_mean")

var = vars[1]
deciles = quantile(df[[var]], (1:9)/10)
deciles = data.frame(x=deciles, xend=deciles, y=0, yend=10)
out = pdp::partial(xgb_fit, pred.var=var, train=train, 
                   type="classification", prob=T, which.class=1,
                   plot=F, progress="text")
out$yhat = out$yhat * 100000
g1 = ggplot(out, aes(x=get(var), y=yhat)) + 
  geom_line() + ylim(0, 300)+ 
  ylab("Predicted risk (/100,000)") + xlab(var) +
  ggtitle(paste0("PDP: ", var)) +
  geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
               inherit.aes=F)
ggsave(paste0(dir_figures, "pdp_A1c_mean.pdf"), width=4, height=4)

var = vars[2]
deciles = quantile(df[[var]], (1:9)/10)
deciles = data.frame(x=deciles, xend=deciles, y=0, yend=10)
out = pdp::partial(xgb_fit, pred.var=var, train=train, 
                   type="classification", prob=T, which.class=1,
                   plot=F, progress="text")
out$yhat = out$yhat * 100000
g2 = ggplot(out, aes(x=get(var), y=yhat)) + 
  geom_line() + ylim(0, 300)+ 
  ylab("Predicted risk (/100,000)") + xlab(var) +
  ggtitle(paste0("PDP: ", var)) +
  geom_segment(data=deciles, aes(x=x, y=y, xend=xend, yend=yend), 
               inherit.aes=F)
ggsave(paste0(dir_figures, "pdp_gluc_mean.pdf"), width=4, height=4)

# bivariate pdp
out2 = pdp::partial(xgb_fit, pred.var=vars, train=train, 
                    type="classification", prob=T, which.class=1,
                    plot=F, progress="text")

out2$yhat = out2$yhat * 100000
g3 = ggplot(out2, aes(x=A1c_mean, y=gluc_mean, fill=yhat)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") + 
  xlim(5,8) + ylim(80, 150)
ggsave(paste0(dir_figures, "pdp_A1c_gluc.pdf"), width=4, height=4)


