setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/classification_metrics.R')

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
# model_path = "R_data/results/models/XGB_nALL_typeANY.rds"
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
# function to get ROC from a df
get_roc = function(df){
  # ensure correct column ordering for xgb model
  df %<>% select(c(id, casecontrol, xgb_fit$feature_names))
  y = df$casecontrol
  # convert to xgb format
  df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                   label=df$casecontrol)
  # get predicted risk and ROC curve
  proba = predict(xgb_fit, newdata=df)
  fg = proba[y==1]; bg = proba[y==0]
  
  
  roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
  roc$curve = roc$curve[seq(1, nrow(roc$curve), by=ceiling(nrow(df)/1000)), ]
  roc$curve %<>% data.frame()
  colnames(roc$curve) = c("fpr", "recall", "tr")
  
  proc = pROC::roc(controls=bg, cases=fg)
  roc$ci = pROC::ci(proc, of="auc")
  roc$display.ci = paste0(
    round(roc$au, 3), " [",
    round(roc$ci[1], 3), ",",
    round(roc$ci[3], 3), "]"
  )
  roc$display = round(roc$au, 3)
  return(roc)
}

roc = get_roc(df)


# =========================================================
# ROC curve
roc$curve$tr = roc$curve$tr*100000
title = paste0("ROC curve ", "(AUC: ", roc$display.ci, ")")

g = ggplot(data=roc$curve, 
           aes(x=fpr, y=recall, color=tr)) + 
  geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  labs(color="Threshold\n(/100,000)") + 
  scale_color_gradientn(trans="log", colors=rainbow(6), breaks=c(1, 10, 100, 1000, 10000)) + 
  ggtitle(title)
filepath = paste0(dir_figures, "roc.pdf")
ggsave(filepath, g, width=5, height=4)

# =========================================================
# calibration plot

calib_50  = calibration_curve(xgb_fit, xgb_df, nbins=50)


calib_df = calib_50
calib_df = calib_df*100000
# calib_df$mid = calib_df$mid / 419

log = T
g = ggplot(data=calib_df, aes(x=mid, y=propcase)) + theme(aspect.ratio=1) + 
  geom_point()  +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ylab("Observed (/100,000)") + xlab("Predicted (/100,000)")+ 
  ggtitle(paste0("Calibration", ifelse(log, " (log-log)", "")))
if(log){
  g = g + scale_x_log10(limits=c(1, 3000)) + scale_y_log10(limits=c(1, 3000))
}else{
  g = g + xlim(0, 500) + ylim(0, 500)
}
g
filename = paste0(dir_figures, "calibration50_all",  ifelse(log, "_log", "_zoom"), ".pdf")
ggsave(filename, g, width=4, height=4)



log = F
g = ggplot(data=calib_df) + theme(aspect.ratio=1) + 
  geom_rect(aes(ymin=0, ymax=propcase, xmin=L, xmax=pmin(U, 100000)))  +
  geom_point(aes(x=mid, y=propcase, color="red")) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ylab("Observed (/100,000)") + xlab("Predicted (/100,000)") + 
  ggtitle(paste0("Calibration", ifelse(log, " (log-log)", ""))) +
  theme(legend.position="none")
if(log){
  g = g + scale_x_log10(limits=c(1, 100000)) + scale_y_log10(limits=c(1, 100000))
}else{
  g = g + xlim(0, 6000) + ylim(0, 6000)
}
g
filename = paste0(dir_figures, "calibration1000_bar",  ifelse(log, "_log", "_zoom"), ".pdf")
ggsave(filename, g, width=4, height=4)



# HL
calib_df = calibration_curve(xgb_fit, xgb_df, nbins=50)
nbins = nrow(calib_df) - 1
pred = calib_df$mid
obs = calib_df$propcase
n = calib_df$N
n1 = calib_df$Ncase

o1 = n1
o0 = n-n1
e1 = pred*n
e0 = n-e1

H51 = sum(((o1-e1)^2/e1 + (o0-e0)^2/e0))
df51 = nbins - 1
p51 = pchisq(H51, df51, lower.tail=F)
print(paste0("H=", H51, ", df=", df51, ", p=", p51))
which = 6:(nbins+1-5)
H50 = sum(((o1-e1)^2/e1 + (o0-e0)^2/e0)[which])
df50 = nbins - 10 - 1
p50 = pchisq(H50, df50, lower.tail=F)
print(paste0("H=", H50, ", df=", df50, ", p=", p50))


dff = df %>% select(c(id, casecontrol, xgb_fit$feature_names))
y = dff$casecontrol
# convert to xgb format
dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                 label=dff$casecontrol)
# get predicted risk and ROC curve
proba = predict(xgb_fit, newdata=dff)

pvals = sapply(c(10, 30, 50, 100, 200, 300, 500, 1000), 
               function(i) ResourceSelection::hoslem.test(y, proba, g=i)$p.value)




# =========================================================
# Threshold table
thresholds = sort(unique(c(
  seq(10, 50, 10), seq(50, 200, 5),
  seq(200, 300, 10), seq(300, 500, 25),
  seq(500, 2000, 100)
))) /100000
tr_df = classification_metrics(xgb_fit, xgb_df, thresholds)
tr_df = tr_df %>% select(one_of("tpr", "ppv", "detection_prevalance"))
# tr_df = tr_df[46:225, ]
rownames(tr_df) = format(round(as.numeric(rownames(tr_df)), 5)*100000, digits=5)
tr_df = tr_df*100
cat(print(xtable::xtable(tr_df)), 
    file=paste0(dir_figures, "all_calibration.tex"))
write.csv(tr_df, paste0(dir_figures, "all_calibration.csv"))

# =========================================================
# Quantile table
dff = df %>% select(c(id, casecontrol, xgb_fit$feature_names))
y = dff$casecontrol
# convert to xgb format
dff = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                  label=dff$casecontrol)
# get predicted risk and ROC curve
proba = predict(xgb_fit, newdata=dff)
qs = c(seq(0, 90, 1), seq(90.1, 99, 0.1), seq(99.01, 99.99, 0.01)) / 100
thresholds = quantile(proba, qs)

tr_df = classification_metrics(xgb_fit, xgb_df, thresholds)
tr_df %<>% select(one_of("tpr", "ppv", "detection_prevalance"))
tr_df = tr_df*100

out = tr_df %>% round(2)
rownames(out) = round(100*(1-qs), 2)
out$risk_threshold = round(thresholds * 100000)

write.csv(out, paste0(dir_figures, "calibration_quantiles.csv"))


# =========================================================
# Quantile table per age and sex

df0 = data.frame(
  casecontrol=df$casecontrol,
  age=df$age,
  sex=ifelse(df$gender==1, "M", "F"),
  risk=proba*100000
)

means = df0 %>% group_by(sex, age) %>% summarise(risk=mean(risk))
g = ggplot() + 
  geom_line(
    data=means, 
    mapping=aes(x=age, y=risk, color=sex)
  ) + 
  scale_y_log10(breaks=c(1, 2, 5, 10, 20, 50, 100, 200)) +
  xlab("Age") + ylab("Avg. pred. risk (/100,000)") + 
  ggtitle("Average predicted risk by sex and age")

filename = paste0(dir_figures, "avg_pred_risk_age_sex.pdf")
ggsave(filename, g, width=8, height=4)

df0 %>% group_by(sex) %>% summarise(risk=mean(risk))
