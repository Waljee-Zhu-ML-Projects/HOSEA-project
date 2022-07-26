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
dir_model = "R_data/results/models/"
dir_figures = "R_code/hosea-project/figures/"
dir_results = "R_data/results/analyses/"
outcome_path = "R_data/processed_records/outcome.rds"

# =========================================================
# read in models
model_names = c(
  ANY="XGB_all_ANY.rds",
  EAC="XGB_all_EAC.rds",
  EGJAC="XGB_all_EGJAC.rds"
)
models = lapply(model_names, function(file){
  results = readRDS(paste0(dir_model, file))
  return(results)
})

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df

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

outcome = "EGJAC"

# =========================================================
# subset to test set
df %<>% filter(id %in% models[[outcome]]$test_ids)
# imputation
set.seed(0)
df = impute_srs(df, models[[outcome]]$quantiles)

df %<>% mutate(casecontrol = get(outcome))

xgb_df = xgb.DMatrix(as.matrix(df %>% select(xgb_fit$feature_names)),
                     label=df$casecontrol)




# =========================================================
# calibration plot

calib_50  = calibration_curve(xgb_fit, xgb_df, nbins=50)
calib_df = calib_50
calib_df = calib_df*100000

# calib_df$propcase = calib_df$propcase*9 #EGJAC
# calib_df$propcase = calib_df$propcase*2.25 #EAC

log = T
g = ggplot(data=calib_df) + theme(aspect.ratio=1) + 
  geom_rect(aes(ymin=0, ymax=propcase, xmin=L, xmax=pmin(U, 100000)))  +
  geom_point(aes(x=mid, y=propcase, color="red")) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ylab("Observed (/100,000)") + xlab("Predicted (/100,000)") + 
  ggtitle(paste0("Calibration ", outcome, ifelse(log, " (log-log)", ""))) +
  theme(legend.position="none")
if(log){
  g = g + scale_x_log10(limits=c(1, 100000)) + scale_y_log10(limits=c(1, 100000))
}else{
  g = g + xlim(0, 500) + ylim(0, 500)
}
g
filename = paste0(dir_figures, outcome, "_calibration50_bar",  ifelse(log, "_log", "_zoom"), ".pdf")
ggsave(filename, g, width=4, height=4)



# HL
calib_df = calib_50
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
nskip = 20
H50 = sum(((o1-e1)^2/e1 + (o0-e0)^2/e0)[1:(nbins+1-nskip)])
df50 = nbins - nskip - 1
p50 = pchisq(H50, df50, lower.tail=F)
print(paste0("H=", H50, ", df=", df50, ", p=", p50))


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

