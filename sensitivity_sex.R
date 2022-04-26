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
model_path = "R_data/results/models/XGB_nALL_typeANY.rds"
model_path = "R_data/results/models/XGB_n7M_typeANY.rds"
complete = T

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
master %<>% filter(ID %in% test_ids)

# =========================================================
# subset to EAC only (ie remove EGJAC)

egjac = master$CancerType == "EGJAC"
df %<>% filter(!egjac)

# =========================================================
# find complete cases for other methods

df$White = !(df$Asian | df$Black | df$HawaiianPacific | df$IndianAlaskan)
df$smoke_ever = df$smoke_current | df$smoke_former

complete_cases = !(
  is.na(df$Gender) |
    is.na(df$ageatindex) |
    is.na(df$smoke_former) |
    is.na(df$smoke_current) |
    is.na(df$bmi) |
    is.na(df$White) | # will be NA iff same for other columns
    is.na(df$GerdAtIndex) |
    is.na(df$h2r_max) |
    is.na(df$ppi_max)
)

if(complete){
  cdf = df %>% filter(complete_cases)
}else{ # need to impute
  set.seed(0)
  cdf = impute_srs(df, quantiles)
}

with(cdf, table(CaseControl, Gender))

# create the representative sample

ids_cases_female = cdf%>%filter(Gender==0,CaseControl==1)%>%pull(ID)
ids_cases_male = cdf%>%filter(Gender==1,CaseControl==1)%>%pull(ID)
ids_controls_female = cdf%>%filter(Gender==0,CaseControl==0)%>%pull(ID)
ids_controls_male = cdf%>%filter(Gender==1,CaseControl==0)%>%pull(ID)

n_cases_female = length(ids_cases_female)
n_controls_female = length(ids_controls_female)
n_cases_male = round(8.33*n_cases_female)
n_controls_male = n_controls_female + n_cases_female - n_cases_male

set.seed(0)
ids_cases_male = sample(ids_cases_male, n_cases_male, replace=F)
ids_controls_male = sample(ids_controls_male, n_controls_male, replace=F)

ids = c(
  ids_cases_female,
  ids_cases_male,
  ids_controls_female,
  ids_controls_male
)

cdf_sample = cdf %>% filter(ID %in% ids)
with(cdf_sample, table(CaseControl, Gender))

# test_complete is in the format for our xgboost model (w/ extra variables)
# the following keeps only the required info for kunzmann or hunt

# variables
cdf_sample2 = cdf_sample %>% select(c(
  "ID", "CaseControl", "smoke_current", "smoke_former", "Gender", 
  "GerdAtIndex", "White", "smoke_ever", "bmi", "ageatindex")
)

# check
cdf_sample2 %>% filter(CaseControl==1) %>% select(c(
  "ID", "CaseControl", "smoke_current", "smoke_former", "Gender", 
  "GerdAtIndex", "White", "smoke_ever", "bmi", "ageatindex"
)) %>% is.na() %>% colSums()

# compute various bins
cdf_sample2$k_age_bin = 
  relevel(cut(cdf_sample$ageatindex, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
cdf_sample2$k_bmi_bin = cut(cdf_sample$bmi, breaks=c(0, 25, 30, 35, 100))
cdf_sample2$h_age_bin = cut(cdf_sample$ageatindex, breaks=c(0, 50, 60, 70, 100))
cdf_sample2$h_bmi_bin = cut(cdf_sample$bmi, breaks=c(0, 30, 100))
cdf_sample2$h_smoke_any = pmax(
  cdf_sample$smoke_former,
  cdf_sample$smoke_current
)
cdf_sample2$k_ec = pmax(
  cdf_sample$h2r_max>0, 
  cdf_sample$ppi_max>0, 
  cdf_sample$GerdAtIndex
) > 0

# drop extra columns just in case
cdf_sample$White = NULL
cdf_sample$smoke_ever = NULL
# Kunzmann score
kunzmann_score = function(df){
  score = 0
  score = score + (df$k_age_bin == "(55,60]") * 1.5
  score = score + (df$k_age_bin == "(60,65]") * 2.5
  score = score + (df$k_age_bin == "(65,100]") * 3.5
  score = score + df$Gender * 4
  score = score + (df$k_bmi_bin == "(25,30]") * 1
  score = score + (df$k_bmi_bin == "(30,35]") * 1.5
  score = score + (df$k_bmi_bin == "(35,100]") * 2.5
  score = score + df$smoke_former * 2.5
  score = score + df$smoke_current * 3.5
  score = score + df$k_ec * 1.5
  return(score)
}
kscores = kunzmann_score(cdf_sample2)

hunt_score = function(df){
  score = 3.6
  score = score * ifelse(df$Gender, 1.9, 1.)
  score = score * ifelse(df$h_age_bin == "(50,60]", 2.1, 1.)
  score = score * ifelse(df$h_age_bin == "(60,70]", 3.2, 1.)
  score = score * ifelse(df$h_age_bin == "(70,100]", 3.1, 1.)
  score = score * ifelse(df$h_bmi_bin == "(30,100]", 1.8, 1.)
  score = score * ifelse(df$GerdAtIndex, 3.7, 1.)
  score = score * ifelse(df$h_smoke_any, 2.1, 1.)
  return(score/100000)
}
hscores = hunt_score(cdf_sample2)


# screening guidelines
scores = cdf_sample2 %>% mutate( 
  ACG2016     = Gender & GerdAtIndex & ( ((ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) ) > 1),
  ACG2022     = GerdAtIndex & ( (Gender + (ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) ) > 2),
  ACP2012     = Gender & GerdAtIndex & (ageatindex>=50) & ( ((bmi>30) + (smoke_ever) ) > 0),
  AGA2011     = ( GerdAtIndex + Gender + (ageatindex>=50) + (White) + (bmi>30) ) > 1,
  AGA_CPU2022 = ( GerdAtIndex + Gender + (ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) ) > 2,
  ASGE2019    = GerdAtIndex & ( (Gender + (ageatindex>=50) + (bmi>30) + (smoke_ever) ) > 0),
  BSG2013     = GerdAtIndex & ( (Gender + (ageatindex>=50) + (White) + (bmi>30) ) > 2),
  ESGE2017    = GerdAtIndex & ( (Gender + (ageatindex>=50) + (White) + (bmi>30) ) > 1)
)
guidelines = tail(colnames(scores), 8)

# AUCs
y = cdf_sample$CaseControl

proba = kscores
fg = proba[y==1]; bg = proba[y==0]
k_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)

proba = hscores
fg = proba[y==1]; bg = proba[y==0]
h_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)

xgb_df = xgb.DMatrix(as.matrix(cdf_sample %>% select(xgb_fit$feature_names)),
                     label=cdf_sample$CaseControl)
proba = predict(xgb_fit, newdata=xgb_df)
fg = proba[y==1]; bg = proba[y==0]
x_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)


curves = data.frame(x_roc$curve[, 1:2], 
                    paste0("xgboost (AUC=", round(x_roc$auc, 3),")"))
colnames(curves) = c("fpr", "recall", "method")
cdf = curves
curves = data.frame(k_roc$curve[, 1:2], 
                    paste0("Kunzmann (AUC=", round(k_roc$auc, 3),")"))
colnames(curves) = c("fpr", "recall", "method")
cdf = rbind(cdf, curves)
curves = data.frame(h_roc$curve[, 1:2], 
                    paste0("HUNT (AUC=", round(h_roc$auc, 3),")"))
colnames(curves) = c("fpr", "recall", "method")
cdf = rbind(cdf, curves)
cdf$label = ""

# guideline points
guide_roc = data.frame(t(sapply(guidelines,
                                function(g) PRROC::roc.curve(
                                  scores %>% select(g) %>% subset(y==1) %>% pull(g), 
                                  scores %>% select(g) %>% subset(y==0) %>% pull(g),
                                  curve=TRUE)$curve[2, 1:2]
)))
colnames(guide_roc) = c("fpr", "recall")
guide_roc$method = "guideline"
guide_roc$label = rownames(guide_roc)

guide_roc %<>% arrange(fpr)

guide_roc$xlab = c(.2, .27, .34, .41, .48, .55, .8, .9)
guide_roc$ylab = c(.15, .2, .25, .3, .35, .4, .85, .9)

filepath = paste0("R_code/hosea-project/figures/comparison_sex_", 
                  ifelse(complete, "complete", "imputed"), ".pdf")
g = ggplot(data=cdf, aes(x=fpr, y=recall, color=method)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  geom_point(data=guide_roc) + 
  ggtitle(paste0("Representative sample (sex): ", ifelse(complete, "complete", "imputed")))
# g = g +
#   geom_segment(data=guide_roc, aes(x=fpr, xend=xlab, y=recall, yend=ylab)) + 
#   geom_label(data=guide_roc, aes(label=label, x=xlab, y=ylab), size=3)
ggsave(filepath, g, width=8, height=6)
