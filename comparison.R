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

# =========================================================
# subset to outcomes

outcomes = master %>% select(id, casecontrol, cancertype)
outcomes %<>% mutate(
  ANY=casecontrol,
  EAC=as.integer(cancertype=="EAC"),
  EGJAC=as.integer(cancertype=="EGJAC")
)

outcome = "ANY"

outcomes_ = outcomes%>%select(id, !!outcome)
df = df %>% 
  left_join(outcomes_, by="id") %>%
  select(-casecontrol) %>% 
  rename(casecontrol=!!outcome)

# =========================================================
# find complete cases for other methods

df$white = !(df$asian | df$black | df$hawaiianpacific | df$indianalaskan)
df$smoke_ever = df$smoke_current | df$smoke_former

complete_cases = !(
  is.na(df$gender) |
    is.na(df$age) |
    is.na(df$smoke_former) |
    is.na(df$smoke_current) |
    is.na(df$bmi) |
    is.na(df$white) | # will be NA iff same for other columns
    is.na(df$gerd) |
    (is.na(df$h2r_max) &
    is.na(df$ppi_max))
)

test_complete = df %>% subset(complete_cases)
test_complete$casecontrol %>% table()

# test_complete is in the format for our xgboost model (w/ extra variables)
# the following keeps only the required info for kunzmann or hunt

# variables
df = test_complete %>% select(c(
  "id", "casecontrol", "smoke_current", "smoke_former", "gender", 
  "gerd", "white", "smoke_ever", "bmi", "age")
)

# check
df %>% filter(casecontrol==1) %>% select(c(
  "id", "casecontrol", "smoke_current", "smoke_former", "gender", 
  "gerd", "white", "smoke_ever", "bmi", "age"
  )) %>% is.na() %>% colSums()

# compute various bins
df$k_age_bin = 
  relevel(cut(test_complete$age, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
df$k_bmi_bin = cut(test_complete$bmi, breaks=c(0, 25, 30, 35, 100))
df$h_age_bin = cut(test_complete$age, breaks=c(0, 50, 60, 70, 100))
df$h_bmi_bin = cut(test_complete$bmi, breaks=c(0, 30, 100))
df$h_smoke_any = pmax(
  test_complete$smoke_former,
  test_complete$smoke_current
)
df$k_ec = pmax(
  test_complete$h2r_max>0, 
  test_complete$ppi_max>0, 
  test_complete$gerd, na.rm=T
) > 0

# drop extra columns just in case
test_complete$white = NULL
test_complete$smoke_ever = NULL
# Kunzmann score
kunzmann_score = function(df){
  score = 0
  score = score + (df$k_age_bin == "(55,60]") * 1.5
  score = score + (df$k_age_bin == "(60,65]") * 2.5
  score = score + (df$k_age_bin == "(65,100]") * 3.5
  score = score + df$gender * 4
  score = score + (df$k_bmi_bin == "(25,30]") * 1
  score = score + (df$k_bmi_bin == "(30,35]") * 1.5
  score = score + (df$k_bmi_bin == "(35,100]") * 2.5
  score = score + df$smoke_former * 2.5
  score = score + df$smoke_current * 3.5
  score = score + df$k_ec * 1.5
  return(score)
}
kscores = kunzmann_score(df)

hunt_score = function(df){
  score = 3.6
  score = score * ifelse(df$gender, 1.9, 1.)
  score = score * ifelse(df$h_age_bin == "(50,60]", 2.1, 1.)
  score = score * ifelse(df$h_age_bin == "(60,70]", 3.2, 1.)
  score = score * ifelse(df$h_age_bin == "(70,100]", 3.1, 1.)
  score = score * ifelse(df$h_bmi_bin == "(30,100]", 1.8, 1.)
  score = score * ifelse(df$gerd, 3.7, 1.)
  score = score * ifelse(df$h_smoke_any, 2.1, 1.)
  return(score/100000)
}
hscores = hunt_score(df)


# screening guidelines
scores = df %>% mutate( 
   ACG2016     = gender & gerd & ( ((age>=50) + (white) + (bmi>30) + (smoke_ever) ) > 1),
   ACG2022     = gerd & ( (gender + (age>=50) + (white) + (bmi>30) + (smoke_ever) ) > 2),
   ACP2012     = gender & gerd & (age>=50) & ( ((bmi>30) + (smoke_ever) ) > 0),
   AGA2011     = ( gerd + gender + (age>=50) + (white) + (bmi>30) ) > 1,
   AGA_CPU2022 = ( gerd + gender + (age>=50) + (white) + (bmi>30) + (smoke_ever) ) > 2,
   ASGE2019    = gerd & ( (gender + (age>=50) + (bmi>30) + (smoke_ever) ) > 0),
   BSG2013     = gerd & ( (gender + (age>=50) + (white) + (bmi>30) ) > 2),
   ESGE2017    = gerd & ( (gender + (age>=50) + (white) + (bmi>30) ) > 1)
)
guidelines = tail(colnames(scores), 8)
# scores %>% select(casecontrol, guidelines) %>% group_by(casecontrol) %>% summarise_all(sum)
# 
# Hmisc::describe(df %>% select(-c(id, bmi, age)))
# 
# n_risk_factors = df %>% 
#   mutate(n_risk_factors=gender+gerd+(age>=50) + (white) + (bmi>30) + (smoke_ever)) %>%
#   group_by(casecontrol, n_risk_factors) %>%
#   summarize(n=n())
# 
# scores %>% select(guidelines) %>% colSums()
# 
# rbind(n_risk_factors, p_risk_factors)

# AUCs
y = df$casecontrol

proba = kscores
fg = proba[y==1]; bg = proba[y==0]
k_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)

proba = hscores
fg = proba[y==1]; bg = proba[y==0]
h_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)

xgb_df = xgb.DMatrix(as.matrix(test_complete %>% select(xgb_fit$feature_names)),
                     label=test_complete$casecontrol)
proba = predict(xgb_fit, newdata=xgb_df)
fg = proba[y==1]; bg = proba[y==0]
x_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)


curves = data.frame(x_roc$curve[, 1:2], "xgboost")
colnames(curves) = c("fpr", "recall", "method")
cdf = curves
curves = data.frame(k_roc$curve[, 1:2], "Kunzmann")
colnames(curves) = c("fpr", "recall", "method")
cdf = rbind(cdf, curves)
curves = data.frame(h_roc$curve[, 1:2], "HUNT")
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

guide_roc$xlab = c(.3, .37, .44, .51, .58, .65, .8, .9)
guide_roc$ylab = c(.15, .2, .25, .3, .35, .4, .85, .9)

filepath = paste0("R_code/hosea-project/figures/comparison_", outcome, ".pdf")
g = ggplot(data=cdf, aes(x=fpr, y=recall, color=method)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  geom_point(data=guide_roc) + 
  ggtitle(paste0("Cancer type: ", outcome))
g = g +
  geom_segment(data=guide_roc, aes(x=fpr, xend=xlab, y=recall, yend=ylab)) + 
  geom_label(data=guide_roc, aes(label=label, x=xlab, y=ylab), size=3)
ggsave(filepath, g, width=6, height=5)
