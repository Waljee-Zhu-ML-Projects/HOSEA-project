setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/evaluation_split.R')
source('R_code/hosea-project/screening_guidelines.R')

# import data
complete_data = readRDS('R_data/complete_data_raw.rds')
complete_data$n_visits = NULL

log = function(df) cat(paste("Full data set: ", nrow(df), "observations,", 
                             df$CaseControl%>%sum, "cases", 
                             (df$CaseControl==0)%>%sum, "controls"), fill=T)

# reconstruct test set
set.seed(0)
complete_data = train_test_split(df=complete_data, weights=c(3, 1))[[2]]

# subset to complete cases
complete_data$White = !(complete_data$Asian | complete_data$Black | 
                          complete_data$HawaiianPacific | complete_data$IndianAlaskan)
complete_data$smoke_ever = complete_data$smoke_current | complete_data$smoke_former

complete_cases = !(
  is.na(complete_data$Gender) |
    is.na(complete_data$ageatindex) |
    is.na(complete_data$smoke_former) |
    is.na(complete_data$smoke_current) |
    is.na(complete_data$bmi) |
    is.na(complete_data$White) |
    is.na(complete_data$GerdAtIndex) |
    is.na(complete_data$h2r_max) |
    is.na(complete_data$ppi_max)
)

test_complete = complete_data %>% subset(complete_cases)
log(test_complete)

# complete_Data is in the format for our xgboost model
# the following keeps only the required info for kunzmann or hunt

# variables
df = test_complete[, c("ID", "CaseControl", "smoke_current", "smoke_former", "Gender", 
                       "GerdAtIndex", "White", "smoke_ever", "bmi", "ageatindex")]
df$k_age_bin = 
  relevel(cut(test_complete$ageatindex, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
df$k_bmi_bin = cut(test_complete$bmi, breaks=c(0, 25, 30, 35, 100))
df$h_age_bin = cut(test_complete$ageatindex, breaks=c(0, 50, 60, 70, 100))
df$h_bmi_bin = cut(test_complete$bmi, breaks=c(0, 30, 100))
df$h_smoke_any = pmax(
  test_complete$smoke_former,
  test_complete$smoke_current
)
df$k_ec = pmax(
  test_complete$h2r_max>0, 
  test_complete$ppi_max>0, 
  test_complete$GerdAtIndex
) > 0

test_complete$White = NULL
test_complete$smoke_ever = NULL
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
kscores = kunzmann_score(df)

kpred = as.numeric(kscores >= 8)
table(kpred, y)

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
hscores = hunt_score(df)


# screening guidelines
scores = df %>% mutate( 
   ACG2016     = Gender & GerdAtIndex & ( ((ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) ) > 1),
   ACG2022     = GerdAtIndex & ( (Gender + (ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) ) > 2),
   ACP2012     = Gender & GerdAtIndex & (ageatindex>=50) & ( ((bmi>30) + (smoke_ever) ) > 0),
   AGA2011     = ( GerdAtIndex + Gender + (ageatindex>=50) + (White) + (bmi>30) ) > 1,
   AGA_CPU2022 = ( GerdAtIndex + Gender + (ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) ) > 2,
   ASGE2019    = GerdAtIndex & ( (Gender + (ageatindex>=50) + (bmi>30) + (smoke_ever) ) > 0),
   BSG2013     = GerdAtIndex & ( (Gender + (ageatindex>=50) + (White) + (bmi>30) ) > 2),
   ESGE2017    = GerdAtIndex & ( (Gender + (ageatindex>=50) + (White) + (bmi>30) ) > 1)
)
scores %>% select(CaseControl, guidelines) %>% group_by(CaseControl) %>% summarise_all(sum)
guidelines = tail(colnames(scores), 8)

Hmisc::describe(df %>% select(-c(ID, bmi, ageatindex)))

n_risk_factors = df %>% 
  mutate(n_risk_factors=Gender+GerdAtIndex+(ageatindex>=50) + (White) + (bmi>30) + (smoke_ever)) %>%
  group_by(CaseControl, n_risk_factors) %>%
  summarize(n=n())

scores %>% select(guidelines) %>% colSums()

rbind(n_risk_factors, p_risk_factors)

# AUCs
y = df$CaseControl

proba = kscores
fg = proba[y==1]; bg = proba[y==0]
k_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)

proba = hscores
fg = proba[y==1]; bg = proba[y==0]
h_roc = PRROC::roc.curve(fg, bg ,curve=TRUE)

xgb_fit = readRDS("R_data/results/models/resample_nall.rds")$xgb_fit
newdata = xgb.DMatrix(as.matrix(test_complete[-c(1,2)]),
                      label=test_complete$CaseControl)
proba = predict(xgb_fit, newdata=newdata)
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

guide_roc$xlab = c(.1, .2, .3, .4, .5, .6, .8, .9)
guide_roc$ylab = c(.15, .2, .25, .3, .35, .4, .85, .9)

library(ggplot2)
filepath = paste0("R_code/hosea-project/figures/xgb_kunzmann_hunt_screening_roc.pdf")
g = ggplot(data=cdf, aes(x=fpr, y=recall, color=method)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  geom_point(data=guide_roc)
g = g +
  geom_segment(data=guide_roc, aes(x=fpr, xend=xlab, y=recall, yend=ylab)) + 
  geom_label(data=guide_roc, aes(label=label, x=xlab, y=ylab), size=3)
ggsave(filepath, g, width=8, height=7)
