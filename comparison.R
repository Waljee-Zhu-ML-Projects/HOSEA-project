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
dir_model = "R_data/results/models/"

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
# subset to outcomes

outcomes = master %>% select(id, casecontrol, cancertype)
outcomes %<>% mutate(
  ANY=casecontrol,
  EAC=as.integer(cancertype=="EAC"),
  EGJAC=as.integer(cancertype=="EGJAC")
)


# =========================================================
# get all AUCs
aucs = data.frame(matrix(NA, nrow=3, ncol=3))
colnames(aucs) = names(model_names)
rownames(aucs) = names(model_names)

# columns are training (ie model name)
# rows are testing (ie outcome)

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

for(model in names(model_names)){
  for(outcome in names(model_names)){
    xgb_fit = models[[outcome]]$xgb_fit
    quantiles = models[[outcome]]$quantiles
    test_ids = models[[outcome]]$test_ids
    
    outcomes_ = outcomes%>%select(id, !!outcome)
    
    dff = df %>% filter(id %in% test_ids)
    dff %<>% 
      left_join(outcomes_, by="id") %>%
      select(-casecontrol) %>% 
      rename(casecontrol=!!outcome)
    
    set.seed(0)
    dff = impute_srs(dff, quantiles)
    
    roc = get_roc(dff)
    
    aucs[outcome, model] = roc$display.ci
  }
}


# =========================================================
# ROC curves
outcome = "EGJAC"

xgb_fit = models[[outcome]]$xgb_fit
quantiles = models[[outcome]]$quantiles
test_ids = models[[outcome]]$test_ids

outcomes_ = outcomes%>%select(id, !!outcome)

dff = df %>% filter(id %in% test_ids)
dff %<>% 
  left_join(outcomes_, by="id") %>%
  select(-casecontrol) %>% 
  rename(casecontrol=!!outcome)

# =========================================================
# find complete cases for other methods

dff$white = !(dff$asian | dff$black | dff$hawaiianpacific | dff$indianalaskan)
dff$smoke_ever = dff$smoke_current | dff$smoke_former

complete_cases = !(
  is.na(dff$gender) |
    is.na(dff$age) |
    is.na(dff$smoke_former) |
    is.na(dff$smoke_current) |
    is.na(dff$bmi) |
    is.na(dff$white) | # will be NA iff same for other columns
    is.na(dff$gerd) |
    (is.na(dff$h2r_max) &
    is.na(dff$ppi_max))
)

test_complete = dff %>% subset(complete_cases)
test_complete$casecontrol %>% table()

# test_complete is in the format for our xgboost model (w/ extra variables)
# the following keeps only the required info for kunzmann or hunt

# variables
dffk = test_complete %>% select(c(
  "id", "casecontrol", "smoke_current", "smoke_former", "gender", 
  "gerd", "white", "smoke_ever", "bmi", "age")
)

# check
dffk %>% filter(casecontrol==1) %>% select(c(
  "id", "casecontrol", "smoke_current", "smoke_former", "gender", 
  "gerd", "white", "smoke_ever", "bmi", "age"
  )) %>% is.na() %>% colSums()

# compute various bins
dffk$k_age_bin = 
  relevel(cut(test_complete$age, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
dffk$k_bmi_bin = cut(test_complete$bmi, breaks=c(0, 25, 30, 35, 100))
dffk$h_age_bin = cut(test_complete$age, breaks=c(0, 50, 60, 70, 100))
dffk$h_bmi_bin = cut(test_complete$bmi, breaks=c(0, 30, 100))
dffk$h_smoke_any = pmax(
  test_complete$smoke_former,
  test_complete$smoke_current
)
dffk$k_ec = pmax(
  test_complete$h2r_max>0, 
  test_complete$ppi_max>0, 
  test_complete$gerd, na.rm=T
) > 0

# drop extra columns just in case
test_complete$white = NULL
test_complete$smoke_ever = NULL
# Kunzmann score
kunzmann_score = function(dff){
  score = 0
  score = score + (dff$k_age_bin == "(55,60]") * 1.5
  score = score + (dff$k_age_bin == "(60,65]") * 2.5
  score = score + (dff$k_age_bin == "(65,100]") * 3.5
  score = score + dff$gender * 4
  score = score + (dff$k_bmi_bin == "(25,30]") * 1
  score = score + (dff$k_bmi_bin == "(30,35]") * 1.5
  score = score + (dff$k_bmi_bin == "(35,100]") * 2.5
  score = score + dff$smoke_former * 2.5
  score = score + dff$smoke_current * 3.5
  score = score + dff$k_ec * 1.5
  return(score)
}
kscores = kunzmann_score(dffk)

hunt_score = function(dff){
  score = 3.6
  score = score * ifelse(dff$gender, 1.9, 1.)
  score = score * ifelse(dff$h_age_bin == "(50,60]", 2.1, 1.)
  score = score * ifelse(dff$h_age_bin == "(60,70]", 3.2, 1.)
  score = score * ifelse(dff$h_age_bin == "(70,100]", 3.1, 1.)
  score = score * ifelse(dff$h_bmi_bin == "(30,100]", 1.8, 1.)
  score = score * ifelse(dff$gerd, 3.7, 1.)
  score = score * ifelse(dff$h_smoke_any, 2.1, 1.)
  return(score/100000)
}
hscores = hunt_score(dffk)


# screening guidelines
scores = dffk %>% mutate( 
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
# Hmisc::describe(dffk %>% select(-c(id, bmi, age)))
# 
# n_risk_factors = dffk %>% 
#   mutate(n_risk_factors=gender+gerd+(age>=50) + (white) + (bmi>30) + (smoke_ever)) %>%
#   group_by(casecontrol, n_risk_factors) %>%
#   summarize(n=n())
# 
# scores %>% select(guidelines) %>% colSums()
# 
# rbind(n_risk_factors, p_risk_factors)

add_ci = function(roc){
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

# AUCs
y = dffk$casecontrol

proba = kscores
fg = proba[y==1]; bg = proba[y==0]
k_roc = PRROC::roc.curve(fg, bg ,curve=TRUE) %>% add_ci()

proba = hscores
fg = proba[y==1]; bg = proba[y==0]
h_roc = PRROC::roc.curve(fg, bg ,curve=TRUE) %>% add_ci()

xgb_dff = xgb.DMatrix(as.matrix(test_complete %>% select(xgb_fit$feature_names)),
                     label=test_complete$casecontrol)
proba = predict(xgb_fit, newdata=xgb_dff)
fg = proba[y==1]; bg = proba[y==0]
x_roc = PRROC::roc.curve(fg, bg ,curve=TRUE) %>% add_ci()


curves = data.frame(x_roc$curve[, 1:2], paste0("HOSEA (", x_roc$display.ci, ")"))
colnames(curves) = c("fpr", "recall", "Method")
cdf = curves
curves = data.frame(k_roc$curve[, 1:2], paste0("Kunzmann (", k_roc$display.ci, ")"))
colnames(curves) = c("fpr", "recall", "Method")
cdf = rbind(cdf, curves)
curves = data.frame(h_roc$curve[, 1:2], paste0("HUNT (", h_roc$display.ci, ")"))
colnames(curves) = c("fpr", "recall", "Method")
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
guide_roc$Method = "Guideline"
guide_roc$label = rownames(guide_roc)

guide_roc %<>% arrange(fpr)

guide_roc$xlab = c(.3, .37, .44, .51, .58, .65, .8, .9)
guide_roc$ylab = c(.15, .2, .25, .3, .35, .4, .85, .9)

filepath = paste0("R_code/hosea-project/figures/comparison_", outcome, ".pdf")
g = ggplot(data=cdf, aes(x=fpr, y=recall, color=Method)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  geom_point(data=guide_roc) + 
  ggtitle(paste0("Cancer type: ", outcome))
g = g +
  geom_segment(data=guide_roc, aes(x=fpr, xend=xlab, y=recall, yend=ylab)) + 
  geom_label(data=guide_roc, aes(label=label, x=xlab, y=ylab), size=3, show.legend=FALSE)
g
ggsave(filepath, g, width=8, height=5)
