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
dir_figures = "R_code/hosea-project/figures/risk_distributions/"
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
outcomes = master %>% select(id, cancertype)
df %<>% left_join(outcomes, by="id")

# =========================================================
# get all XGBoost risks

risk_df = lapply(names(models), function(mname){
  
  xgb_fit = models[[mname]]$xgb_fit
  quantiles = models[[mname]]$quantiles
  test_ids = models[[mname]]$test_ids
  
  # ensure correct column ordering for xgb model
  dff = df %>% select(c(id, casecontrol, xgb_fit$feature_names))
  y = dff$casecontrol
  # convert to xgb format
  xgb_df = xgb.DMatrix(as.matrix(dff %>% select(xgb_fit$feature_names)),
                   label=dff$casecontrol)
  # get predicted risk and ROC curve
  proba = predict(xgb_fit, newdata=xgb_df)
  
  risk = data.frame(dff$id, proba*100000)
  colnames(risk) = c("id", mname)
  
  return(risk)
  
})

names(risk_df) = names(models)

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


# =========================================================
# get Kunzmann risks

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

risk_df[["Kunzmann"]] = data.frame(id=dffk$id, Kunzmann=kscores)


# =========================================================
# get HUNT risks

hunt_score = function(dff){
  score = 3.6
  score = score * ifelse(dff$gender, 1.9, 1.)
  score = score * ifelse(dff$h_age_bin == "(50,60]", 2.1, 1.)
  score = score * ifelse(dff$h_age_bin == "(60,70]", 3.2, 1.)
  score = score * ifelse(dff$h_age_bin == "(70,100]", 3.1, 1.)
  score = score * ifelse(dff$h_bmi_bin == "(30,100]", 1.8, 1.)
  score = score * ifelse(dff$gerd, 3.7, 1.)
  score = score * ifelse(dff$h_smoke_any, 2.1, 1.)
  return(score)
}
hscores = hunt_score(dffk)

risk_df[["HUNT"]] = data.frame(id=dffk$id, HUNT=hscores)


# =========================================================
# join all risks
risks = purrr::reduce(risk_df, full_join, by="id")

risks %<>% left_join(df %>% select(id, cancertype, casecontrol), by="id")

risks %<>% mutate(cancertype=ifelse(cancertype=="", "Control", cancertype))


# =========================================================
# plot risks

for(mname in c(names(models), "Kunzmann", "HUNT")){
  g = ggplot(
    data=risks,
    mapping=aes_string(x=mname, color="cancertype", fill="cancertype")
  ) + 
    geom_density(alpha=0.5, bw=switch(mname, "HUNT"=0.1, "Kunzmann"=0.5, "nrd0")) + 
    ylab("Density") + xlab("Predicted risk") + ggtitle(mname)
  if(mname!="Kunzmann") g = g + scale_x_log10()
  filename = paste0(dir_figures, mname, ".pdf")
  g
  ggsave(filename, g, width=8, height=4)
}


# =========================================================
# threshold tables

metrics = function(proba, y, threshold=0.5){
  threshold = c(threshold)
  threshold = sort(threshold)
  yhat = lapply(threshold, function(tr) as.numeric(proba>tr))
  # metrics depending on threshold
  cols = c("tpr", "tnr", "ppv", "npv", "detection_prevalance", 
           "accuracy", "balanced_accuracy", "F1_score", "mcc",
           "jaccard_index")
  metrics = data.frame(matrix(NA, nrow=length(threshold), ncol=length(cols)))
  rownames(metrics) = threshold
  colnames(metrics) = cols
  for(i in seq(length(threshold))){
    yy = yhat[[i]]
    N = sum(!is.na(yy))
    p = sum(yy, na.rm=T)
    n = sum(1-yy, na.rm=T)
    tp = sum((yy*y), na.rm=T)
    tn = sum((1-yy)*(1-y), na.rm=T)
    fp = sum(yy*(1-y), na.rm=T)
    fn = sum((1-yy)*y, na.rm=T)
    metrics$ppv[i]                  = tp/p
    metrics$npv[i]                  = tn/n
    metrics$tpr[i]                  = tp/(tp+fn)
    metrics$tnr[i]                  = tn/(fp+tn)
    metrics$detection_prevalance[i] = p/N
    metrics$F1_score[i]             = 2 *tp/(2*tp+fp+fn)
    metrics$accuracy[i]             = (tp+tn)/N
    metrics$balanced_accuracy[i]    = (metrics$tpr[i] + metrics$tnr[i]) / 2
    metrics$mcc[i]                  = (tp*tn-fp*fn) / sqrt(p*n*(tp+fn)*(tn+fp))
    metrics$jaccard_index[i]        = (tp+tn) / (2*N-tp-tn)
    
  }
  
  return(metrics)
}

# Kunzmann
thresholds = seq(0, 15, 0.5)
tr_df = metrics(risks$Kunzmann, risks$casecontrol, thresholds)
tr_df = tr_df %>% select(one_of("tpr", "ppv", "detection_prevalance"))
# rownames(tr_df) = format(round(as.numeric(rownames(tr_df)), 5)*100000, digits=5)
tr_df = tr_df*100
cat(print(xtable::xtable(tr_df)), 
    file=paste0(dir_figures, "Kunzmann_threshold.tex"))

# HUNT
thresholds = seq(0, 320, 10)
tr_df = metrics(risks$HUNT, risks$casecontrol, thresholds)
tr_df = tr_df %>% select(one_of("tpr", "ppv", "detection_prevalance"))
# rownames(tr_df) = format(round(as.numeric(rownames(tr_df)), 5)*100000, digits=5)
tr_df = tr_df*100
cat(print(xtable::xtable(tr_df)), 
    file=paste0(dir_figures, "HUNT_threshold.tex"))

