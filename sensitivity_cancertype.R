setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(pdp)
library(HOSEA)
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
  ANY="1M_ANY.rds",
  EAC="1M_EAC.rds",
  EGJAC="1M_EGJAC.rds"
)
models = lapply(model_names, function(file){
  results = readRDS(paste0(dir_model, file))
  return(results)
})

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1.rds")
df = readRDS(file_path)
master = df$master
df = df$df

# =========================================================
# get outcomes
outcomes = master %>% select(ID, CancerType)
outcomes %<>% mutate(
  ANY=ifelse(CancerType=="", 0, 1),
  EAC=ifelse(CancerType=="EAC", 1, 0),
  EGJAC=ifelse(CancerType=="EGJAC", 1, 0)
)
# merge into df
df %<>% left_join(outcomes, by="ID")
outcome_names = c("ANY", "EAC", "EGJAC")

# =========================================================
# get ROC curves
probas = list()
rocs = list()
for(outcome in outcome_names){
  rocs[[outcome]] = list()
  probas[[outcome]] = list()
  # move outcome into CaseControl
  df %<>% mutate(CaseControl:=!!sym(outcome))
  # XGBoost
  for(name in names(models)){
    dff = df %>% filter(ID %in% models[[name]]$test_ids)
    set.seed(0)
    dff %<>% impute_srs(models[[name]]$quantiles)
    dff %<>% select(c(ID, CaseControl, models[[name]]$xgb_fit$feature_names))
    y = dff$CaseControl
    dff = xgb.DMatrix(as.matrix(dff%>% select(models[[name]]$xgb_fit$feature_names)),
                     label=dff$CaseControl)
    # get predicted risk and ROC curve
    proba = predict(models[[name]]$xgb_fit, newdata=dff)
    probas[[outcome]][[name]] = data.frame(proba=proba, y=y)
    fg = proba[y==1]; bg = proba[y==0]
    roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
    roc$curve = roc$curve[seq(1, nrow(roc$curve), by=1000), ]
    rocs[[outcome]][[name]] = roc
  }
}


# =========================================================
# postprocessing
aucs = data.frame(lapply(rocs, function(roc) sapply(roc, function(r) r$auc)))
xtable::xtable(aucs, digits=3)
# columns = outcome, rows = models
curves = list()
for(outcome in outcome_names){
  for(name in names(rocs[[outcome]])){
    curve = data.frame(rocs[[outcome]][[name]]$curve)
    colnames(curve) = c("fpr", "recall", "threshold")
    curve$outcome = outcome
    curve$model = name
    curves[[paste0(outcome, "_", name)]] = curve
  }
}
curves %<>% bind_rows()


# =========================================================
# plot
for(outcome in outcome_names){
  filepath = paste0(dir_figures, paste0("roc_", outcome, ".pdf"))
  tmp = curves %>% filter(outcome==!!outcome)
  g = ggplot(data=tmp, 
             aes(x=fpr, y=recall, color=model)) + 
    geom_line() +
    theme(aspect.ratio=1) +
    xlab("1 - Specificity") + ylab("Sensitivity") + 
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    ggtitle(paste0("Sensitivity analysis: Cancer type ", outcome))
  ggsave(filepath, g, width=9, height=7)
}

# =========================================================
# plot risk distribution
for(outcome in outcome_names){
  filepath = paste0(dir_figures, paste0("risk_test", outcome, ".pdf"))
  tmpprobas = lapply(names(probas[[outcome]]), function(nm) {
    tmpdf = probas[[outcome]][[nm]]
    tmpdf$model = nm
    tmpdf
  })
  tmp = tmpprobas %>% bind_rows()
  tmp %<>% filter(y==1)
  g = ggplot(data=tmp, aes(x=proba*100000, colour=model, fill=model)) + 
    geom_density(alpha=0.2) + xlab("Predicted risk (/100,000)") +
    ylab("Density") + scale_x_continuous(trans="log10") + 
    ggtitle(paste0("Risk distribution: Test ", outcome))
  ggsave(filepath, g, width=9, height=7)
}


for(trainnm in outcome_names){
  filepath = paste0(dir_figures, paste0("risk_train", trainnm, ".pdf"))
  tmpprobas = lapply(outcome_names, function(nm) {
    tmpdf = probas[[nm]][[trainnm]]
    tmpdf$outcome = nm
    tmpdf
  })
  tmp = tmpprobas %>% bind_rows()
  tmp %<>% filter(y==1)
  g = ggplot(data=tmp, aes(x=proba*100000, colour=outcome, fill=outcome)) + 
    geom_density(alpha=0.2) + xlab("Predicted risk (/100,000)") +
    ylab("Density") + scale_x_continuous(trans="log10") + 
    ggtitle(paste0("Risk distribution: Train ", trainnm))
  ggsave(filepath, g, width=9, height=7)
}




