setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(ggridges)
theme_set(theme_minimal())
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
file_path = paste0(dir_path, "5-1_test_merged.rds")
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


# =========================================================
# get predictions
df %<>% filter(id %in% models[["ANY"]]$test_ids)
pred = predict.HOSEA(df, 10)
pred_long = pred %>% tidyr::pivot_longer(c(ANY, EAC, EGJAC))
pred_long %<>% left_join(df %>% select(id, casecontrol), by="id")
pred_long$name2 = paste0(pred_long$name, " ", ifelse(pred_long$casecontrol==1, "case", "control"))

g = ggplot(pred_long, aes(y=name, x=value*100000)) + 
  geom_density_ridges(quantile_lines = TRUE) +
  scale_x_log10() +
  xlab("Predicted risk (/100,000)") + 
  ylab("Density")


filepath = paste0(dir_figures, paste0("risk_distribution.pdf"))
ggsave(filepath, g, width=6, height=4)


g = ggplot(pred_long, aes(y=name2, x=value*100000, fill=casecontrol)) + 
  geom_density_ridges(quantile_lines = TRUE) +
  scale_x_log10() +
  xlab("Predicted risk (/100,000)") + 
  ylab("Density") + 
  theme(legend.position="none")


filepath = paste0(dir_figures, paste0("risk_distribution_casecontrol.pdf"))
ggsave(filepath, g, width=6, height=4)

# =========================================================
# get ROC curves
probas = list()
rocs = list()
for(outcome in outcome_names){
  rocs[[outcome]] = list()
  probas[[outcome]] = list()
  # move outcome into casecontrol
  df %<>% mutate(casecontrol:=!!sym(outcome))
  # XGBoost
  for(name in names(models)){
    dff = df %>% filter(id %in% models[[name]]$test_ids)
    set.seed(0)
    dff %<>% impute_srs(models[[name]]$quantiles)
    dff %<>% select(c(id, casecontrol, models[[name]]$xgb_fit$feature_names))
    y = dff$casecontrol
    dff = xgb.DMatrix(as.matrix(dff%>% select(models[[name]]$xgb_fit$feature_names)),
                     label=dff$casecontrol)
    # get predicted risk and ROC curve
    proba = predict(models[[name]]$xgb_fit, newdata=dff)
    probas[[outcome]][[name]] = data.frame(proba=proba, y=y)
    fg = proba[y==1]; bg = proba[y==0]
    roc = PRROC::roc.curve(fg, bg ,curve=TRUE)
    roc$curve = roc$curve[seq(1, nrow(roc$curve), by=1000), ]
    
    proc = pROC::roc(controls=bg, cases=fg)
    roc$ci = pROC::ci(proc, of="auc")
    roc$display.ci = paste0(
      round(roc$au, 3), " [",
      round(roc$ci[1], 3), ",",
      round(roc$ci[3], 3), "]"
    )
    roc$display = round(roc$au, 3)
    
    rocs[[outcome]][[name]] = roc
  }
}


# =========================================================
# postprocessing
aucs = data.frame(lapply(rocs, function(roc) sapply(roc, function(r) r$display.ci)))
xtable::xtable(aucs)
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
  ggsave(filepath, g, width=6, height=4)
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
  ggsave(filepath, g, width=6, height=4)
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
  ggsave(filepath, g, width=6, height=4)
}




