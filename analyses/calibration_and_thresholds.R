# ==============================================================================
# CALILBRATION CURVE AND THRESHOLD PERFORMANCE TABLE
# Author: Simon Fontaine (simfont@umich.edu)
# ------------------------------------------------------------------------------




# ==============================================================================
# REQUIRED PACKAGES
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(HOSEA)
theme_set(theme_minimal())
# ------------------------------------------------------------------------------




# ==============================================================================
# PATHS
imputation = "srs"
setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/calibration/")
dir_tables = paste0("./R_code/hosea-project/tables/", imputation, "/calibration/")
imputed_data = paste0("5-1test_", imputation, "_any.rds")
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------


# ==============================================================================
# MODELS
models = load_models(
  files_meta=list(
    ANY=paste0("xgb_", imputation, "_any.meta"), 
    EAC=paste0("xgb_", imputation, "_eac.meta"), 
    EGJAC=paste0("xgb_", imputation, "_egjac.meta")
  ),
  files_models=list(
    ANY=paste0("xgb_", imputation, "_any.model"), 
    EAC=paste0("xgb_", imputation, "_eac.model"), 
    EGJAC=paste0("xgb_", imputation, "_egjac.model")
  )
)
# ------------------------------------------------------------------------------



# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------




# ==============================================================================
# PARAMETERS
seed = 0
outcome = "ANY"
missing_which = "all"
representative = F
# ------------------------------------------------------------------------------

# to get everything
for(outcome in c("ANY", "EAC", "EGJAC")){
for(missing_which in c("all", "complete", "incomplete")){
for(representative in c(T, F)){

      
# ==============================================================================
# PREPARE DATA
# get complete cases for Kunzmann, HUNT and guidelines
if(missing_which == "complete"){
  ids = complete_for_comparison(raw_df$df)
}
if(missing_which == "all"){
  ids = imputed_df %>% pull(id)
}
if(missing_which == "incomplete"){
  complete_ids = complete_for_comparison(raw_df$df)
  ids = imputed_df %>% filter(!(id %in% complete_ids)) %>% pull(id)
}
# working df
imputed_wdf = imputed_df %>% filter(id %in% ids)
imputed_wdf %<>% patch_outcome(master=raw_df$master, outcome=outcome)
# representative sample
if(representative){
  set.seed(seed)
  imputed_wdf %<>% representative_sample()
}

n_cases = imputed_wdf$casecontrol %>% sum()
n_patients = imputed_wdf %>% nrow()
# ------------------------------------------------------------------------------




# ==============================================================================
# GET SCORES
models = load_models(
    files_meta=list(
      ANY=paste0("xgb_", imputation, "_any.meta"), 
      EAC=paste0("xgb_", imputation, "_eac.meta"), 
      EGJAC=paste0("xgb_", imputation, "_egjac.meta")
      ),
    files_models=list(
      ANY=paste0("xgb_", imputation, "_any.model"), 
      EAC=paste0("xgb_", imputation, "_eac.model"), 
      EGJAC=paste0("xgb_", imputation, "_egjac.model")
      )
)
proba = predict.HOSEA(imputed_wdf, imputer=NULL, models=models) %>% 
  select(id, !!outcome) %>% rename(HOSEA=!!outcome)
y = imputed_wdf %>% pull(casecontrol)
# ------------------------------------------------------------------------------




# ==============================================================================
# CALIBRATION PLOT
calibration = calibration_curve(proba$HOSEA, y, 50)
hl = hosmer_lemeshow(
  calibration$mid,
  calibration$prop_cases,
  calibration$N,
  calibration$N_cases
)

for(log in c(T, F)){
  
filepath = paste0(dir_figures, outcome, "_", 
                  missing_which, 
                  ifelse(representative, "_representative", ""),
                  ifelse(log, "_log", ""),
                  ".pdf")
main = paste0("Cancer type: ", outcome, "\n",
               "Dataset: test, ", missing_which, 
               ifelse(representative, ", representative", ""), "\n",
               "Cases: ", n_cases, "/", n_patients, "\n",
               "HL: ", hl)
g = ggplot(data=calibration, aes(x=mid*100000, y=prop_cases*100000)) + theme(aspect.ratio=1) + 
  geom_point()  +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ylab("Observed (/100,000)") + xlab("Predicted (/100,000)") + 
  ggtitle(main)
if(log){
  m = max(max(calibration$mid), max(calibration$prop_cases))*100000
  g = g + scale_x_log10(limits=c(1, m)) + scale_y_log10(limits=c(1, m))
}else{
  g = g + xlim(0, ifelse(outcome=="EGJAC", 200, 400)) + ylim(0, ifelse(outcome=="EGJAC", 200, 400))
}
g
ggsave(filepath, g, width=5, height=6, bg="white")
ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=5, height=6, bg="white")

}
# ------------------------------------------------------------------------------




# ==============================================================================
# THRESHOLD TABLE
ttable = classification_metrics(proba$HOSEA, y)

filepath = paste0(dir_tables, "thresholds_", outcome, "_", 
                  missing_which, 
                  ifelse(representative, "_representative", ""),
                  ".tex")

textable = ttable %>% 
  select(threshold, tpr, ppv, det_prev) %>%
  mutate(tpr=100*tpr, ppv=100*ppv, det_prev=100*det_prev) %>%
  rename(Threshold=threshold, TPR=tpr, PPV=ppv, DetPrevalence=det_prev) %>% 
  mutate(Threshold=as.character(Threshold))

nrows = nrow(textable)
every5 = seq(5, nrows-1, 5)

textable = xtable::xtable(textable, align="llccc", digits=2)
textable %<>% print(
  table.placement="ht",
  booktabs=T,
  include.rownames=F,
  add.to.row=list(pos=as.list(every5), command=rep("\\addlinespace\n", length(every5)))
)
cat(textable, file=filepath)

# ------------------------------------------------------------------------------


}}}

      
      
  
