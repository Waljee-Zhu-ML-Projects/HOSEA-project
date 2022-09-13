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
setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = "./R_code/hosea-project/figures/calibration/"
dir_tables = "./R_code/hosea-project/tables/calibration/"
imputed_data = "test_mice_any.rds"
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------




# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------




# ==============================================================================
# PARAMETERS
outcome = "EAC"
complete = F # F: uses imputed data, T: only use complete records wrt other methods
representative = F # F: uses everything, T: downsamples males so get a more representative sample 
seed = 0
# ------------------------------------------------------------------------------

# to get everything
for(outcome in c("ANY", "EAC", "EGJAC")){
for(complete in c(T, F)){
for(representative in c(T, F)){

      
# ==============================================================================
# PREPARE DATA
# get complete cases for Kunzmann, HUNT and guidelines
if(complete){
  ids = complete_for_comparison(raw_df$df)
}else{
  ids = imputed_df %>% pull(id)
}
# working df
imputed_wdf = imputed_df %>% filter(id %in% ids)
imputed_wdf %<>% patch_outcome(raw_df$master, outcome=outcome)
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
proba = predict.HOSEA(imputed_wdf, imputer=NULL) %>% 
  select(id, !!outcome) %>% rename(HOSEA=!!outcome)
y = imputed_wdf %>% pull(casecontrol)
# ------------------------------------------------------------------------------




# ==============================================================================
# CALIBRATION PLOT
calibration = calibration_curve(proba$HOSEA, y)
hl = hosmer_lemeshow(
  calibration$mid,
  calibration$prop_cases,
  calibration$N,
  calibration$N_cases
)

for(log in c(T, F)){
  
filepath = paste0(dir_figures, outcome, "_", 
                  ifelse(complete, "complete", "imputed"), 
                  ifelse(representative, "_representative", ""),
                  ifelse(log, "_log", ""),
                  ".pdf")
main = paste0("Cancer type: ", outcome, "\n",
               "Dataset: test, ", ifelse(complete, "complete", "imputed"), 
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
  g = g + xlim(0, 400) + ylim(0, 400)
}
g
ggsave(filepath, g, width=5, height=6)

}
# ------------------------------------------------------------------------------




# ==============================================================================
# THRESHOLD TABLE
ttable = classification_metrics(proba$HOSEA, y)

filepath = paste0(dir_tables, "thresholds_", outcome, "_", 
                  ifelse(complete, "complete", "imputed"), 
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

      
      
  
