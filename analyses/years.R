# ==============================================================================
# COMPARISON OF PERFORMANCE BY NUMBER OF YEARS OF DATA
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
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/years/")
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
# PARAMETERS
outcome = "ANY"
missing_which = "all"
representative = F # F: uses everything, T: downsamples males so get a more representative sample 
seed = 0

imputed_data = list(
  "[5-1]"=paste0("5-1test_", imputation, "_any.rds"),
  "[5-3]"=paste0("5-3test_", imputation, "_any.rds"),
  "[4-2]"=paste0("4-2test_", imputation, "_any.rds"),
  "[3-1]"=paste0("3-1test_", imputation, "_any.rds")
)

raw_data = list(
  "[5-1]"="5-1_merged.rds",
  "[5-3]"="5-3_merged.rds",
  "[4-2]"="4-2_merged.rds",
  "[3-1]"="3-1_merged.rds"
)

subsets = list(
  "2y" = c("[5-1]", "[5-3]", "[4-2]", "[3-1]")
)
# ------------------------------------------------------------------------------

for(outcome in c("ANY", "EAC", "EGJAC")){
model = models[[outcome]]

for(nsubset in names(subsets)){
dfs = subsets[[nsubset]]
rocs = lapply(dfs, function(ndf){
  cat(ndf, "\n")
  
  # ==============================================================================
  # READ IN DATA
  imputed_df = readRDS(paste0(dir_imputed_data, imputed_data[[ndf]]))
  raw_df = readRDS(paste0(dir_raw_data, raw_data[[ndf]]))
  # ------------------------------------------------------------------------------
  
  
  # ==============================================================================
  # PREPARE DATA
  # get complete cases for Kunzmann, HUNT and guidelines
  imputed_wdf = prepare_test_set(
    df=imputed_df,
    raw_df=raw_df$df,
    master=raw_df$master,
    missing_which=missing_which,
    outcome=outcome,
    representative=representative,
    seed=seed
  )
  # ------------------------------------------------------------------------------
  
  
  # ==============================================================================
  # GET SCORES
  scores = predict.HOSEA(imputed_wdf, imputer=NULL, models=models) %>% 
    select(id, !!outcome) %>% 
    rename(HOSEA=!!outcome)
  scores %<>% 
    left_join(imputed_wdf %>% select(id, casecontrol), by="id") %>% 
    rename(y=casecontrol)
  # ------------------------------------------------------------------------------
  
  
  
  # ==============================================================================
  # GET ROC
  roc_ = roc(scores)
  # ------------------------------------------------------------------------------
  
  return(roc_$HOSEA)
})
  

# ==============================================================================
# PRODUCE PLOT
names(rocs) = dfs

# build dfs
df_curves = sapply(names(rocs), function(name){
  curve_df = data.frame(rocs[[name]]$curve[, 1:2], paste0(name, " (", rocs[[name]]$display.ci, ")"))
  colnames(curve_df) = c("fpr", "recall", "Years")
  return(curve_df)
}, simplify=F) %>% bind_rows()

filepath = paste0(dir_figures, nsubset, "_", outcome, "_", 
                  missing_which, 
                  ifelse(representative, "_representative", ""),
                  ".pdf")
gtitle = paste0("Cancer type: ", outcome, "\n",
                "Dataset: test, ", missing_which, 
                ifelse(representative, ", representative", ""))
g = ggplot(data=df_curves, aes(x=fpr, y=recall, color=Years)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  ggtitle(gtitle)
g
ggsave(filepath, g, width=8, height=6, bg="white")
ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=8, height=6, bg="white")
# ------------------------------------------------------------------------------

  
}}
