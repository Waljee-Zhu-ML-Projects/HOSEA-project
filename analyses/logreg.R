# ==============================================================================
# LOGISTIC REGRESSION
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
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/labs/")
imputed_data = paste0("5-1test_", imputation, "_any.rds")
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------






# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------







# ==============================================================================
# PREPARE DATA
# df = imputed_df
set.seed(0)
# subsample for SHAP/PDPs
# df = imputed_df %>% sample_n(1000000)
# subsample controls
df = imputed_df %>% filter(casecontrol==1) %>% bind_rows(imputed_df%>%filter(casecontrol==0) %>% sample_n(100000))
# compute some variables
df %<>% mutate(
  white=(black+asian+hawaiianpacific+indianalaskan)==0,
  obesity=bmi>30,
  smoke_any=pmax(smoke_former, smoke_current),
  anion_gap_mean=na_mean + k_mean - co2_mean - chlor_mean
  )
df %<>% left_join(
  raw_df$df %>% select(id, ppi_max) %>%
    mutate(
      ppi=factor(
        ifelse(is.na(ppi_max), "NONE", ifelse(ppi_max<30, "LOW", "HIGH")), 
        levels=c("NONE", "LOW", "HIGH"))
      ),
  by="id"
)
# ------------------------------------------------------------------------------






# ==============================================================================
# MODELS
models = list(
  "age"=c(age=T),
  "sex"=c(gender=F),
  "gerd"=c(gerd=F),
  "white"=c(white=F),
  "obesity"=c(obesity=F),
  "smoking"=c(smoke_any=F),
  "known_predictors"=c(age=T, gender=F, gerd=F, white=F, obesity=F, smoke_any=F),
  "anion_gap"=c(anion_gap_mean=T),
  "known_predictors_anion_gap"=c(age=T, gender=F, gerd=F, white=F, obesity=F, smoke_any=F, anion_gap_mean=T),
  "ppi"=c(ppi=F),
  "ppi_anion_gap"=c(ppi=F, anion_gap_mean=T)
)
# ------------------------------------------------------------------------------






# ==============================================================================
# FORMULA BUILDER
to_formula = function(x){
  vnames = names(x)
  vterms = sapply(vnames, function(vname) ifelse(x[vname], paste0("s(", vname, ")"), vname))
  f = paste(vterms, collapse=" + ")
  f = paste0("casecontrol ~ ", f)
  return(f)
}
# ------------------------------------------------------------------------------






# ==============================================================================
# FIT
mname = "known_predictors_anion_gap"
model = models[[mname]]
fit = mgcv::gam(formula=to_formula(model) %>% formula, family=binomial, data=df)
summary_fit = summary(fit)
plot(fit)
coef(fit)
# ------------------------------------------------------------------------------



























