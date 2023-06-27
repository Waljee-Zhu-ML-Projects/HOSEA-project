# ==============================================================================
# INSPECTING FEATURE DISTRIBUTION AND IMPUTATION
# Author: Simon Fontaine (simfont@umich.edu)
# ------------------------------------------------------------------------------




# ==============================================================================
# REQUIRED PACKAGES
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(HOSEA)
library(rigr)
theme_set(theme_minimal())
# ------------------------------------------------------------------------------




# ==============================================================================
# PATHS
setwd('/nfs/turbo/umms-awaljee-secure/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = "./R_code/hosea-project/figures/feature_comparison/"
dir_tables = "./R_code/hosea-project/tables/srs/feature_comparison/"
raw_data = "5-1_merged.rds"
imputed_test_data = "5-1test_srs_any.rds"
imputed_train_data = "5-1train_srs_any.rds"
imputed_valid_data = "5-1valid_srs_any.rds"

# Binary variables
bin_vars = c(
  "gender", "casecontrol",
  "asian", "black", "hawaiianpacific", "indianalaskan",
  "smoke_current", "smoke_former",
  HOSEA:::charlson_vars
)
# ------------------------------------------------------------------------------




# ==============================================================================
# READ IN DATA
raw_df = readRDS(paste0(dir_raw_data, raw_data))
train_ids = readRDS(paste0(dir_imputed_data, imputed_train_data))$id
valid_ids = readRDS(paste0(dir_imputed_data, imputed_valid_data))$id
test_ids = readRDS(paste0(dir_imputed_data, imputed_test_data))$id
raw_df$df %<>%
  mutate(
    subset=ifelse(
      id %in% train_ids, 
      "training",
      ifelse(id %in% valid_ids, "validation", "testing")
    )
  )
# ------------------------------------------------------------------------------







# ==============================================================================
# COMPUTE STATISTICS
# df = raw_df$df %>% sample_n(10000)
feat_dist = lapply(
  df %>% select(-id, -subset, -casecontrol) %>% colnames,
  function(cname){
    bin_var = cname %in% bin_vars
    xtrain_controls = df %>% filter(subset=="training", casecontrol==0) %>% pull(!!cname)
    xvalid_controls = df %>% filter(subset=="validation", casecontrol==0) %>% pull(!!cname)
    xtest_controls = df %>% filter(subset=="testing", casecontrol==0) %>% pull(!!cname)
    xdev_controls = df %>% filter(subset!="testing", casecontrol==0) %>% pull(!!cname)
    xtrain_cases = df %>% filter(subset=="training", casecontrol==1) %>% pull(!!cname)
    xvalid_cases = df %>% filter(subset=="validation", casecontrol==1) %>% pull(!!cname)
    xtest_cases = df %>% filter(subset=="testing", casecontrol==1) %>% pull(!!cname)
    xdev_cases = df %>% filter(subset!="testing", casecontrol==1) %>% pull(!!cname)
    out = data.frame(
      feature=cname,
      bin_var=bin_var
      ,
      
      ntrain_controls=sum(!is.na(xtrain_controls)),
      natrain_controls=sum(is.na(xtrain_controls)),
      mtrain_controls=mean(xtrain_controls, na.rm=T),
      strain_controls=sd(xtrain_controls, na.rm=T)
      ,
      nvalid_controls=sum(!is.na(xvalid_controls)),
      navalid_controls=sum(is.na(xvalid_controls)),
      mvalid_controls=mean(xvalid_controls, na.rm=T),
      svalid_controls=sd(xvalid_controls, na.rm=T)
      ,
      ntest_controls=sum(!is.na(xtest_controls)),
      natest_controls=sum(is.na(xtest_controls)),
      mtest_controls=mean(xtest_controls, na.rm=T),
      stest_controls=sd(xtest_controls, na.rm=T)
      ,
      ndev_controls=sum(!is.na(xdev_controls)),
      nadev_controls=sum(is.na(xdev_controls)),
      mdev_controls=mean(xdev_controls, na.rm=T),
      sdev_controls=sd(xdev_controls, na.rm=T)
      ,
      
      ntrain_cases=sum(!is.na(xtrain_cases)),
      natrain_cases=sum(is.na(xtrain_cases)),
      mtrain_cases=mean(xtrain_cases, na.rm=T),
      strain_cases=sd(xtrain_cases, na.rm=T)
      ,
      nvalid_cases=sum(!is.na(xvalid_cases)),
      navalid_cases=sum(is.na(xvalid_cases)),
      mvalid_cases=mean(xvalid_cases, na.rm=T),
      svalid_cases=sd(xvalid_cases, na.rm=T)
      ,
      ntest_cases=sum(!is.na(xtest_cases)),
      natest_cases=sum(is.na(xtest_cases)),
      mtest_cases=mean(xtest_cases, na.rm=T),
      stest_cases=sd(xtest_cases, na.rm=T)
      ,
      ndev_cases=sum(!is.na(xdev_cases)),
      nadev_cases=sum(is.na(xdev_cases)),
      mdev_cases=mean(xdev_cases, na.rm=T),
      sdev_cases=sd(xdev_cases, na.rm=T)
    )
    out
  }
) %>% bind_rows()
write.csv(feat_dist, paste0(dir_tables, "subset_stats.csv"))
# ------------------------------------------------------------------------------







# ==============================================================================
# PERFORM TESTS
vttest = Vectorize(ttesti, vectorize.args=c("obs", "mean", "sd", "obs2", "mean2", "sd2"))
vptest = Vectorize(proptesti, vectorize.args=c("x1", "n1", "x2", "n2"))
feat_dist %>% mutate(
  pval_controls_train_v_test=tryCatch(ifelse(
    bin_var,
    {
      vptest(
        x1=round(ntrain_controls * mtrain_controls),
        n1=ntrain_controls,
        x2=round(ntest_controls * mtest_controls),
        n2=ntest_controls
      )$p
    },
    {
      vttest(
        obs=ntrain_controls,
        mean=mtrain_controls,
        sd=strain_controls,
        obs2=ntest_controls,
        mean2=mtest_controls,
        sd2=stest_controls
      )$p
    }
  ), error=function(err) NA)
)

# ------------------------------------------------------------------------------

