setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(HOSEA)
library(dplyr)
library(magrittr)

y0=3; y1=1

dir_icd10 = "./R_data/processed_records/icd10cohort/" 
dir_og = "./R_data/imputed_records/"
dir_out = "./R_data/imputed_records/icd10cohort/" 

icd10_raw = readRDS(paste0(dir_icd10, y0, "-", y1, ".rds"))$df
og_train = readRDS(paste0(dir_og, "5-1train_srs_any.rds"))
og_valid = readRDS(paste0(dir_og, "5-1valid_srs_any.rds"))
og_test = readRDS(paste0(dir_og, "5-1test_srs_any.rds"))

# figure out which observations to keep
master = bind_rows(
  og_train %>% distinct(id, .keep_all=T) %>% select(id, casecontrol) %>% mutate(sset="train", cohort="og"),
  og_valid %>% select(id, casecontrol) %>% mutate(sset="valid", cohort="og"),
  og_test %>% select(id, casecontrol) %>% mutate(sset="test", cohort="og"),
  icd10_raw %>% select(id, casecontrol) %>% mutate(sset="all", cohort="icd10")
)


sample_df = HOSEA:::load_sas("./unzipped_data/icd10cohort/sampleext.rds", "sample")

# - exclude prior cases from being cases or controls
# - exclude controls that were in training or validation sets
ids_to_exclude = c(
  og_train$id,
  og_valid$id
)

icd10_select = icd10_raw %>% filter(!(id %in% ids_to_exclude))
saveRDS(icd10_select, paste0(dir_icd10, y0, "-", y1, "filtered.rds"))


imputer = HOSEA::load_imputer()

icd10_imputed = HOSEA::impute(imputer, icd10_select, 1)
saveRDS(icd10_imputed, paste0(dir_out, y0, "-", y1, ".rds"))

