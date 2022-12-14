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
# get a sense of how many are use multiple times
master %>% group_by(casecontrol, sset) %>%
  summarise(n=n(), n_unique=n_distinct(id)) %>%
  mutate(repeated=n-n_unique)

# NB: some repeated IDs in the new data?
# # A tibble: 8 × 5
# # Groups:   casecontrol [2]
# casecontrol sset        n n_unique repeated
# <dbl> <chr>   <int>    <int>    <int>
# 1           0 all   6609865  6588877    20988
# 2           0 test  2564221  2564221        0
# 3           0 train 5128443  5128443        0
# 4           0 valid 2564221  2564221        0
# 5           1 all      2094     2082       12
# 6           1 test     2848     2848        0
# 7           1 train    5697     5697        0
# 8           1 valid    2848     2848        0

master %>% group_by(casecontrol) %>%
  summarise(n=n(), n_unique=n_distinct(id)) %>%
  mutate(repeated=n-n_unique)
# most new control are repeated
# # A tibble: 2 × 4
#   casecontrol         n n_unique repeated
# <dbl>    <int>    <int>    <int>
# 1           0 16866750 11248566  5618184
# 2           1    13487    13467       20

master %>% group_by(sset) %>%
  summarise(n=n(), n_unique=n_distinct(id)) %>%
  mutate(repeated=n-n_unique)
# # A tibble: 4 × 4
#   sset        n n_unique repeated
# <chr>   <int>    <int>    <int>
# 1 all   6611959  6590959    21000
# 2 test  2567069  2567069        0
# 3 train 5134140  5134140        0
# 4 valid 2567069  2567069        0

# check if any case->case, case->control, control->case
master %>% arrange(cohort, sset) %>% group_by(id) %>%
  summarise(n=n(), n_cases=sum(casecontrol==1)) %>%
  group_by(n, n_cases) %>% summarize(count=n())
# A few are repeated multiple times!?
# A tibble: 13 × 3
# Groups:   n [8]
#       n     n_cases   count
# <int>   <int>   <int>
# 1     1       0 5646552
# 2     1       1   11490
# 3     2       0 5593352
# 4     2       1    1965
# 5     2       2       8
# 6     4       0    2832
# 7     4       4       2
# 8     5       0    3779
# 9     5       4       2
# 10     9       0       3
# 11    10       0      12
# 12    16       0      16
# 13    17       0      53

# - exclude prior cases from being cases or controls
# - exclude controls that were in training or validation sets
ids_to_exclude = c(
  og_train$id,
  og_valid$id
)

icd10_select = icd10_raw %>% filter(!(id %in% ids_to_exclude))
saveRDS(icd10_select, paste0(dir_icd10, y0, "-", y1, "filtered.rds"))

icd10_select$casecontrol %>% table()
# 0       1 
# 2402295     638 

imputer = HOSEA::load_imputer()

icd10_imputed = HOSEA::impute(imputer, icd10_select, 1)
saveRDS(icd10_imputed, paste0(dir_out, y0, "-", y1, ".rds"))


icd10_imputed %>% group_by(id) %>% summarise(count=n()) %>% group_by(count) %>% summarise(n=n())
# # A tibble: 4 × 2
#   count       n
# <int>   <int>
# 1     1 2387412
# 2     4    3765
# 3     9       5
# 4    16      26