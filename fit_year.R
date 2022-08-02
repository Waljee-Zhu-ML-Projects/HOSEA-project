setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/evaluation_split.R')
library(HOSEA)
library(dplyr)
library(magrittr)

# import data
complete_data = readRDS('R_data/processed_records/5-1_merged.rds')
master = complete_data$master
complete_data = complete_data$df

# drop columns
complete_data %<>% select(-c(starts_with("chol"), starts_with("rbc"), starts_with("hgb")))

# add in year
master$indexdate = master$end + 365
master %<>% left_join(complete_data %>% select(id, age), by="id")
year_end = 19055
master %<>% mutate(indexyear = (indexdate - year_end)/365 + 2019)
hist(master$indexyear)
master %<>% mutate(birthyear = indexyear - age)
hist(master$birthyear)
master %<>% mutate(birth35 = birthyear >= 1935)

complete_data %<>% left_join(master%>%select(id,birthyear,birth35), by="id")


# details
set.seed(0)
n_quantiles = 10000
n_controls = 1e6
nname = "1M"
method = "resample"
n0all = (complete_data$casecontrol==0)%>%sum
n1all = (complete_data$casecontrol==1)%>%sum

param_xg = list(
  max_depth = 5,
  subsample = 0.2,
  eta = .05,
  objective = 'binary:logistic',
  eval_metric = 'auc',
  nthread=-1
)

# train-test split
complete_data = subsample_controls(complete_data, n_controls)
set.seed(0)
out = train_test_split(df=complete_data, weights=c(3, 1))
train = out[[1]]
test = out[[2]]
rm(complete_data)
gc()

# train-valid split
set.seed(0)
out = train_test_split(df=train, weights=c(2, 1))
rm(train);gc()
train_n = out[[1]]
valid_n = out[[2]]
rm(out);gc()
n0 = (train_n$casecontrol==0)%>%sum
n1 = (train_n$casecontrol==1)%>%sum
timestamp()
# produce training set
train_ = balanced_resample(train_n)
# imputation
quantiles = compute_quantiles(train_n, n_quantiles)
rm(train_n);gc()
set.seed(0)
train_ = impute_srs(train_, quantiles)
test_ = impute_srs(test, quantiles)
valid_ = impute_srs(valid_n, quantiles)

cases = list(
  "NONE", "DROP", 
  "NUM", "IND"
)

for(case in cases){
  if(case == "DROP"){
    train_s = train_ %>% filter(birth35) %>% select(-c(birthyear, birth35))
    test_s = test_ %>% filter(birth35) %>% select(-c(birthyear, birth35))
    valid_s = valid_ %>% filter(birth35) %>% select(-c(birthyear, birth35))
  }
  if(case == "NUM"){
    train_s = train_ %>% select(-c(birth35))
    test_s = test_ %>% select(-c(birth35))
    valid_s = valid_ %>% select(-c(birth35))
  }
  if(case == "IND"){
    train_s = train_ %>% select(-c(birthyear))
    test_s = test_ %>% select(-c(birthyear))
    valid_s = valid_ %>% select(-c(birthyear))
  }
  if(case == "NONE"){
    train_s = train_  %>% select(-c(birthyear, birth35))
    test_s = test_  %>% select(-c(birthyear, birth35))
    valid_s = valid_  %>% select(-c(birthyear, birth35))
  }
  
  dwatchlist = xgb_prep(train=train_s,
                        test=test_s,
                        valid=valid_s)
  dwatchlist$test = NULL
  
  # fit
  set.seed(0)
  n0all = n_controls
  param_xg$scale_pos_weight = (n1all/n0all)*(n0/n0)
  timestamp()
  xgb_fit = xgb.train(param_xg,
                      dwatchlist$train,
                      nrounds=2000,
                      dwatchlist,
                      verbose=1,print_every_n=10,
                      early_stopping_rounds=100)
  timestamp()
  cat("nfeatures: ")
  cat(xgb_fit$nfeatures, fill=T)
  # save
  out = list(
    xgb_fit=xgb_fit,
    quantiles=quantiles,
    test_ids=test_$id
  )
  filepath = paste0("R_data/results/models/year/XGB_", nname, "_ANY_year", case, ".rds")
  saveRDS(out, filepath)
  rm(xgb_fit, out, dwatchlist, train_s, test_s, valid_s)
}

