# overnight script, 07/05

# helper functions (mmean, mmin, mmax, simple imputation)
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_xgb.R')

print('Start')
timestamp()

# # begin with demographic table
# sample_tab <- readRDS('R_data/sample.rds')
# demo_vars <- c('ID','CaseControl',
#                'ageatindex','Gender','bmi','weight',
#                'Asian','Black','HawaiianPacific','IndianAlaskan',
#                'SmokeStatus',
#                'agentorange',
#                'GerdAtIndex',
#                'CHF','CTD','DEM','DIAB_C','HIV','MLD','MSLD','PARA','RD',
#                'cd','copd','diab_nc','mi','pud','pvd')
# 
# complete_data <- sample_tab[,demo_vars] # 28 columns
# 
# # convert factors to integer and expand multiple smoking classes
# # Gender
# complete_data$Gender[complete_data$Gender==''] <- NA
# complete_data$Gender <- as.integer(complete_data$Gender=='M')
# # agentorange
# complete_data$agentorange <- as.integer(complete_data$agentorange=='YES')
# # SmokeStatus (keep current and former)
# complete_data$smoke_current <- as.integer(complete_data$SmokeStatus==1)
# complete_data$smoke_former <- as.integer(complete_data$SmokeStatus==2)
# complete_data <- select(complete_data,-SmokeStatus) # adds one extra variable (29 columns)
# 
# # event tables
# event_tables <- c('colonoscopy',
#                   'labs_fobt')
# # removed EGD and other procedures
# 
# # join summary tables (33 columns)
# for(tab in event_tables){
#   print(paste0('For table ',tab,':'))
#   table <- readRDS(paste0('R_data/',tab,'_summary.rds'))
#   complete_data <- left_join(complete_data,table,by="ID")
#   print('read in and join')
#   rm(table)
# }
# 
# # medication table (43 columns)
# print('For table: allmeds')
# table <- readRDS('R_data/allmeds_summary.rds')
# complete_data <- left_join(complete_data,table,by="ID")
# print('read in and join')
# rm(table)
# 
# # valued tables
# value_tables <- c('labs_a1c',
#                   'labs_bmp',
#                   'labs_cbc',
#                   'labs_crp',
#                   'labs_lft',
#                   'labs_lipid')
# 
# # join summary tables (241 columns)
# for(tab in value_tables){
#   print(paste0('For table ',tab,':'))
#   table <- readRDS(paste0('R_data/',tab,'_summary.rds'))
#   complete_data <- left_join(complete_data,table,by="ID")
#   print('read in and join')
#   rm(table)
# }
# 
# # save complete data (241 columns)
# saveRDS(complete_data,
#         file='R_data/complete_data_raw.rds')

print('data linked')
timestamp()

# variable blocks
# demographic (age,bmi,weight numeric, SmokeStatus multiclass, remaining 2-class)
demo_vars <- c('ageatindex','Gender','bmi','weight',
               'Asian','Black','HawaiianPacific','IndianAlaskan',
               'agentorange',
               'GerdAtIndex',
               'CHF','CTD','DEM','DIAB_C','HIV','MLD','MSLD','PARA','RD',
               'cd','copd','diab_nc','mi','pud','pvd')
demo_vars_num <- demo_vars[c(1,3,4)]
demo_vars_ind <- demo_vars[-c(1,3,4)]
smoke_vars <- c('smoke_current','smoke_former')
# colonoscopy vars
other_vars <- c('colonoscopy_n','colonoscopy_max_diff',
                'labs_fobt_n','labs_fobt_max_diff',
                'h2r_int','h2r_mean','h2r_max','h2r_maxdiff','h2r_tv',
                'ppi_int','ppi_mean','ppi_max','ppi_maxdiff','ppi_tv')

# blood lab measurements (all numeric)
lab_vars <- list(c('A1c'),
                 c('bun','calc','chlor','co2','creat','gluc','k','na'),
                 c('baso','eos','hct','hgb','lymph','mch','mchc','mcv','mono','mpv','neut','platelet','rbc','rdw','wbc'),
                 c('CRP'),
                 c('alkphos','alt','ast','totprot'),
                 c('chol','hdl','ldl','trig'))
# blood lab longitudinal summaries
lab_suffix <- c('_mean',
                '_max','_min',
                '_maxdiff','_mindiff','_tv')

# impute some vars with zero, save
complete_data_clean <- complete_data

# impute colonoscopies, fobt labs and meds with zeros
for(var in other_vars){
  complete_data_clean[[var]] <- fill_by_zero(complete_data_clean[[var]])
  #print(paste0('Impute ',var))
}

saveRDS(complete_data_clean,
        file='R_data/complete_data_clean.rds')

# impute remaining vars with sample, save
complete_data_init <- complete_data_clean
# original complete data has NAs to determine original missing vals

# blood labs
set.seed(1999)
for(lab in lab_vars){
  for(v in lab){
    for(summ in lab_suffix){
      varname <- paste0(v,summ)
      complete_data_init[[varname]] <- fill_by_sample(complete_data_init[[varname]])
      #print(paste0('Impute ',varname))
    }
  }
}

# demographic variables
for(varname in demo_vars){
  complete_data_init[[varname]] <- fill_by_sample(complete_data_init[[varname]])
  #print(paste0('Impute ',varname))
}

# smoking indicator
smoke_miss <- is.na(complete_data_init[['smoke_current']])
complete_data_init[smoke_miss,smoke_vars] <- complete_data_init[sample(which(!smoke_miss),sum(smoke_miss),replace=T),smoke_vars]
#print('Impute smoking status')

#saveRDS(complete_data_init,
#        file='R_data/complete_data_samp.rds')

print('Data prepped')
timestamp()

# test/train/validation split
# take a test set, 10% of observations
set.seed(1994)
n <- nrow(complete_data_init)
ntest <- as.integer(.1*n)
itest <- sample(1:n,ntest) 
test <- rep(FALSE,n)
test[itest] <- TRUE
train <- !test

# select a validation set (10% of observations)
set.seed(1997)
ivalid <- sample(which(train),size=as.integer(.1*n),replace=FALSE)
valid <- rep(FALSE,n)
valid[ivalid] <- TRUE
traintrain <- as.logical(train*(!valid))

# xgb formatting for each set
# weights would go in here as well
dtrain <- xgb.DMatrix(as.matrix(complete_data_init[traintrain,-c(1,2)]),
                      label=complete_data_init$CaseControl[traintrain])
dvalid <-  xgb.DMatrix(as.matrix(complete_data_init[valid,-c(1,2)]),
                       label=complete_data_init$CaseControl[valid])
dtest <-  xgb.DMatrix(as.matrix(complete_data_init[test,-c(1,2)]),
                      label=complete_data_init$CaseControl[test])
# combine as a watchlist
dwatchlist <- list(train=dtrain,test=dtest,valid=dvalid)

# xgboost parameters
param_xg <- list(max_depth=3,eta=.1,objective='binary:logistic',
                 eval_metric='logloss')


print('Begin fitting')
timestamp()

# fit initial xgboost model
xgb_fit_old <- xgb.train(param_xg,
                          dwatchlist$train, # training set
                          nrounds=10,
                          dwatchlist, # data watchlist
                          verbose=1,print_every_n=5,
                          early_stopping_rounds=5)

print('Initial fit')
timestamp()

for(ii in 1:24){
  xgb_fit_new <- xgb.train(param_xg,
                           dwatchlist$train, # training set
                           nrounds=10,
                           dwatchlist, # data watchlist
                           verbose=1,print_every_n=5,
                           early_stopping_rounds=5,
                           xgb_model = xgb_fit_old)
  
  print(paste0('Fit up to ',10*(ii+1),' trees'))
  timestamp()
  
  xgb_fit_old <- xgb_fit_new
}

saveRDS(xgb_fit_new,
        file='R_data/models/xgb_samp.rds')
