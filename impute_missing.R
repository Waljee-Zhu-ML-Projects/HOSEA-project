# imputation script (subsample)

# helper functions
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/impute_models.R')

# import data
complete_data <- readRDS('R_data/subsample/sub_complete_data_raw.rds')
complete_data_y <- select(complete_data,all_of(c('ID','CaseControl')))
complete_data <- select(complete_data,-all_of(c('ID','CaseControl')))

print(mean(is.na(complete_data)))

# variable blocks
# demographic (age,BMI,weight numeric, SmokeStatus multiclass, remaining 2-class)
demo_vars <- c('ageatindex','Gender','BMI','weight',
                'Asian','Black','HawaiianPacific','IndianAlaskan',
                'agentorange',
                'GerdAtIndex',
                'CHF','CTD','DEM','DIAB_C','HIV','MLD','MSLD','PARA','RD',
                'cd','copd','diab_nc','mi','pud','pvd')
demo_vars_num <- demo_vars[c(1,3,4)]
demo_vars_ind <- demo_vars[-c(1,3,4)]
smoke_vars <- c('smoke_current','smoke_former')
# colonoscopy vars
other_vars <- c('colonscopy_n','colonoscopy_max_diff',
                'labs_fobt_n','labs_fobt_max_diff',
                'H2R_int','H2R_mean','H2R_max','H2R_maxdiff','H2R_tv',
                'PPI_int','PPI_mean','PPI_max','PPI_maxdiff','PPI_tv')

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

#### impute zeros ####

# impute colonoscopies, fobt labs and meds with zeros
for(var in colnames(complete_data)){
  # allmeds variables
  if(substr(var,nchar(var)-3,nchar(var))=='_int'){
    complete_data[[var]] <- fill_by_zero(complete_data[[var]])
    print(paste0('Impute ',var))
  }
  if(substr(var,nchar(var)-4,nchar(var))=='_mean'){
    if(substr(var,0,3) %in% c('H2R','PPI')){
      complete_data[[var]] <- fill_by_zero(complete_data[[var]])
      print(paste0('Impute ',var))
    }
  }
  if(substr(var,nchar(var)-3,nchar(var))=='_max'){
    if(substr(var,0,3) %in% c('H2R','PPI')){
      complete_data[[var]] <- fill_by_zero(complete_data[[var]])
      print(paste0('Impute ',var))
    }
  }
  if(substr(var,nchar(var)-7,nchar(var))=='_maxdiff'){
    if(substr(var,0,3) %in% c('H2R','PPI')){
      complete_data[[var]] <- fill_by_zero(complete_data[[var]])
      print(paste0('Impute ',var))
    }
  }
  if(substr(var,nchar(var)-2,nchar(var))=='_tv'){
    if(substr(var,0,3) %in% c('H2R','PPI')){
      complete_data[[var]] <- fill_by_zero(complete_data[[var]])
      print(paste0('Impute ',var))
    }
  }
  if(substr(var,nchar(var)-3,nchar(var))=='_max'){
    if(substr(var,0,3) %in% c('H2R','PPI')){
      complete_data[[var]] <- fill_by_zero(complete_data[[var]])
      print(paste0('Impute ',var))
    }
  }
  # events data
  if(substr(var,nchar(var)-4,nchar(var))=='_diff'){
    complete_data[[var]] <- fill_by_zero(complete_data[[var]])
    print(paste0('Impute ',var))
  }
  if(substr(var,nchar(var)-1,nchar(var))=='_n'){
    complete_data[[var]] <- fill_by_zero(complete_data[[var]])
    print(paste0('Impute ',var))
  }
}

# save data with no imputation
saveRDS(complete_data,
        file='R_data/subsample/sub_complete_data.rds')

print(mean(is.na(complete_data)))

#### initial imputation ####

complete_data_init <- complete_data
# original complete data has NAs to determine original missing vals

# blood labs
for(lab in lab_vars){
  for(v in lab){
    for(summ in lab_suffix){
      varname <- paste0(v,summ)
      complete_data_init[[varname]] <- fill_by_sample(complete_data_init[[varname]])
      print(paste0('Impute ',varname))
    }
  }
}

# demographic variables
for(varname in demo_vars){
    complete_data_init[[varname]] <- fill_by_sample(complete_data_init[[varname]])
    print(paste0('Impute ',varname))
}

# smoking indicator
smoke_miss <- is.na(complete_data_init[['smoke_current']])
complete_data_init[smoke_miss,smoke_vars] <- complete_data_init[sample(which(!smoke_miss),sum(smoke_miss),replace=T),smoke_vars]

# save initial (random sample) imputation
saveRDS(complete_data_init,
        file='R_data/subsample/sub_complete_data_impsamp.rds')

#### model-based imputation ####

# order:
# lab means (BMP, Lipid, CBC, LFT, A1c, CRP)
# lab max/min
# lab maxdiff/mindiff
# lab TV
# demo vars (smoking, charlson, race, then age/BMI/weight)

complete_data_impute <- complete_data_init

# lab order
lab_order <- c(2,6,3,5,1,4)
# number of imputation cycles
ncycles <- 4

for(cc in 1:ncycles){
  # impute lab means
  # temporary data for imputing
  temp_data <- select(complete_data_impute,
                      all_of(c(demo_vars,other_vars)))
  # impute
  print(paste0('Cycle ',cc,', lab means'))
  for(ll in lab_order){
    for(v in lab_vars[[ll]]){
      var <- paste0(v,'_mean')
      # missing indices from complete_data
      imiss <- is.na(complete_data[[var]])
      # impute
      temp <- impute_reg_numeric(temp_data,
                                 var,
                                 imiss)
      # update complete_data_impute
      complete_data_impute[[var]] <- temp
      # print progress 
      print(paste0('Impute ',var))
    }
  }
  
  # impute lab max/min
  # temporary data for imputing
  temp_data <- select(complete_data_impute,
                      all_of(c(demo_vars,other_vars,
                               paste0(unlist(lab_vars),'_mean'))))
  # impute
  print(paste0('Cycle ',cc,', lab max/min'))
  for(ll in lab_order){
    for(v in lab_vars[[ll]]){
      for(summ in c('_max','_min')){
        var <- paste0(v,summ)
        # missing indices from complete_data
        imiss <- is.na(complete_data[[var]])
        # impute
        temp <- impute_reg_numeric(temp_data,
                                   var,
                                   imiss)
        # update complete_data_impute
        complete_data_impute[[var]] <- temp
        # print progress 
        print(paste0('Impute ',var))
      }
    }
  }
  
  # impute lab maxdiff/mindiff
  # temporary data for imputing
  temp_data <- select(complete_data_impute,
                      all_of(c(demo_vars,other_vars,
                               paste0(unlist(lab_vars),'_mean'),
                               paste0(unlist(lab_vars),'_max'),
                               paste0(unlist(lab_vars),'_min'))))
  # impute
  print(paste0('Cycle ',cc,', lab max/min slopes'))
  for(ll in lab_order){
    for(v in lab_vars[[ll]]){
      for(summ in c('_maxdiff','_mindiff')){
        var <- paste0(v,summ)
        # missing indices from complete_data
        imiss <- is.na(complete_data[[var]])
        # impute
        temp <- impute_reg_numeric(temp_data,
                                   var,
                                   imiss)
        # update complete_data_impute
        complete_data_impute[[var]] <- temp
        # print progress 
        print(paste0('Impute ',var))
      }
    }
  }
  
  # impute lab total variation
  # temporary data for imputing
  temp_data <- select(complete_data_impute,
                      all_of(c(demo_vars,other_vars,
                               paste0(unlist(lab_vars),'_mean'),
                               paste0(unlist(lab_vars),'_max'),
                               paste0(unlist(lab_vars),'_min'),
                               paste0(unlist(lab_vars),'_maxdiff'),
                               paste0(unlist(lab_vars),'_mindiff'))))
  # impute
  print(paste0('Cycle ',cc,', lab total variation'))
  for(ll in lab_order){
    for(v in lab_vars[[ll]]){
      var <- paste0(v,'_tv')
      # missing indices from complete_data
      imiss <- is.na(complete_data[[var]])
      # impute
      temp <- impute_reg_numeric(temp_data,
                                 var,
                                 imiss)
      # update complete_data_impute
      complete_data_impute[[var]] <- temp
      # print progress 
      print(paste0('Impute ',var))
    }
  }
  
  # impute smoking status using all other vars
  print(paste0('Cycle ',cc,', demographic variables'))
  imiss <- is.na(complete_data[['smoke_current']])
  temp <- impute_reg_multi(complete_data_impute,
                           smoke_vars,
                           imiss)
  complete_data_impute[,smoke_vars] <- temp
  print('Impute SmokeStatus')
  
  # impute demographic variables (binary)
  for(var in rev(demo_vars_ind)){
    imiss <- is.na(complete_data[[var]])
    temp <- impute_reg_binary(complete_data_impute,
                              var,
                              imiss)
    complete_data_impute[[var]] <- temp
    print(paste0('Impute ',var))
  }
  
  # impute demographic variables (numeric)
  for(var in demo_vars_num){
    imiss <- is.na(complete_data[[var]])
    temp <- impute_reg_numeric(complete_data_impute,
                               var,
                               imiss)
    complete_data_impute[[var]] <- temp
    print(paste0('Impute ',var))
  }
  
}

# save completed imputed data
saveRDS(complete_data_impute,
        file='R_data/subsample/sub_complete_data_impreg.rds')

