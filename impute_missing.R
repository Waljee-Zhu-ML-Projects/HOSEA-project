# imputation script (subsample)

#### global variables ####

# demographic (age,bmi,weight numeric, SmokeStatus multiclass, remaining 2-class)
demo_vars <- c(
  'ageatindex','Gender','bmi','weight',
  'Asian','Black','HawaiianPacific','IndianAlaskan',
  'agentorange','GerdAtIndex'
)
demo_vars_num <- demo_vars[c(1,3,4)]
demo_vars_ind <- demo_vars[-c(1,3,4)]
charlson_vars = c(
  'CHF','CTD','DEM','DIAB_C','HIV','MLD','MSLD','PARA','RD',
  'cd','copd','diab_nc','mi','pud','pvd'
)
smoke_vars <- c('smoke_current','smoke_former')
# colonoscopy vars
other_vars <- c('colonoscopy_n','colonoscopy_maxdiff',
                'labs_fobt_n','labs_fobt_maxdiff',
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


# specify seed for sampling
# specify ncycles
# [DEPRECATED, kept for consistency] hybrid_reg 
# if test is missing, it is assume we are imputing the training set directly
impute_missing_hosea <- function(train,test=train,ncycles=5,seed=1,hybrid_reg=FALSE){
  
  #### pre-processing ####
  # import data
  train_y <- select(train,all_of(c('ID','CaseControl')))
  train_x <- select(train,-all_of(c('ID','CaseControl')))
  test_y <- select(test,all_of(c('ID','CaseControl')))
  test_x <- select(test,-all_of(c('ID','CaseControl')))
  print("Original data")
  logging(train_x, test_x)
  
  #### impute zeros ####
  
  print('Zero imputation for events')
  # impute colonoscopies, fobt labs and meds with zeros
  for(var in other_vars){
    train_x[[var]] <- fill_by_zero(train_x[[var]])
    test_x[[var]] <- fill_by_zero(test_x[[var]])
  }
  logging(train_x, test_x)
  
  #### sampling imputation ####
  
  set.seed(seed)
  # print('Sampling imputation')
  # df = sampling_imputation(train_x, test_x)
  # test_sample = lab_consistency(lab_vars, lab, df)
  # logging(test_sample)
  
  print('Sampling imputation (MICE)')
  df = sampling_imputation_mice(train_x, test_x)
  test_sample = lab_consistency(lab_vars, lab, df)
  logging(test_sample)
  
  #### median-based imputation ####
  
  print('Median imputation')
  df = median_imputation(train_x, test_x)
  test_med = lab_consistency(lab_vars, lab, df)
  logging(test_med)
  
  #### model-based imputation ####
  
  set.seed(seed)
  print('Model-based imputation using MICE')
  df = model_imputation(train_x, test_x)
  test_model = lab_consistency(lab_vars, lab, df)
  logging(test_model)
  
  #### return complete imputed data ####
  return(list(clean=bind_cols(test_y,test_x),
              impmed=bind_cols(test_y,test_med),
              impsamp=bind_cols(test_y,test_samp),
              impreg=bind_cols(test_y,test_model)))
}


logging = function(train_x, test_x=NULL) {
  print(paste("Missing proportion in training set:", mean(is.na(train_x))))
  if(!missing(test_x)) print(paste("Missing proportion in training set:", mean(is.na(test_x))))
}

median_imputation = function(train_x, test_x) {
  test_med = test_x
  # blood labs
  for(lab in lab_vars){
    for(v in lab){
      for(summ in lab_suffix){
        varname <- paste0(v,summ)
        test_med[[varname]] = fill_by_median(train_x[[varname]], test_med[[varname]])
      }
    }
  }
  # other variables
  for(varname in c(demo_vars, other_vars, smoke_vars, charlson_vars)){
    test_med[[varname]] = fill_by_median(train_x[[varname]], test_med[[varname]])
  }
  # TODO: no missing values in other columns?
  return(test_med)
}

model_imputation = function(train_x, test_x) {
  # merged to pass to mice
  merged_x = bind_rows(train_x, test_x)
  train_n = nrow(train_x)
  test_n = nrow(test_x)
  ignore = rep(T, train_n + test_n)
  ignore[1:train_n] = F
  
  # call mice
  # ignore tells mice to only fit the model on the training set
  # but still outputs imputation for the testing set
  mice_result = mice::mice(merged_x, method=NULL, m=1, maxit=10, 
                           defaultMethod = c("norm", "logreg", "polyreg", "polr"),
                           visitSequence="monotone", ignore=ignore)
  merged_imputed_x = mice::complete(mice_result)
  
  # get testing set
  df = merged_imputed_x[-(1:train_n), ]
  
  return(tibble(df))
}

sampling_imputation_mice = function(train_x, test_x) {
  # merged to pass to mice
  merged_x = bind_rows(train_x, test_x)
  train_n = nrow(train_x)
  test_n = nrow(test_x)
  ignore = rep(T, train_n + test_n)
  ignore[1:train_n] = F
  
  # call mice
  # ignore tells mice to only fit the model on the training set
  # but still outputs imputation for the testing set
  mice_result = mice::mice(merged_x, m=1, maxit=10, 
                           method="sample",
                           visitSequence="monotone", ignore=ignore)
  merged_imputed_x = mice::complete(mice_result)
  
  # get testing set
  df = merged_imputed_x[-(1:train_n), ]
  
  return(tibble(df))
}

sampling_imputation = function(train_x, test_x){
  test_samp = test_x
  
  # blood labs
  for(lab in lab_vars){
    for(v in lab){
      for(summ in lab_suffix){
        varname <- paste0(v,summ)
        test_samp[[varname]] =
          fill_by_sample(train_x[[varname]], test_samp[[varname]])
      }
    }
  }
  
  # demographic and smoke variables
  for(varname in c(demo_vars, other_vars, smoke_vars, charlson_vars)){
    test_samp[[varname]] = 
      fill_by_sample(train_x[[varname]], test_samp[[varname]])
  }
}

lab_consistency = function(lab_vars, lab, df) {
  for(lab in lab_vars){
    for(v in lab){
      v_mean = paste0(v,"_mean")
      v_max = paste0(v,"_max")
      v_min = paste0(v,"_min")
      v_maxdiff = paste0(v,"_maxdiff")
      v_mindiff = paste0(v,"_mindiff")
      v_tv = paste0(v,"_tv")
      # flip min and max if incorrect order
      flip = df[[v_min]] > df[[v_max]]
      df[[v_min]] = ifelse(flip, df[[v_max]], df[[v_min]])
      df[[v_max]] = ifelse(flip, df[[v_min]], df[[v_max]])
      # if mean is not in between, replace by average
      between = (df[[v_mean]] >= df[[v_min]]) %AND% (df[[v_max]] >= df[[v_mean]])
      df[[v_mean]] = ifelse(between, df[[v_mean]], (df[[v_max]]-df[[v_min]])/2)
      # flip diff if incorrect order
      flip = df[[v_mindiff]] > df[[v_maxdiff]]
      df[[v_mindiff]] = ifelse(flip, df[[v_maxdiff]], df[[v_mindiff]])
      df[[v_maxdiff]] = ifelse(flip, df[[v_mindiff]], df[[v_maxdiff]])
      # tv should be between abs(min) and abs(max), otherwise replace by average
      tv = (abs(df[[v_mindiff]]) + abs(df[[v_maxdiff]])) / 2.
      between = (df[[v_tv]] >= abs(df[[v_mindiff]])) %AND% (abs(df[[v_maxdiff]]) >= df[[v_tv]])
      df[[v_tv]] = ifelse(between, df[[v_tv]], tv)
    }
  }
  return(df)
}
