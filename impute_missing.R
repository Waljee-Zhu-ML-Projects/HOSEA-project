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
impute_missing_hosea <- function(train,test=NULL,valid=NULL,
                                 ncycles=5,seed=1,hybrid_reg=FALSE){
  
  #### pre-processing ####
  # import data
  train_y <- select(train,all_of(c('ID','CaseControl')))
  train_x <- select(train,-all_of(c('ID','CaseControl')))
  valid_y <- select(valid,all_of(c('ID','CaseControl')))
  valid_x <- select(valid,-all_of(c('ID','CaseControl')))
  test_y <- select(test,all_of(c('ID','CaseControl')))
  test_x <- select(test,-all_of(c('ID','CaseControl')))
  cat("=== Original data ===", fill=T)
  logging(train_x, valid_x, test_x)
  
  #### impute zeros ####
  
  cat('\n=== Zero imputation for events ===', fill=T)
  # impute colonoscopies, fobt labs and meds with zeros
  for(var in other_vars){
    train_x[[var]] <- fill_by_zero(train_x[[var]])
    test_x[[var]] <- fill_by_zero(test_x[[var]])
    valid_x[[var]] <- fill_by_zero(valid_x[[var]])
  }
  logging(train_x, valid_x, test_x)
  
  #### sampling imputation ####
  
  set.seed(seed)
  # print('Sampling imputation')
  # df = sampling_imputation(train_x, test_x)
  # test_sample = lab_consistency(lab_vars, lab, df)
  # logging(test_sample)
  
  cat('\n=== Sampling imputation (MICE) ===', fill=T)
  df = mice_imputation(train_x, valid_x, test_x, method="sample")
  train_sample = lab_consistency(lab_vars, lab, df$train)
  valid_sample = lab_consistency(lab_vars, lab, df$valid)
  test_sample = lab_consistency(lab_vars, lab, df$test)
  logging(train_sample, valid_sample, test_sample)
  
  #### median-based imputation ####
  
  cat('\n=== Median imputation ===', fill=T)
  df = median_imputation(train_x, valid_x, test_x)
  train_med = lab_consistency(lab_vars, lab, df$train)
  valid_med = lab_consistency(lab_vars, lab, df$valid)
  test_med = lab_consistency(lab_vars, lab, df$test)
  logging(train_med, valid_med, test_med)
  
  #### model-based imputation ####
  
  set.seed(seed)
  cat('\n=== Model-based imputation (MICE w/ CART) ===', fill=T)
  df = mice_imputation(train_x, valid_x, test_x, method="cart")
  train_model = lab_consistency(lab_vars, lab, df$train)
  valid_model = lab_consistency(lab_vars, lab, df$valid)
  test_model = lab_consistency(lab_vars, lab, df$test)
  logging(train_model, valid_model, test_model)
  
  #### prepare output ####
  train = list(clean=bind_cols(train_y,train_x),
               impmed=bind_cols(train_y,train_med),
               impsamp=bind_cols(train_y,train_sample),
               impreg=bind_cols(train_y,train_model))
  valid = list(clean=bind_cols(valid_y,valid_x),
               impmed=bind_cols(valid_y,valid_med),
               impsamp=bind_cols(valid_y,valid_sample),
               impreg=bind_cols(valid_y,valid_model))
  test = list(clean=bind_cols(test_y,test_x),
              impmed=bind_cols(test_y,test_med),
              impsamp=bind_cols(test_y,test_sample),
              impreg=bind_cols(test_y,test_model))
  
  #### return complete imputed data ####
  return(list(train=train, valid=valid, test=test))
}


logging = function(train_x, valid_x=NULL, test_x=NULL) {
  cat(paste("Missing proportion in training set:", mean(is.na(train_x))), fill=T)
  if(!missing(valid_x)) cat(paste("Missing proportion in validation set:", mean(is.na(valid_x))), fill=T)
  if(!missing(test_x)) cat(paste("Missing proportion in training set:", mean(is.na(test_x))), fill=T)
}

median_imputation = function(train_x, valid_x, test_x) {
  train_med = train_x
  valid_med = valid_x
  test_med = test_x
  # blood labs
  for(lab in lab_vars){
    for(v in lab){
      for(summ in lab_suffix){
        varname <- paste0(v,summ)
        train_med[[varname]] = fill_by_median(train_x[[varname]], train_med[[varname]])
        valid_med[[varname]] = fill_by_median(train_x[[varname]], valid_med[[varname]])
        test_med[[varname]] = fill_by_median(train_x[[varname]], test_med[[varname]])
      }
    }
  }
  # other variables
  for(varname in c(demo_vars, other_vars, smoke_vars, charlson_vars)){
    train_med[[varname]] = fill_by_median(train_x[[varname]], train_med[[varname]])
    valid_med[[varname]] = fill_by_median(train_x[[varname]], valid_med[[varname]])
    test_med[[varname]] = fill_by_median(train_x[[varname]], test_med[[varname]])
  }
  # TODO: no missing values in other columns?
  return(list(train=train_med, valid=valid_med, test=test_med))
}

mice_imputation = function(train_x, valid_x, test_x, method, ...) {
  # merged to pass to mice
  merged_x = bind_rows(train_x, valid_x, test_x)
  train_n = nrow(train_x)
  valid_n = nrow(valid_x)
  test_n = nrow(test_x)
  ignore = rep(T, train_n + test_n + valid_n)
  ignore[1:train_n] = F
  
  # call mice
  # ignore tells mice to only fit the model on the training set
  # but still outputs imputation for the testing set
  mice_result = mice::mice(merged_x, m=1, method=method,
                           maxit=ifelse(method=="sample", 1, 5),
                           visitSequence="monotone", ignore=ignore, ...)
  merged_imputed_x = mice::complete(mice_result)
  
  out = list()
  # get training set
  df = merged_imputed_x[1:train_n, ]
  out$train = tibble(df)
  # get validation set
  df = merged_imputed_x[(train_n+1):(train_n+valid_n), ]
  out$valid = tibble(df)
  # get testing set
  df = merged_imputed_x[(train_n+valid_n+1):(train_n+valid_n+test_n), ]
  out$test = tibble(df)
  
  return(out)
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
      between = (df[[v_mean]] >= df[[v_min]]) && (df[[v_max]] >= df[[v_mean]])
      df[[v_mean]] = ifelse(between, df[[v_mean]], (df[[v_max]]-df[[v_min]])/2)
      # flip diff if incorrect order
      flip = df[[v_mindiff]] > df[[v_maxdiff]]
      df[[v_mindiff]] = ifelse(flip, df[[v_maxdiff]], df[[v_mindiff]])
      df[[v_maxdiff]] = ifelse(flip, df[[v_mindiff]], df[[v_maxdiff]])
      # tv should be between abs(min) and abs(max), otherwise replace by average
      tv = (abs(df[[v_mindiff]]) + abs(df[[v_maxdiff]])) / 2.
      between = (df[[v_tv]] >= abs(df[[v_mindiff]])) && (abs(df[[v_maxdiff]]) >= df[[v_tv]])
      df[[v_tv]] = ifelse(between, df[[v_tv]], tv)
    }
  }
  return(df)
}
