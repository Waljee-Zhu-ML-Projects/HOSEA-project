# imputation script (subsample)

# specify seed for sampling
# specify ncycles
impute_missing_hosea <- function(data_raw,ncycles=4,seed=1){
  
  # import data
  complete_data_y <- select(data_raw,all_of(c('ID','CaseControl')))
  complete_data <- select(data_raw,-all_of(c('ID','CaseControl')))
  
  print(mean(is.na(complete_data)))
  
  # variable blocks
  # demographic (age,BMI,weight numeric, SmokeStatus multiclass, remaining 2-class)
  # will need BMI -> bmi
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
  other_vars <- c('colonoscopy_n','colonoscopy_max_diff',
                  'labs_fobt_n','labs_fobt_max_diff',
                  'H2R_int','H2R_mean','H2R_max','H2R_max_diff','H2R_tv',
                  'PPI_int','PPI_mean','PPI_max','PPI_max_diff','PPI_tv')
  
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
                  '_max_diff','_min_diff','_tv')
  
  #### impute zeros ####
  
  # impute colonoscopies, fobt labs and meds with zeros
  for(var in other_vars){
    complete_data[[var]] <- fill_by_zero(complete_data[[var]])
    #print(paste0('Impute ',var))
  }
  
  print('zero imputation for events')
  print(mean(is.na(complete_data)))
  
  #### initial imputation ####
  
  complete_data_init <- complete_data
  # original complete data has NAs to determine original missing vals
  
  # blood labs
  set.seed(seed)
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
  
  print('Initial imputation')
  print(mean(is.na(complete_data_init)))
  
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
  # ncycles <- 4
  
  for(cc in 1:ncycles){
    # impute lab means
    # temporary data for imputing
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars,smoke_vars,other_vars)))
    # impute
    print(paste0('Cycle ',cc,', lab means'))
    for(ll in lab_order){
      for(v in lab_vars[[ll]]){
        var <- paste0(v,'_mean')
        # missing indices from complete_data
        imiss <- is.na(complete_data[[var]])
        # impute
        temp <- impute_simple_numeric(temp_data,
                                      complete_data_impute[[var]],
                                      imiss,
                                      ridge=1)
        # update complete_data_impute
        complete_data_impute[[var]] <- temp
        # print progress 
        #print(paste0('Impute ',var))
      }
    }
    
    # impute lab max/min
    # temporary data for imputing
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars,smoke_vars,other_vars,
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
          temp <- impute_simple_numeric(temp_data,
                                        complete_data_impute[[var]],
                                        imiss,
                                        nneg=TRUE,
                                        ridge=1)
          # update complete_data_impute
          complete_data_impute[[var]] <- temp
          # print progress 
          #print(paste0('Impute ',var))
        }
      }
    }
    
    # impute lab maxdiff/mindiff
    # temporary data for imputing
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars,smoke_vars,other_vars,
                                 paste0(unlist(lab_vars),'_mean'),
                                 paste0(unlist(lab_vars),'_max'),
                                 paste0(unlist(lab_vars),'_min'))))
    # impute
    print(paste0('Cycle ',cc,', lab max/min slopes'))
    for(ll in lab_order){
      for(v in lab_vars[[ll]]){
        for(summ in c('_max_diff','_min_diff')){
          var <- paste0(v,summ)
          # missing indices from complete_data
          imiss <- is.na(complete_data[[var]])
          # impute
          temp <- impute_simple_numeric(temp_data,
                                        complete_data_impute[[var]],
                                        imiss,
                                        ridge=1)
          # update complete_data_impute
          complete_data_impute[[var]] <- temp
          # print progress 
          #print(paste0('Impute ',var))
        }
      }
    }
    
    # impute lab total variation
    # temporary data for imputing
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars,smoke_vars,other_vars,
                                 paste0(unlist(lab_vars),'_mean'),
                                 paste0(unlist(lab_vars),'_max'),
                                 paste0(unlist(lab_vars),'_min'),
                                 paste0(unlist(lab_vars),'_max_diff'),
                                 paste0(unlist(lab_vars),'_min_diff'))))
    # impute
    print(paste0('Cycle ',cc,', lab total variation'))
    for(ll in lab_order){
      for(v in lab_vars[[ll]]){
        var <- paste0(v,'_tv')
        # missing indices from complete_data
        imiss <- is.na(complete_data[[var]])
        # impute
        temp <- impute_simple_numeric(temp_data,
                                      complete_data_impute[[var]],
                                      imiss,
                                      nneg=TRUE,
                                      ridge=1)
        # update complete_data_impute
        complete_data_impute[[var]] <- temp
        # print progress 
        #print(paste0('Impute ',var))
      }
    }
    # impute demographic vars using lab means
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars,other_vars,
                                 paste0(unlist(lab_vars),'_mean'))))
    # impute smoking status using all other vars
    print(paste0('Cycle ',cc,', demographic variables'))
    imiss <- is.na(complete_data[['smoke_current']])
    temp <- impute_simple_multi(temp_data,
                                complete_data_impute[,smoke_vars],
                                imiss,
                                ridge=1)
    complete_data_impute[,smoke_vars] <- temp
    #print('Impute SmokeStatus')
    
    # impute demographic variables (binary) using lab means
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars_num,smoke_vars,other_vars,
                                 paste0(unlist(lab_vars),'_mean'))))
    # impute
    for(var in rev(demo_vars_ind)){
      imiss <- is.na(complete_data[[var]])
      temp <- impute_simple_numeric(temp_data,
                                    complete_data_impute[[var]],
                                    imiss,
                                    binary=TRUE,
                                    ridge=1)
      complete_data_impute[[var]] <- temp
      #print(paste0('Impute ',var))
    }
    
    # impute demographic variables (numeric)
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars_ind,smoke_vars,other_vars,
                                 paste0(unlist(lab_vars),'_mean'))))
    # impute
    for(var in demo_vars_num){
      imiss <- is.na(complete_data[[var]])
      temp <- impute_simple_numeric(temp_data,
                                    complete_data_impute[[var]],
                                    imiss,
                                    ridge=1)
      complete_data_impute[[var]] <- temp
      #print(paste0('Impute ',var))
    }
  }
  # return complete imputed data
  return(list(clean=bind_cols(complete_data_y,complete_data),
              impsamp=bind_cols(complete_data_y,complete_data_init),
              impreg=bind_cols(complete_data_y,complete_data_impute)))
}
