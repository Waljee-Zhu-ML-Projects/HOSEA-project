# imputation script (subsample)

# specify seed for sampling
# specify ncycles
impute_missing_hosea <- function(data_raw,ncycles=5,seed=1){
  
  # import data
  complete_data_y <- select(data_raw,all_of(c('ID','CaseControl')))
  complete_data_raw <- select(data_raw,-all_of(c('ID','CaseControl')))
  
  print(mean(is.na(complete_data_raw)))
  
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
    complete_data_raw[[var]] <- fill_by_zero(complete_data_raw[[var]])
    #print(paste0('Impute ',var))
  }
  
  print('zero imputation for events')
  print(mean(is.na(complete_data_raw)))
  
  #### initial imputation ####
  
  complete_data_init <- complete_data_raw
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
  
  #### median-based imputation ####
  
  complete_data_med <- complete_data_init
  
  # blood labs
  for(lab in lab_vars){
    for(v in lab){
      for(summ in lab_suffix){
        # variable name
        varname <- paste0(v,summ)
        # missing indices (from raw data)
        imiss <- is.na(complete_data_raw[[varname]])
        # fill with median
        complete_data_med[[varname]][imiss] <- NA
        complete_data_med[[varname]] <- fill_by_median(complete_data_med[[varname]])
        #print(paste0('Impute ',varname))
      }
    }
  }
  
  # demographic variables
  for(varname in demo_vars){
    # missing indices (from raw data)
    imiss <- is.na(complete_data_raw[[varname]])
    # fill with median
    complete_data_med[[varname]][imiss] <- NA
    complete_data_med[[varname]] <- fill_by_median(complete_data_med[[varname]])
    #print(paste0('Impute ',varname))
  }
  
  print('Median imputation')
  
  #### model-based imputation ####
  
  # order:
  # lab means (BMP, Lipid, CBC, LFT, A1c, CRP)
  # lab max/min
  # lab maxdiff/mindiff
  # lab TV
  # demo vars (smoking, charlson, race, then age/bmi/weight)
  
  #complete_data_impute <- complete_data_init
  complete_data_impute <- data.frame(lapply(complete_data_init,quantile_normalize))
  
  # store missing values to check convergence
  all_miss_old <- c(complete_data_impute[is.na(complete_data_raw)])
  # print(sqrt(mean(all_miss_old^2)))
  
  # lab order
  lab_order <- c(2,6,3,5,1)
  # leave CRP with intial sample imputation, not enough data to
  # impute cleanly
  
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
        imiss <- is.na(complete_data_raw[[var]])
        # impute
        temp <- impute_reg(temp_data,
                           complete_data_impute[[var]],
                           imiss)
        # update complete_data_impute
        complete_data_impute[[var]] <- temp
        # print progress 
        # print(paste0('Impute ',var))
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
          imiss <- is.na(complete_data_raw[[var]])
          # impute
          temp <- impute_reg(temp_data,
                             complete_data_impute[[var]],
                             imiss,
                             nneg=F)
          # update complete_data_impute
          complete_data_impute[[var]] <- temp
          # print progress 
          # print(paste0('Impute ',var))
        }
      }
    }
    
    # impute lab maxdiff/mindiff
    # temporary data for imputing
    # temp_data <- select(complete_data_impute,
    #                     all_of(c(demo_vars,smoke_vars,other_vars,
    #                              paste0(unlist(lab_vars),'_mean'),
    #                              paste0(unlist(lab_vars),'_max'),
    #                              paste0(unlist(lab_vars),'_min'))))
    # impute
    print(paste0('Cycle ',cc,', lab max/min slopes'))
    for(ll in lab_order){
      for(v in lab_vars[[ll]]){
        for(summ in c('_max_diff','_min_diff')){
          var <- paste0(v,summ)
          # missing indices from complete_data
          imiss <- is.na(complete_data_raw[[var]])
          # impute
          temp <- impute_reg(temp_data,
                             complete_data_impute[[var]],
                             imiss)
          # update complete_data_impute
          complete_data_impute[[var]] <- temp
          # print progress 
          # print(paste0('Impute ',var))
        }
      }
    }
    
    # impute lab total variation
    # temporary data for imputing
    # temp_data <- select(complete_data_impute,
    #                     all_of(c(demo_vars,smoke_vars,other_vars,
    #                              paste0(unlist(lab_vars),'_mean'),
    #                              paste0(unlist(lab_vars),'_max'),
    #                              paste0(unlist(lab_vars),'_min'),
    #                              paste0(unlist(lab_vars),'_max_diff'),
    #                              paste0(unlist(lab_vars),'_min_diff'))))
    # impute
    print(paste0('Cycle ',cc,', lab total variation'))
    for(ll in lab_order){
      for(v in lab_vars[[ll]]){
        var <- paste0(v,'_tv')
        # missing indices from complete_data
        imiss <- is.na(complete_data_raw[[var]])
        # impute
        temp <- impute_reg(temp_data,
                           complete_data_impute[[var]],
                           imiss,
                           nneg=F)
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
    imiss <- is.na(complete_data_raw[['smoke_current']])
    temp <- impute_reg_multi(temp_data,
                             complete_data_impute[,smoke_vars],
                             imiss)
    complete_data_impute[,smoke_vars] <- temp
    #print('Impute SmokeStatus')
    
    # impute demographic variables (binary) using lab means
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars_num,smoke_vars,other_vars,
                                 paste0(unlist(lab_vars),'_mean'))))
    # impute
    for(var in rev(demo_vars_ind)){
      imiss <- is.na(complete_data_raw[[var]])
      temp <- impute_reg(temp_data,
                         complete_data_impute[[var]],
                         imiss)
      complete_data_impute[[var]] <- temp
      #print(paste0('Impute ',var))
    }
    
    # impute demographic variables (numeric)
    temp_data <- select(complete_data_impute,
                        all_of(c(demo_vars_ind,smoke_vars,other_vars,
                                 paste0(unlist(lab_vars),'_mean'))))
    # impute
    for(var in demo_vars_num){
      imiss <- is.na(complete_data_raw[[var]])
      temp <- impute_reg(temp_data,
                         complete_data_impute[[var]],
                         imiss)
      complete_data_impute[[var]] <- temp
      #print(paste0('Impute ',var))
    }
    
    # check convergence
    all_miss_new <- c(complete_data_impute[is.na(complete_data_raw)])
    print(paste0('RMSE change after cycle ',cc))
    print(sqrt(mean((all_miss_new - all_miss_old)^2)))
    print(max(abs(all_miss_new - all_miss_old)))
    all_miss_old <- all_miss_new
  }
  
  # ensure smoke status is a multiclass indicator
  smoke_wmiss <- which(is.na(complete_data_raw[['smoke_current']]))
  for(ii in smoke_wmiss){
    temp <- rep(-Inf,length(smoke_vars))
    temp[which.max(c(complete_data_impute[ii,smoke_vars]))] <- max(complete_data_impute[ii,smoke_vars])
    complete_data_impute[ii,smoke_vars] <- temp  
  }

  complete_data_impfinal <- complete_data_impute
  for(var in colnames(complete_data_impfinal)){
    complete_data_impfinal[[var]] <- quantile_unnormalize(complete_data_impfinal[[var]],
                                                          complete_data_init[[var]])
  }
  
  # return complete imputed data
  return(list(clean=bind_cols(complete_data_y,complete_data_raw),
              impmed=bind_cols(complete_data_y,complete_data_med),
              impsamp=bind_cols(complete_data_y,complete_data_init),
              impreg=bind_cols(complete_data_y,complete_data_impfinal)))
}
