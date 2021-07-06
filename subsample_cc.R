# subsample a complete test set

# load complete data
sub_master <- readRDS('R_data/subsample/sub_master.rds')
sub_IDs <- sub_master$ID

complete_data <- readRDS('R_data/complete_data_clean.rds')

nonsub_complete_data <- complete_data %>% 
  filter(!(ID %in% sub_IDs)) %>%
  filter(!is.na(CHF)) %>%
  filter(!is.na(A1c_mean)) %>%
  filter(!is.na(alkphos_mean)) %>%
  filter(!is.na(baso_mean)) %>%
  filter(!is.na(chol_mean)) %>%
  filter(!is.na(bun_mean)) %>%
  filter(!is.na(smoke_current)) %>%
  filter(!is.na(bmi)) %>%
  filter(!is.na(Black))
  
set.seed(300)
n2 <- 5000
i_case <- sample(which(nonsub_complete_data$CaseControl==1),n2)
i_control <- sample(which(nonsub_complete_data$CaseControl==0),n2)
  
complete_test <- nonsub_complete_data[c(i_case,i_control),]

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

# impute colonoscopies, fobt labs and meds with zeros
for(var in other_vars){
  complete_test[[var]] <- fill_by_zero(complete_test[[var]])
  #print(paste0('Impute ',var))
}

# impute remaining vars with sample, save

# blood labs
set.seed(1999)
for(lab in lab_vars){
  for(v in lab){
    for(summ in lab_suffix){
      varname <- paste0(v,summ)
      complete_test[[varname]] <- fill_by_sample(complete_test[[varname]])
      #print(paste0('Impute ',varname))
    }
  }
}

# demographic variables
for(varname in demo_vars){
  complete_test[[varname]] <- fill_by_sample(complete_test[[varname]])
  #print(paste0('Impute ',varname))
}

# smoking indicator
smoke_miss <- is.na(complete_test[['smoke_current']])
complete_test[smoke_miss,smoke_vars] <- complete_test[sample(which(!smoke_miss),sum(smoke_miss),replace=T),smoke_vars]
#print('Impute smoking status')

saveRDS(complete_test,file='R_data/subsample/cc_test.rds')
  
  
  