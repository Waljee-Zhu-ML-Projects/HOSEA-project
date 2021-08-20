# subsample a complete test set
source('R_code/hosea-project/utils.R')


# load complete data
sub_master <- readRDS('R_data/subsample/sub_master.rds')
sub_IDs <- sub_master$ID

complete_data <- readRDS('R_data/complete_data_raw.rds')

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

# check sample sizes
table(nonsub_complete_data$CaseControl) # 290 complete cases, leaves 710 controls
  
set.seed(300)
n2 <- 500
i_case <- which(nonsub_complete_data$CaseControl==1)
i_control <- sample(which(nonsub_complete_data$CaseControl==0),(2*n2 - length(i_case)))
  
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
other_vars <- c('colonoscopy_n','colonoscopy_maxdiff',
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
                '_maxdiff','_mindiff',
                '_tv')

# impute remaining vars with sample, save
# blood labs
set.seed(1999)
for(lab in lab_vars){
  for(v in lab){
    for(summ in lab_suffix){
      varname <- paste0(v,summ)
      complete_test[[varname]] <- fill_by_sample(complete_test[[varname]])
    }
  }
}

# demographic variables
for(varname in demo_vars){
  complete_test[[varname]] <- fill_by_sample(complete_test[[varname]])
}

# smoking indicator
smoke_miss <- is.na(complete_test[['smoke_current']])
complete_test[smoke_miss,smoke_vars] <- complete_test[sample(which(!smoke_miss),sum(smoke_miss),replace=T),smoke_vars]

saveRDS(complete_test,file='R_data/cc_complete_data.rds')
  
  
  