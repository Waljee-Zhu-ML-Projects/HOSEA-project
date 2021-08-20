# helper functions (mmean, mmin, mmax, simple imputation)
source('R_code/hosea-project/utils.R')

# begin with demographic table
sample_tab <- readRDS('R_data/sample.rds')
demo_vars <- c('ID','CaseControl',
             'ageatindex','Gender','bmi','weight',
             'Asian','Black','HawaiianPacific','IndianAlaskan',
             'SmokeStatus',
             'agentorange',
             'GerdAtIndex',
             'CHF','CTD','DEM','DIAB_C','HIV','MLD','MSLD','PARA','RD',
             'cd','copd','diab_nc','mi','pud','pvd')

complete_data <- sample_tab[,demo_vars] # 28 columns
print('load sample data')
timestamp()

# convert factors to integer and expand multiple smoking classes
# Gender
complete_data$Gender[complete_data$Gender==''] <- NA
complete_data$Gender <- as.integer(complete_data$Gender=='M')
# agentorange
complete_data$agentorange <- as.integer(complete_data$agentorange=='YES')
# SmokeStatus (keep current and former)
complete_data$smoke_current <- as.integer(complete_data$SmokeStatus==1)
complete_data$smoke_former <- as.integer(complete_data$SmokeStatus==2)
complete_data <- select(complete_data,-SmokeStatus) # adds one extra variable (29 columns)

# event tables
event_tables <- c('colonoscopy',
                  'labs_fobt')
# removed EGD and other procedures

# join summary tables (33 columns)
for(tab in event_tables){
  print(paste0('For table ',tab,':'))
  table <- readRDS(paste0('R_data/',tab,'_summary.rds'))
  complete_data <- left_join(complete_data,table,by="ID")
  print('read in and join')
  timestamp()
  rm(table)
}

# medication table (43 columns)
print('For table: allmeds')
table <- readRDS('R_data/allmeds_summary.rds')
complete_data <- left_join(complete_data,table,by="ID")
print('read in and join')
timestamp()
rm(table)

# valued tables
value_tables <- c('labs_a1c',
                  'labs_bmp',
                  'labs_cbc',
                  'labs_crp',
                  'labs_lft',
                  'labs_lipid')

# join summary tables (241 columns)
for(tab in value_tables){
  print(paste0('For table ',tab,':'))
  table <- readRDS(paste0('R_data/',tab,'_summary.rds'))
  complete_data <- left_join(complete_data,table,by="ID")
  print('read in and join')
  timestamp()
  rm(table)
}

# fill in 0's for missing event, medication data
other_vars <- c('colonoscopy_n','colonoscopy_maxdiff',
                'labs_fobt_n','labs_fobt_maxdiff',
                'h2r_int','h2r_mean','h2r_max','h2r_maxdiff','h2r_tv',
                'ppi_int','ppi_mean','ppi_max','ppi_maxdiff','ppi_tv')

for(varname in other_vars){
  complete_data[[varname]] <- fill_by_zero(complete_data[[varname]])
}
print('impute 0s for missing event data')

# save complete data (241 columns)
saveRDS(complete_data,
        file='R_data/complete_data_raw.rds')
print('save complete table')
