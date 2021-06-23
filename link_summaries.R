# helper functions (mmean, mmin, mmax, simple imputation)
source('R_code/hosea-project/utils.R')

# begin with demographic table
sample_tab <- readRDS('R_data/subsample/sub_sample.rds')
demo_vars <- c('ID','CaseControl',
             'ageatindex','Gender','bmi','weight',
             'Asian','Black','HawaiianPacific','IndianAlaskan',
             'SmokeStatus',
             'agentorange',
             'GerdAtIndex',
             'CHF','CTD','DEM','DIAB_C','HIV','MLD','MSLD','PARA','RD',
             'cd','copd','diab_nc','mi','pud','pvd')

complete_data <- sample_tab[,demo_vars] # 28 columns

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
  table <- readRDS(paste0('R_data/subsample/sub_',tab,'_summary.rds'))
  complete_data <- left_join(complete_data,table,by="ID")
  print('read in and join')
  rm(table)
}

# medication table (43 columns)
print('For table: allmeds')
table <- readRDS('R_data/subsample/sub_allmeds_summary.rds')
complete_data <- left_join(complete_data,table,by="ID")
print('read in and join')
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
  table <- readRDS(paste0('R_data/subsample/sub_',tab,'_summary.rds'))
  complete_data <- left_join(complete_data,table,by="ID")
  print('read in and join')
  rm(table)
}

# save complete data (241 columns)
saveRDS(complete_data,
        file='R_data/subsample/sub_complete_data_raw.rds')
