# source helper functions
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/process_events.R')
source('R_code/hosea-project/process_longitudinal.R')

# import master table
master <- readRDS('R_data/master.rds')

# event tables
event_tables <- c('colonoscopy',
                  'labs_fobt')
event_date_fields <- c('Procdate',
                       'labdate')

ntab <- length(event_tables)
# loop through tables
for(ii in 1:ntab){
  print(paste0('For table ',event_tables[ii],':'))
  process_events(event_tables[ii],
                 event_date_fields[ii])
}

# valued tables
value_tables <- c('labs_a1c',
                  'labs_crp')
value_date_fields <- rep('labdate',2)
value_value_fields <- list(c('A1c'),
                           c('CRP'))

ntab <- length(value_tables)
# loop through tables
for(ii in 1:ntab){
  print(paste0('For table ',value_tables[ii],':'))
  process_longitudinal(value_tables[ii],
                       value_date_fields[ii],
                       value_value_fields[[ii]])
}
