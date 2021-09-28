# script to link new charlson indicators into one big table

# source code
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_missingcharl.R')

# load master table
charlson_complete <- readRDS('R_data/master.rds')

for(charl_name in charl_names){
  # load indicator and join
  temp <- readRDS(paste0('R_data/charlson_',charl_name,'.rds'))
  charlson_complete <- left_join(charlson_complete,temp,by='ID')
  # replace
  charlson_complete[[charl_name]] <- fill_by_zero(charlson_complete[[charl_name]])
}

# load n_visits for controls and cases
n_visits <- readRDS('R_data/n_visits.rds')
n_visits_case <- readRDS('R_data/n_visits_case.rds')

# bind, join tables and impute
charlson_complete <- left_join(charlson_complete,bind_rows(n_visits,n_visits_case),by='ID')
charlson_complete$n_visits <- fill_by_zero(charlson_complete$n_visits)

# save table
saveRDS(charlson_complete,file='R_data/charlson_complete_raw.rds')

# next step:
# subset to cases and controls
# estimate missingness parameters for cases/controls separately for each indicator (list 
# of thetas for cases and controls)
# copy the tables and impute with expected value
# bind rows and save as an imputed version



