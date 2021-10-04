# script to link new charlson indicators into one big table

# source code
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_missingcharl.R')

# load master table
charlson_complete <- readRDS('R_data/master.rds')

for(charl_name in charl_names){
  print(paste0('For disease ',charl_name,':'))
  # load indicator and join
  temp <- readRDS(paste0('R_data/charlson_',charl_name,'.rds'))
  print('Load data')
  charlson_complete <- left_join(charlson_complete,temp,by='ID')
  print('Join')
  # replace
  charlson_complete[[charl_name]] <- fill_by_zero(charlson_complete[[charl_name]])
  print('Fill zeros')
}

# load n_visits for controls and cases
n_visits <- readRDS('R_data/n_visits.rds')
n_visits_case <- readRDS('R_data/n_visits_case.rds')

# bind, join tables and impute
charlson_complete <- left_join(charlson_complete,bind_rows(n_visits,n_visits_case),by='ID')
charlson_complete$n_visits <- fill_by_zero(charlson_complete$n_visits)

# save table
saveRDS(charlson_complete,file='R_data/charlson_complete_raw.rds')

# imputation with 'geometric' model 
charlson_impute <- charlson_complete

theta <- list()

for(charl_name in charl_names){
  print(paste0('For disease ',charl_name,':'))
  timestamp()
  
  # estimate for controls 
  theta[[charl_name]] <- estimate_mm(charlson_impute[[charl_name]],charlson_impute$n_visits)
  # if estimation leaves the constrained region just estimate with initializer
  if(any(theta[[charl_name]] < 0) || any(theta[[charl_name]] > 1)){
    theta[[charl_name]] <- init_mm(charlson_impute[[charl_name]],charlson_impute$n_visits,phi=TRUE)
  }
  print('Estimated')

  # impute for control
  charlson_impute[[charl_name]] <- impute_mm(charlson_impute[[charl_name]],
                                                     charlson_impute$n_visits,
                                                     theta[[charl_name]])
  print('Imputed')
}

# save imputed table
saveRDS(charlson_impute,file='R_data/charlson_complete_impute.rds')
