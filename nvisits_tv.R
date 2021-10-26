# import data
complete_data <- readRDS('R_data/complete_data_raw.rds')
complete_data$n_visits = NULL
n_visits <- readRDS('R_data/n_visits.rds')
n_visits_case <- readRDS('R_data/n_visits_case.rds')
# join tables and impute
complete_data <- left_join(complete_data,bind_rows(n_visits,n_visits_case),by='ID')
complete_data$n_visits <- fill_by_zero(complete_data$n_visits)

library(ggplot2)

ggplot(complete_data, aes(log(n_visits), log(chol_tv))) + 
  geom_bin2d(bins=100) + theme_bw() 
ggplot(complete_data, aes(log(n_visits), log(hdl_tv))) + 
  geom_bin2d(bins=100) + theme_bw() 
ggplot(complete_data, aes(log(n_visits), log(trig_tv))) + 
  geom_bin2d(bins=100) + theme_bw() 
