# analysis of n_visits
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_missingcharl.R')

# import data
complete_data <- readRDS('R_data/complete_data_raw.rds')
n_visits <- readRDS('R_data/n_visits.rds')

# join tables and impute
complete_data <- left_join(filter(complete_data,CaseControl==0),n_visits,by='ID')
complete_data$n_visits <- fill_by_zero(complete_data$n_visits)

# histogram
hist(complete_data$n_visits,1000,probability=F,main='Histogram of number of unique visit dates')

# proportion of observed Charlson codes by n_visits
prop_obs <- complete_data %>% group_by(n_visits) %>% summarize(obs = mean(!is.na(CHF)))

plot(as.matrix(prop_obs),
     main='Proportion of observed Charlson scores',
     xlab='n_visits',
     ylab='P(not NA)')

# restricted to < 200, very nice pattern
plot(as.matrix(prop_obs)[1:401,],
     main='Proportion of observed Charlson scores',
     xlab='n_visits',
     ylab='P(not NA)')

# same code treating 0/NA as indistinguishable
# need to pick a specific disease, here COPD, could do this for any code
prop_one <- complete_data %>% group_by(n_visits) %>% summarize(obs = mean(copd %in% 1))

plot(as.matrix(prop_one),
     main='Proportion with COPD=1',
     xlab='n_visits',
     ylab='P(COPD=1)')

plot(as.matrix(prop_one)[1:401,],
     main='Proportion with COPD=1',
     xlab='n_visits',
     ylab='P(COPD=1)')

# estimation for COPD
theta_copd <- estimate_mm(as.integer(complete_data$copd %in% 1),complete_data$n_visits)

# plot the estimated curve
V <- 0:400
pV <- theta_copd[1]*(1-(1-theta_copd[2])^V)
lines(V,pV,col='red',type='l',lwd=2)

