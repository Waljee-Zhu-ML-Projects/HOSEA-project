# analysis of n_visits
source('R_code/hosea-project/utils.R')
source('R_code/hosea-project/utils_missingcharl.R')

# import data
complete_data <- readRDS('R_data/complete_data_raw.rds')
complete_data$n_visits = NULL
n_visits <- readRDS('R_data/n_visits.rds')
n_visits_case <- readRDS('R_data/n_visits_case.rds')

# join tables and impute
complete_data <- left_join(complete_data,bind_rows(n_visits,n_visits_case),by='ID')
complete_data$n_visits <- fill_by_zero(complete_data$n_visits)

# histogram
hist(complete_data$n_visits,1000,probability=F,main='Histogram of number of unique visit dates')

# compare means
# controls
median(complete_data$n_visits[complete_data$CaseControl==0]) 
# cases
median(complete_data$n_visits[complete_data$CaseControl==1]) # cases have more visits in both mean and median?

# proportion of observed Charlson codes by n_visits
prop_obs <- complete_data %>% group_by(n_visits) %>% summarize(obs = mean(!is.na(CHF)))
# also get sample sizes
prop_ss <- complete_data %>% group_by(n_visits) %>% summarize(ss = n())


# plot(as.matrix(prop_obs),
#      main='Proportion of observed Charlson scores',
#      xlab='n_visits',
#      ylab='P(not NA)')
# 
# # restricted to < 200, very nice pattern among controls
# plot(as.matrix(prop_obs)[1:401,],
#      main='Proportion of observed Charlson scores',
#      xlab='n_visits',
#      ylab='P(not NA)')

# now look at cases and controls
prop_obs_control <- complete_data %>% filter(CaseControl==0) %>%  group_by(n_visits) %>% summarize(obs = mean(!is.na(CHF)))
prop_obs_case <- complete_data %>% filter(CaseControl==1) %>% group_by(n_visits) %>% summarize(obs = mean(!is.na(CHF)))
# also get sample sizes
prop_ss_control <- complete_data %>% filter(CaseControl==0) %>% group_by(n_visits) %>% summarize(ss = n())
prop_ss_case <- complete_data %>% filter(CaseControl==1) %>% group_by(n_visits) %>% summarize(ss = n())

plot(as.matrix(prop_obs_control),
     main='Proportion of observed Charlson scores',
     xlab='n_visits',
     ylab='P(not NA)')
lines(as.matrix(prop_obs_case),type='p',col='red')

par(mfrow=c(2,2))
plot(as.matrix(prop_ss_control),
     main="Sample sizes: controls",
     xlab='n_visits',
     ylab='sample size',type='l')
plot(as.matrix(prop_ss_case),
     type='l',col='red',
     main="Sample sizes: cases",
     xlab='n_visits',
     ylab='sample size')
plot(as.matrix(prop_ss_control)[1:201,],
     main="Sample sizes: controls (<=200 visits)",
     xlab='n_visits',
     ylab='sample size',type='l')
plot(as.matrix(prop_ss_case)[1:201,],
     type='l',col='red',
     main="Sample sizes: cases (<= 200 visits)",
     xlab='n_visits',
     ylab='sample size')

# restricted to n_visits <= 200 to see the pattern?
par(mfrow=c(1,1))
plot(as.matrix(prop_obs_control)[1:201,],
     main='Proportion of observed Charlson scores',
     xlab='n_visits',
     ylab='P(not NA)',
     ylim=c(0,1))
lines(as.matrix(prop_obs_case)[1:201,],type='p',col='red')

# same code treating 0/NA as indistinguishable
# need to pick a specific disease, here COPD, could do this for any code
# prop_one <- complete_data %>% group_by(n_visits) %>% summarize(obs = mean(copd %in% 1))
# 
# plot(as.matrix(prop_one),
#      main='Proportion with COPD=1',
#      xlab='n_visits',
#      ylab='P(COPD=1)')
# 
# plot(as.matrix(prop_one)[1:401,],
#      main='Proportion with COPD=1',
#      xlab='n_visits',
#      ylab='P(COPD=1)')

# code again split by control/case
prop_one_control <- complete_data %>% filter(CaseControl==0) %>% 
   group_by(n_visits) %>% summarize(obs = mean(copd %in% 1))
prop_one_case <- complete_data %>% filter(CaseControl==1) %>% 
   group_by(n_visits) %>% summarize(obs = mean(copd %in% 1))

plot(as.matrix(prop_one_control),
     main='Proportion with COPD=1',
     xlab='n_visits',
     ylab='P(COPD=1)')
lines(as.matrix(prop_one_case),type='p',col='red')

# restict to n_visits <= 200 to see pattern
plot(as.matrix(prop_one_control)[1:201,],
     main='Proportion with COPD=1',
     xlab='n_visits',
     ylab='P(COPD=1)',
     ylim=c(0,1))
lines(as.matrix(prop_one_case)[1:201,],type='p',col='red')

# pooled estimation for COPD
theta_copd <- estimate_mm(as.integer(complete_data$copd %in% 1),complete_data$n_visits)

# plot the estimated curve for n_visits <= 200
V <- 0:200
pV <- theta_copd[1]*(1-(1-theta_copd[2])^V)
lines(V,pV,col='blue',type='l',lwd=2)

#### icd9 vs icd10 ####
# additional analysis

# import data
complete_data <- readRDS('R_data/complete_data_raw.rds')
n_visits_case_icd9 <- readRDS('R_data/n_visits_case_icd9.rds')
n_visits_case_icd10 <- readRDS('R_data/n_visits_case_icd10.rds')

complete_data <- left_join(filter(complete_data,CaseControl==1),n_visits_case_icd9,by='ID')
complete_data <- left_join(complete_data,n_visits_case_icd10,by='ID')
complete_data$n_visits.x <- fill_by_zero(complete_data$n_visits.x)
complete_data$n_visits.y <- fill_by_zero(complete_data$n_visits.y)

# scatter plot of ICD9 vs ICD10 nvisits
plot(complete_data$n_visits.x,complete_data$n_visits.y,
     xlab='ICD9 visits',ylab='ICD10 visits')

# missing charlson analysis/scatter plot with just ICD9
prop_obs9 <- complete_data %>% group_by(n_visits.x) %>% summarize(obs = mean(!is.na(CHF)))
# also get sample sizes
prop_ss9 <- complete_data %>% group_by(n_visits.x) %>% summarize(ss = n())

# plot(as.matrix(prop_obs9),
#      main='Proportion of observed Charlson scores, ICD9',
#      xlab='n_visits ICD9',
#      ylab='P(not NA)',col='blue')

# missing charlson analysis/scatter plot with just ICD10
prop_obs10 <- complete_data %>% group_by(n_visits.y) %>% summarize(obs = mean(!is.na(CHF)))
# also get sample sizes
prop_ss10 <- complete_data %>% group_by(n_visits.y) %>% summarize(ss = n())

# plot(as.matrix(prop_obs10),
#      main='Proportion of observed Charlson scores, ICD10',
#      xlab='n_visits ICD10',
#      ylab='P(not NA)',col='red')

# plot both
plot(as.matrix(prop_obs9),
     main='Proportion of observed Charlson scores, blue:ICD9, red:ICD10',
     xlab='n_visits',
     ylab='P(not NA)',col='blue')
lines(as.matrix(prop_obs10),
      type='p',col='red')


