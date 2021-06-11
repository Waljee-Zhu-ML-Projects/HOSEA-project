# code to create a master table with IDs, case/control
# and start/end dates for prediction window

# load sample 
sample <- readRDS('R_data/sample.rds') # (loads in much less time)

master <- sample[,c('ID')]

master$case <- !is.na(sample$datedx)

# calculate a prediction window in days
master$start <- sample$IndexDate - (2*365 + 364)
master$end <- sample$IndexDate - 364

# save the table
saveRDS(master,file='R_data/master.rds')