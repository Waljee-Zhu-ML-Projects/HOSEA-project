# subsample the master table and all the other tables with the same 
# IDs
# random sample of n2 = 5000 cases and controls

# load in the master table
master <- readRDS('R_data/master.rds')

# get the subsample of IDs
set.seed(200)
n2 <- 5000
i_case <- sample(which(master$case),n2)
i_control <- sample(which(!master$case),n2)
sub_master <- master[c(i_case,i_control),]
# save this table
saveRDS(sub_master,file='R_data/subsample/sub_master.rds')
# IDs to keep
sub_IDs <- sub_master$ID
print('sampled IDs')

temp <- readRDS(paste0('R_data/complete_data_raw.rds'))
sub_temp <- filter(temp,ID %in% sub_IDs)
saveRDS(sub_temp,file='R_data/subsample/sub_complete_data_raw')
print('saved subsampled data')

