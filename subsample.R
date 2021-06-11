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
# IDs to keep in other tables
sub_IDs <- sub_master$ID

# rds file names without master
sample_table <- c('sample')
event_tables <- c('colonoscopy',
                  'labs_fobt')
value_tables <- c('labs_a1c',
                  'labs_bmp',
                  'labs_cbc',
                  'labs_crp',
                  'labs_lft',
                  'labs_lipid')
med_table <- 'allmeds'
rdsfiles <- paste0(c(sample_table,event_tables,value_tables,med_table),'.rds')

for(f in rdsfiles){
  temp <- readRDS(paste0('R_data/',f))
  sub_temp <- filter(temp,ID %in% sub_IDs)
  saveRDS(sub_temp,file=paste0('R_data/subsample/sub_',f))
  print(f)
}



