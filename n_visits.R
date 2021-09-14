##### PART 1: import and process alldx files ##### 

# import master table
master <- readRDS('R_data/master.rds')

for(chunk in 2:5){
  print(paste0('For chunk ',chunk))
  
  # load alldxs chunk
  #timestamp()
  alldxs <- read_sas(paste0('unzipped_data/alldxs',chunk,'.sas7bdat'))
  #timestamp()
  print('Loaded')
  
  # join master
  #timestamp()
  alldxs <- left_join(alldxs,master,by='ID')
  #timestamp()
  print('Joined master table')
  
  # filter by prediction window 
  alldxs_filter <- filter(alldxs,(Dxdate >= start) & (Dxdate <= end))
  print('Filtered')
  
  # remove full table and gc
  rm(alldxs)
  gc()
  
  # save filtered data as rds 
  saveRDS(alldxs_filter,file=paste0('R_data/alldxs',chunk,'_filter.rds'))
  print('Saved')
  
  # remove and gc
  rm(alldxs_filter)
  gc()
  
}

##### PART 2: load rds files and calculate number of visits #####

n_visits_chunk <- list()
spot_check1 <- list()
spot_check2 <- list()
length(n_visits_chunk) <- length(spot_check1) <- length(spot_check2) <- 5

for(chunk in 1:5){
  print(paste0('For chunk ',chunk))
  
  # load all dx codes
  timestamp()
  alldxs <- readRDS(paste0('R_data/alldxs',chunk,'_filter.rds'))
  timestamp()
  print('Loaded')
  
  # spot check IDs 1002 and 3229242
  spot_check1[[chunk]] <- filter(alldxs,ID==1002)
  spot_check2[[chunk]] <- filter(alldxs,ID==3229242)
  # see if dates can overlap between files
  # looks like dates are mutually exclusive between files
  
  # combine into a table
  timestamp()
  n_visits_chunk[[chunk]] <- summarize(group_by(alldxs,ID),n_visits = n_distinct(Dxdate))
  timestamp()
  print('Summarized')
  
  # remove full table
  rm(alldxs)
  gc()
  
}

# combine individual files
n_visits_rbound <- map_dfr(1:5,function(ii){n_visits_chunk[[ii]]})
print('Bind individual files')

n_visits <- summarize(group_by(n_visits_rbound,ID),n_visits = sum(n_visits))
print('Sum over chunks')

# save as an rds file
saveRDS(n_visits,file='R_data/n_visits.rds')
print('Saved')

# most likely is around 7-10 visits, some people have more, some very extreme cases
# one person has 1431/1461 days with a code



