# # working directory
# setwd('~')
#
# # libraries
# library(dplyr)
#
# # names
# table_name <- 'labs_fobt'
# date_field <- 'labdate'
# # pre-imported tables
# master <- readRDS('R_data/subsample/sub_master.rds')

process_events <- function(table_name,date_field){
  # import current data
  table <- readRDS(paste0('R_data/',table_name,'.rds'))
  #table <- readRDS(paste0('R_data/subsample/sub_',table_name,'.rds'))
  
  print('load data')

  # only for event data
  table <- distinct(table[,c('ID',date_field)])

  # left join master to get windowing dates
  table <- left_join(table,master,by='ID')
  print('link prediction window')
  # filter only dates in each window
  table_keep <- filter(table,(!!as.symbol(date_field) >= start) & (!!as.symbol(date_field) <= end))
  print('filter lab dates')

  # helper to process groups for number and min diff
  # note diff assumes entries are already ordered by date within each ID
  event_stats <- function(tab,key){
    temp <- as.data.frame(c(nrow(tab)))
    colnames(temp) <- paste0(table_name,'_n')
    if(nrow(tab) > 1){
      diffs <- 1/pmax(diff(tab[[date_field]]),1)
      temp[[paste0(table_name,'_maxdiff')]] <- max(diffs)
    }
    else{
      temp[[paste0(table_name,'_maxdiff')]] <- NA
    }
    return(temp)
  }

  # summarize with count and most recent
  table_summary <- group_modify(group_by(table_keep,ID),event_stats)
  print('group and summarize')

  # save result
  saveRDS(table_summary,
          file=paste0('R_data/',table_name,'_summary.rds'))
  #saveRDS(table_summary,
  #        file=paste0('R_data/subsample/sub_',table_name,'_summary.rds'))
  print('save result')
  
  # garbage collection
  gc()
}

