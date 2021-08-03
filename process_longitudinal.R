# # working directory
# setwd('~')
#
# # libraries
# library(dplyr)
#
# # names
# table_name <- 'labs_crp'
# date_field <- 'labdate'
# value_fields <- c('CRP')
#
# # names
# table_name <- 'labs_lft'
# date_field <- 'labdate'
# value_fields <- c('alkphos','alt','ast','totprot')

# pre-imported table
# master <- readRDS('R_data/subsample/sub_master.rds')

process_longitudinal <- function(table_name,date_field,value_fields){
  
  # import current table
  table <- readRDS(paste0('R_data/',table_name,'.rds'))
  #table <- readRDS(paste0('R_data/subsample/sub_',table_name,'.rds'))
  print('load data')
  
  # left join master to get windowing dates
  table <- left_join(table,master,by='ID')
  print('link prediction window')
  # filter only dates in each window
  table_keep <- filter(table,(!!as.symbol(date_field) >= start) & (!!as.symbol(date_field) <= end))
  print('filter lab dates')
  # trivial amount of time to here
  
  # initialize table
  table_summary <- distinct(table_keep[,c('ID')])
  for(var_name in value_fields){
    print(paste0('for variable ',var_name,':'))
    # subset columns
    table_temp <- table_keep[,c('ID',date_field,var_name)]
    # function to process longitudinal stats
    # redefine for current (table,variable) pair
    long_stats <- function(tab,key){
      temp <- as.data.frame(c(mmean(tab[[var_name]])))
      colnames(temp) <- paste0(var_name,'_mean')
      if(nrow(tab) > 1){
        # max/min predictors
        temp[[paste0(var_name,'_max')]] <- mmax(tab[[var_name]])
        temp[[paste0(var_name,'_min')]] <- mmin(tab[[var_name]])
        # max/min slope predictors
        slopes <- diff(tab[[var_name]])/pmax(diff(tab[[date_field]]),1)
        temp[[paste0(var_name,'_maxdiff')]] <- mmax(slopes)
        temp[[paste0(var_name,'_mindiff')]] <- mmin(slopes)
        if(nrow(tab) > 2){
          # total variation predictor
          temp[[paste0(var_name,'_tv')]] <- mmean(abs(slopes))
        }
        else{
          temp[[paste0(var_name,'_tv')]] <- NA
        }
      }
      else{
        temp[[paste0(var_name,'_max')]] <- NA
        temp[[paste0(var_name,'_min')]] <- NA
        temp[[paste0(var_name,'_maxdiff')]] <- NA
        temp[[paste0(var_name,'_mindiff')]] <- NA
        temp[[paste0(var_name,'_tv')]] <- NA
      }
      return(temp)
    }
    # summarize grouped table
    summary_temp <- group_modify(group_by(table_temp,ID),long_stats)
    print('group and summarize')
    # link to current summary table
    table_summary <- left_join(table_summary,summary_temp,by='ID')
    print('link to existing summary')
  }
  
  # save summary table
  saveRDS(table_summary,
          file=paste0('R_data/',table_name,'_summary.rds'))
  #saveRDS(table_summary,
  #        file=paste0('R_data/subsample/sub_',table_name,'_summary.rds'))
  print('save result')
  
  # garbage collection
  gc()
}
