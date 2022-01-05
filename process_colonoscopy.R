# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/colonoscopy.rds')
print('load data')

# left join master to get windowing dates
table <- left_join(table,master,by='ID')
print('link prediction window')
# filter only dates in each window
table_keep <- filter(table,(Procdate >= start) & (Procdate <= end))
print('filter lab dates')

# lagged variables (ID, Procdate and value fields)
table_keep <- table_keep %>% 
              mutate(ID_lag = lag(ID),
                     Procdate_lag = lag(Procdate))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$Procdate_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dProcdate = ifelse(Procdate - Procdate_lag==0,NA,Procdate - Procdate_lag)) %>%
             mutate(sProcdate = 1 / dProcdate)
print('calculate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats
table_colonoscopy <- table_group %>%
  summarise(colonoscopy_n = n(),
            colonoscopy_maxdiff = mmax(sProcdate))
print('summarize colonoscopy')

# save result
saveRDS(table_colonoscopy,file='R_data/y45/colonoscopy_summary.rds')
print('save result')



