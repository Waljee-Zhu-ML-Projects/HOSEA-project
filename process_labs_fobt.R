# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/labs_fobt.rds')
print('load data')

# left join master to get windowing dates
table <- left_join(table,master,by='ID')
print('link prediction window')
# filter only dates in each window
table_keep <- filter(table,(labdate >= start) & (labdate <= end))
print('filter lab dates')

# lagged variables (ID, labdate and value fields)
table_keep <- table_keep %>% 
              mutate(ID_lag = lag(ID),
                     labdate_lag = lag(labdate))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$labdate_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dlabdate = ifelse(labdate - labdate_lag==0,NA,labdate - labdate_lag)) %>%
             mutate(slabdate = 1 / dlabdate)
print('calculate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats
table_FOBT <- table_group %>%
  summarise(labs_fobt_n = n(),
            labs_fobt_maxdiff = mmax(slabdate))
print('summarize FOBT')

# save result
saveRDS(table_FOBT,file='R_data/y45/labs_fobt_summary.rds')
print('save result')



