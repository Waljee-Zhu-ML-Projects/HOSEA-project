# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/labs_a1c.rds')
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
                     labdate_lag = lag(labdate),
                     A1c_lag = lag(A1c))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$labdate_lag[skip_rows] <- NA
table_keep$A1c_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dlabdate = pmax(1,labdate - labdate_lag),
                    dA1c = A1c - A1c_lag) %>%
             mutate(sA1c = dA1c / dlabdate)
print('calculate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats for each variable
# A1c
table_A1c <- table_group %>%
  summarise(A1c_mean = mmean(A1c),
            A1c_max = mmax(A1c),
            A1c_min = mmin(A1c),
            A1c_maxdiff = mmax(sA1c),
            A1c_mindiff = mmin(sA1c),
            A1c_tv = mmean(abs(sA1c)))
print('summarize variable A1c')

# save result
saveRDS(table_A1c,file='R_data/labs_a1c_summary.rds')
print('save result')



