# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/labs_crp.rds')
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
                     CRP_lag = lag(CRP))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$labdate_lag[skip_rows] <- NA
table_keep$CRP_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dlabdate = pmax(1,labdate - labdate_lag),
                    dCRP = CRP - CRP_lag) %>%
             mutate(sCRP = dCRP / dlabdate)
print('calculate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats for each variable
# CRP
table_CRP <- table_group %>%
  summarise(CRP_mean = mmean(CRP),
            CRP_max = mmax(CRP),
            CRP_min = mmin(CRP),
            CRP_maxdiff = mmax(sCRP),
            CRP_mindiff = mmin(sCRP),
            CRP_tv = mmean(abs(sCRP)))
print('summarize variable CRP')

# save result
saveRDS(table_CRP,file='R_data/labs_crp_summary.rds')
print('save result')



