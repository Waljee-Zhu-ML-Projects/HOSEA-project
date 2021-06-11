# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/labs_lipid.rds')
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
                     chol_lag = lag(chol),
                     hdl_lag = lag(hdl),
                     ldl_lag = lag(ldl),
                     trig_lag = lag(trig))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$labdate_lag[skip_rows] <- NA
table_keep$chol_lag[skip_rows] <- NA
table_keep$hdl_lag[skip_rows] <- NA
table_keep$ldl_lag[skip_rows] <- NA
table_keep$trig_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dlabdate = pmax(1,labdate - labdate_lag),
                    dchol = chol - chol_lag,
                    dhdl = hdl - hdl_lag,
                    dldl = ldl - ldl_lag,
                    dtrig = trig - trig_lag) %>%
             mutate(schol = dchol / dlabdate,
                    shdl = dhdl / dlabdate,
                    sldl = dldl / dlabdate,
                    strig = dtrig / dlabdate)
print('calculate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats for each variable
# chol
table_chol <- table_group %>%
  summarise(chol_mean = mmean(chol),
            chol_max = mmax(chol),
            chol_min = mmin(chol),
            chol_maxdiff = mmax(schol),
            chol_mindiff = mmin(schol),
            chol_tv = mmean(abs(schol)))
print('summarize variable chol')

# hdl
table_hdl <- table_group %>%
  summarise(hdl_mean = mmean(hdl),
            hdl_max = mmax(hdl),
            hdl_min = mmin(hdl),
            hdl_maxdiff = mmax(shdl),
            hdl_mindiff = mmin(shdl),
            hdl_tv = mmean(abs(shdl)))
print('summarize variable hdl')

# ldl
table_ldl <- table_group %>%
  summarise(ldl_mean = mmean(ldl),
            ldl_max = mmax(ldl),
            ldl_min = mmin(ldl),
            ldl_maxdiff = mmax(sldl),
            ldl_mindiff = mmin(sldl),
            ldl_tv = mmean(abs(sldl)))
print('summarize variable ldl')


# trig
table_trig <- table_group %>%
  summarise(trig_mean = mmean(trig),
            trig_max = mmax(trig),
            trig_min = mmin(trig),
            trig_maxdiff = mmax(strig),
            trig_mindiff = mmin(strig),
            trig_tv = mmean(abs(strig)))
print('summarize variable trig')

# bind columns (remove ID columns)
table_summary <- bind_cols(table_chol,
                           table_hdl[,-1],
                           table_ldl[,-1],
                           table_trig[,-1])
print('bind tables')

# save result
saveRDS(table_summary,file='R_data/labs_lipid_summary.rds')
print('save result')



