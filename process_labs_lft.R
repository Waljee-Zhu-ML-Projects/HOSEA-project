# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/labs_lft.rds')
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
                     alkphos_lag = lag(alkphos),
                     alt_lag = lag(alt),
                     ast_lag = lag(ast),
                     totprot_lag = lag(totprot))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$labdate_lag[skip_rows] <- NA
table_keep$alkphos_lag[skip_rows] <- NA
table_keep$alt_lag[skip_rows] <- NA
table_keep$ast_lag[skip_rows] <- NA
table_keep$totprot_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dlabdate = pmax(1,labdate - labdate_lag),
                    dalkphos = alkphos - alkphos_lag,
                    dalt = alt - alt_lag,
                    dast = ast - ast_lag,
                    dtotprot = totprot - totprot_lag) %>%
             mutate(salkphos = dalkphos / dlabdate,
                    salt = dalt / dlabdate,
                    sast = dast / dlabdate,
                    stotprot = dtotprot / dlabdate)
print('calculate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats for each variable
# alkphos
table_alkphos <- table_group %>%
  summarise(alkphos_mean = mmean(alkphos),
            alkphos_max = mmax(alkphos),
            alkphos_min = mmin(alkphos),
            alkphos_maxdiff = mmax(salkphos),
            alkphos_mindiff = mmin(salkphos),
            alkphos_tv = mmean(abs(salkphos)))
print('summarize variable alkphos')

# alt
table_alt <- table_group %>%
  summarise(alt_mean = mmean(alt),
            alt_max = mmax(alt),
            alt_min = mmin(alt),
            alt_maxdiff = mmax(salt),
            alt_mindiff = mmin(salt),
            alt_tv = mmean(abs(salt)))
print('summarize variable alt')


# ast
table_ast <- table_group %>%
  summarise(ast_mean = mmean(ast),
            ast_max = mmax(ast),
            ast_min = mmin(ast),
            ast_maxdiff = mmax(sast),
            ast_mindiff = mmin(sast),
            ast_tv = mmean(abs(sast)))
print('summarize variable ast')


# totprot
table_totprot <- table_group %>%
  summarise(totprot_mean = mmean(totprot),
            totprot_max = mmax(totprot),
            totprot_min = mmin(totprot),
            totprot_maxdiff = mmax(stotprot),
            totprot_mindiff = mmin(stotprot),
            totprot_tv = mmean(abs(stotprot)))
print('summarize variable totprot')


# bind columns (remove ID columns)
table_summary <- bind_cols(table_alkphos,
                           table_alt[,-1],
                           table_ast[,-1],
                           table_totprot[,-1])
print('bind tables')

# save result
saveRDS(table_summary,file='R_data/labs_lft_summary.rds')
print('save result')



