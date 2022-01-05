# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/labs_bmp.rds')
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
                     bun_lag = lag(bun),
                     calc_lag = lag(calc),
                     chlor_lag = lag(chlor),
                     co2_lag = lag(co2),
                     creat_lag = lag(creat),
                     gluc_lag = lag(gluc),
                     k_lag = lag(k),
                     na_lag = lag(na))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$labdate_lag[skip_rows] <- NA
table_keep$bun_lag[skip_rows] <- NA
table_keep$calc_lag[skip_rows] <- NA
table_keep$chlor_lag[skip_rows] <- NA
table_keep$co2_lag[skip_rows] <- NA
table_keep$creat_lag[skip_rows] <- NA
table_keep$gluc_lag[skip_rows] <- NA
table_keep$k_lag[skip_rows] <- NA
table_keep$na_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dlabdate = pmax(1,labdate - labdate_lag),
                    dbun = bun - bun_lag,
                    dcalc = calc - calc_lag,
                    dchlor = chlor - chlor_lag,
                    dco2 = co2 - co2_lag,
                    dcreat = creat - creat_lag,
                    dgluc = gluc - gluc_lag,
                    dk = k - k_lag,
                    dna = na - na_lag) %>%
             mutate(sbun = dbun / dlabdate,
                    scalc = dcalc / dlabdate,
                    schlor = dchlor / dlabdate,
                    sco2 = dco2 / dlabdate,
                    screat = dcreat / dlabdate,
                    sgluc = dgluc / dlabdate,
                    sk = dk / dlabdate,
                    sna = dna / dlabdate)
print('calculate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats for each variable
# bun
table_bun <- table_group %>%
  summarise(bun_mean = mmean(bun),
            bun_max = mmax(bun),
            bun_min = mmin(bun),
            bun_maxdiff = mmax(sbun),
            bun_mindiff = mmin(sbun),
            bun_tv = mmean(abs(sbun)))
print('summarize variable bun')

# calc
table_calc <- table_group %>%
  summarise(calc_mean = mmean(calc),
            calc_max = mmax(calc),
            calc_min = mmin(calc),
            calc_maxdiff = mmax(scalc),
            calc_mindiff = mmin(scalc),
            calc_tv = mmean(abs(scalc)))
print('summarize variable calc')

# chlor
table_chlor <- table_group %>%
  summarise(chlor_mean = mmean(chlor),
            chlor_max = mmax(chlor),
            chlor_min = mmin(chlor),
            chlor_maxdiff = mmax(schlor),
            chlor_mindiff = mmin(schlor),
            chlor_tv = mmean(abs(schlor)))
print('summarize variable chlor')


# co2
table_co2 <- table_group %>%
  summarise(co2_mean = mmean(co2),
            co2_max = mmax(co2),
            co2_min = mmin(co2),
            co2_maxdiff = mmax(sco2),
            co2_mindiff = mmin(sco2),
            co2_tv = mmean(abs(sco2)))
print('summarize variable co2')

# creat
table_creat <- table_group %>%
  summarise(creat_mean = mmean(creat),
            creat_max = mmax(creat),
            creat_min = mmin(creat),
            creat_maxdiff = mmax(screat),
            creat_mindiff = mmin(screat),
            creat_tv = mmean(abs(screat)))
print('summarize variable creat')

# gluc
table_gluc <- table_group %>%
  summarise(gluc_mean = mmean(gluc),
            gluc_max = mmax(gluc),
            gluc_min = mmin(gluc),
            gluc_maxdiff = mmax(sgluc),
            gluc_mindiff = mmin(sgluc),
            gluc_tv = mmean(abs(sgluc)))
print('summarize variable gluc')

# k
table_k <- table_group %>%
  summarise(k_mean = mmean(k),
            k_max = mmax(k),
            k_min = mmin(k),
            k_maxdiff = mmax(sk),
            k_mindiff = mmin(sk),
            k_tv = mmean(abs(sk)))
print('summarize variable k')

# na
table_na <- table_group %>%
  summarise(na_mean = mmean(na),
            na_max = mmax(na),
            na_min = mmin(na),
            na_maxdiff = mmax(sna),
            na_mindiff = mmin(sna),
            na_tv = mmean(abs(sna)))
print('summarize variable na')

# bind columns (remove ID columns)
table_summary <- bind_cols(table_bun,
                           table_calc[,-1],
                           table_chlor[,-1],
                           table_co2[,-1],
                           table_creat[,-1],
                           table_gluc[,-1],
                           table_k[,-1],
                           table_na[,-1])
print('bind tables')

# save result
saveRDS(table_summary,file='R_data/y45/labs_bmp_summary.rds')
print('save result')



