# helper functions
source('R_code/hosea-project/utils.R')

# load master
master <- readRDS('R_data/master.rds')

# load table
table <- readRDS('R_data/labs_cbc.rds')
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
                     baso_lag = lag(baso),
                     eos_lag = lag(eos),
                     hct_lag = lag(hct),
                     hgb_lag = lag(hgb),
                     lymph_lag = lag(lymph),
                     mch_lag = lag(mch),
                     mchc_lag = lag(mchc),
                     mcv_lag = lag(mcv),
                     mono_lag = lag(mono),
                     mpv_lag = lag(mpv),
                     neut_lag = lag(neut),
                     platelet_lag = lag(platelet),
                     rbc_lag = lag(rbc),
                     rdw_lag = lag(rdw),
                     wbc_lag = lag(wbc))
print('compute lagged variables')

# find where IDs change
skip_rows <- !(table_keep$ID == table_keep$ID_lag)
# make lagged values NAs
table_keep$ID_lag[skip_rows] <- NA
table_keep$labdate_lag[skip_rows] <- NA
table_keep$baso_lag[skip_rows] <- NA
table_keep$eos_lag[skip_rows] <- NA
table_keep$hct_lag[skip_rows] <- NA
table_keep$hgb_lag[skip_rows] <- NA
table_keep$lymph_lag[skip_rows] <- NA
table_keep$mch_lag[skip_rows] <- NA
table_keep$mchc_lag[skip_rows] <- NA
table_keep$mcv_lag[skip_rows] <- NA
table_keep$mono_lag[skip_rows] <- NA
table_keep$mpv_lag[skip_rows] <- NA
table_keep$neut_lag[skip_rows] <- NA
table_keep$platelet_lag[skip_rows] <- NA
table_keep$rbc_lag[skip_rows] <- NA
table_keep$rdw_lag[skip_rows] <- NA
table_keep$wbc_lag[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
table_all <- table_keep %>%
             mutate(dlabdate = pmax(1,labdate - labdate_lag),
                    dbaso = baso - baso_lag,
                    deos = eos - eos_lag,
                    dhct = hct - hct_lag,
                    dhgb = hgb - hgb_lag,
                    dlymph = lymph - lymph_lag,
                    dmch = mch - mch_lag,
                    dmchc = mchc - mchc_lag,
                    dmcv = mcv - mcv_lag,
                    dmono = mono - mono_lag,
                    dmpv = mpv - mpv_lag,
                    dneut = neut - neut_lag,
                    dplatelet = platelet - platelet_lag,
                    drbc = rbc - rbc_lag,
                    drdw = rdw - rdw_lag,
                    dwbc = wbc - wbc_lag) %>%
             mutate(sbaso = dbaso / dlabdate,
                    seos = deos / dlabdate,
                    shct = dhct / dlabdate,
                    shgb = dhgb / dlabdate,
                    slymph = dlymph / dlabdate,
                    smch = dmch / dlabdate,
                    smchc = dmchc / dlabdate,
                    smcv = dmcv / dlabdate,
                    smono = dmono / dlabdate,
                    smpv = dmpv / dlabdate,
                    sneut = dneut / dlabdate,
                    splatelet = dplatelet / dlabdate,
                    srbc = drbc / dlabdate,
                    srdw = drdw / dlabdate,
                    swbc = dwbc / dlabdate)
print('eosulate slopes')

# group by ID
table_group <- table_all %>%
                 group_by(ID) 
print('group by ID')

# summary stats for each variable
# baso
table_baso <- table_group %>%
  summarise(baso_mean = mmean(baso),
            baso_max = mmax(baso),
            baso_min = mmin(baso),
            baso_maxdiff = mmax(sbaso),
            baso_mindiff = mmin(sbaso),
            baso_tv = mmean(abs(sbaso)))
print('summarize variable baso')

# eos
table_eos <- table_group %>%
  summarise(eos_mean = mmean(eos),
            eos_max = mmax(eos),
            eos_min = mmin(eos),
            eos_maxdiff = mmax(seos),
            eos_mindiff = mmin(seos),
            eos_tv = mmean(abs(seos)))
print('summarize variable eos')

# hct
table_hct <- table_group %>%
  summarise(hct_mean = mmean(hct),
            hct_max = mmax(hct),
            hct_min = mmin(hct),
            hct_maxdiff = mmax(shct),
            hct_mindiff = mmin(shct),
            hct_tv = mmean(abs(shct)))
print('summarize variable hct')


# hgb
table_hgb <- table_group %>%
  summarise(hgb_mean = mmean(hgb),
            hgb_max = mmax(hgb),
            hgb_min = mmin(hgb),
            hgb_maxdiff = mmax(shgb),
            hgb_mindiff = mmin(shgb),
            hgb_tv = mmean(abs(shgb)))
print('summarize variable hgb')

# lymph
table_lymph <- table_group %>%
  summarise(lymph_mean = mmean(lymph),
            lymph_max = mmax(lymph),
            lymph_min = mmin(lymph),
            lymph_maxdiff = mmax(slymph),
            lymph_mindiff = mmin(slymph),
            lymph_tv = mmean(abs(slymph)))
print('summarize variable lymph')

# mch
table_mch <- table_group %>%
  summarise(mch_mean = mmean(mch),
            mch_max = mmax(mch),
            mch_min = mmin(mch),
            mch_maxdiff = mmax(smch),
            mch_mindiff = mmin(smch),
            mch_tv = mmean(abs(smch)))
print('summarize variable mch')

# mchc
table_mchc <- table_group %>%
  summarise(mchc_mean = mmean(mchc),
            mchc_max = mmax(mchc),
            mchc_min = mmin(mchc),
            mchc_maxdiff = mmax(smchc),
            mchc_mindiff = mmin(smchc),
            mchc_tv = mmean(abs(smchc)))
print('summarize variable mchc')

# mcv
table_mcv <- table_group %>%
  summarise(mcv_mean = mmean(mcv),
            mcv_max = mmax(mcv),
            mcv_min = mmin(mcv),
            mcv_maxdiff = mmax(smcv),
            mcv_mindiff = mmin(smcv),
            mcv_tv = mmean(abs(smcv)))
print('summarize variable mcv')

# mono
table_mono <- table_group %>%
  summarise(mono_mean = mmean(mono),
            mono_max = mmax(mono),
            mono_min = mmin(mono),
            mono_maxdiff = mmax(smono),
            mono_mindiff = mmin(smono),
            mono_tv = mmean(abs(smono)))
print('summarize variable mono')

# mpv
table_mpv <- table_group %>%
  summarise(mpv_mean = mmean(mpv),
            mpv_max = mmax(mpv),
            mpv_min = mmin(mpv),
            mpv_maxdiff = mmax(smpv),
            mpv_mindiff = mmin(smpv),
            mpv_tv = mmean(abs(smpv)))
print('summarize variable mpv')

# neut
table_neut <- table_group %>%
  summarise(neut_mean = mmean(neut),
            neut_max = mmax(neut),
            neut_min = mmin(neut),
            neut_maxdiff = mmax(sneut),
            neut_mindiff = mmin(sneut),
            neut_tv = mmean(abs(sneut)))
print('summarize variable neut')

# platelet
table_platelet <- table_group %>%
  summarise(platelet_mean = mmean(platelet),
            platelet_max = mmax(platelet),
            platelet_min = mmin(platelet),
            platelet_maxdiff = mmax(splatelet),
            platelet_mindiff = mmin(splatelet),
            platelet_tv = mmean(abs(splatelet)))
print('summarize variable platelet')

# rbc
table_rbc <- table_group %>%
  summarise(rbc_mean = mmean(rbc),
            rbc_max = mmax(rbc),
            rbc_min = mmin(rbc),
            rbc_maxdiff = mmax(srbc),
            rbc_mindiff = mmin(srbc),
            rbc_tv = mmean(abs(srbc)))
print('summarize variable rbc')

# rdw
table_rdw <- table_group %>%
  summarise(rdw_mean = mmean(rdw),
            rdw_max = mmax(rdw),
            rdw_min = mmin(rdw),
            rdw_maxdiff = mmax(srdw),
            rdw_mindiff = mmin(srdw),
            rdw_tv = mmean(abs(srdw)))
print('summarize variable rdw')

# wbc
table_wbc <- table_group %>%
  summarise(wbc_mean = mmean(wbc),
            wbc_max = mmax(wbc),
            wbc_min = mmin(wbc),
            wbc_maxdiff = mmax(swbc),
            wbc_mindiff = mmin(swbc),
            wbc_tv = mmean(abs(swbc)))
print('summarize variable wbc')

# bind columns (remove ID columns)
table_summary <- bind_cols(table_baso,
                           table_eos[,-1],
                           table_hct[,-1],
                           table_hgb[,-1],
                           table_lymph[,-1],
                           table_mch[,-1],
                           table_mchc[,-1],
                           table_mcv[,-1],
                           table_mono[,-1],
                           table_mpv[,-1],
                           table_neut[,-1],
                           table_platelet[,-1],
                           table_rbc[,-1],
                           table_rdw[,-1],
                           table_wbc[,-1])
print('bind tables')

# save result
saveRDS(table_summary,file='R_data/labs_cbc_summary.rds')
print('save result')



