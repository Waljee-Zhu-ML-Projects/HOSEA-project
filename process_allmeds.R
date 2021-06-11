# processing allMeds table

# helper functions (mmean, mmin, mmax)
source('R_code/hosea-project/utils.R')

# load master
#master <- readRDS('R_data/subsample/sub_master.rds')
master <- readRDS('R_data/master.rds')

# load table
#allmeds <- readRDS('R_data/subsample/sub_allmeds.rds')
allmeds <- readRDS('R_data/allmeds.rds')
print('load data')

# left join master to get windowing dates
allmeds <- left_join(allmeds,master,by='ID')
print('link prediction window')
# filter only dates in each window
allmeds_keep <- filter(allmeds,(newenddate >= start) & (Filldate <= end))
print('filter lab dates')

# filter two tables based on Med_type
allmeds_h2r_keep <- filter(allmeds_keep,Med_Type=='H2R')
allmeds_ppi_keep <- filter(allmeds_keep,Med_Type=='PPI')

# processing H2R

# threshold dates and get next fill
allmeds_h2r <- allmeds_h2r_keep %>%
  mutate(fill_thresh = pmax(start,Filldate),
         end_thresh = pmin(end,newenddate),
         next_ID = lead(ID),
         next_fill = lead(Filldate,default=Inf))
rm(allmeds_h2r_keep)

# Set skips to NA or Inf
skip_rows <- !(allmeds_h2r$ID == allmeds_h2r$next_ID)
allmeds_h2r$next_ID[skip_rows] <- NA
allmeds_h2r$next_fill[skip_rows] <- Inf

# create new tables
keep_enddate <- (allmeds_h2r$end_thresh < allmeds_h2r$next_fill)

# if end_thresh==end want to make the corresponding dose the same as
# the previous dose

allmeds_h2r_clean <- bind_rows(data.frame(ID=allmeds_h2r$ID,
                                          meddate=allmeds_h2r$fill_thresh,
                                          dose=allmeds_h2r$dd),
                               data.frame(ID=allmeds_h2r$ID[keep_enddate],
                                          meddate=allmeds_h2r$end_thresh[keep_enddate],
                                          dose=ifelse(allmeds_h2r$end_thresh==allmeds_h2r$end,allmeds_h2r$dd,0)[keep_enddate]))
rm(allmeds_h2r)

# arrange full table by ID then meddate
allmeds_h2r_clean <- arrange(allmeds_h2r_clean,ID,meddate)

# # lagged variables (ID, labdate and value fields)
allmeds_h2r_lead <- allmeds_h2r_clean %>%
  mutate(ID_next = lead(ID),
         meddate_next = lead(meddate),
         dose_next = lead(dose))
rm(allmeds_h2r_clean)
print('compute leading variables')

# find where IDs change
skip_rows <- !(allmeds_h2r_lead$ID == allmeds_h2r_lead$ID_next)
# make lagged values NAs
allmeds_h2r_lead$ID_next[skip_rows] <- NA
allmeds_h2r_lead$meddate_next[skip_rows] <- NA
allmeds_h2r_lead$dose_next[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
allmeds_h2r_all <- allmeds_h2r_lead %>%
  mutate(dmeddate = pmax(1,meddate_next - meddate),
         ddose = dose_next - dose) %>%
  mutate(sdose = ddose / dmeddate,
         pdose = dose*dmeddate)
rm(allmeds_h2r_lead)
print('calculate slopes')

# group by ID
allmeds_h2r_group <- allmeds_h2r_all %>%
  group_by(ID)
print('group by ID')

# summary stats for dosage variable
summary_h2r <- allmeds_h2r_group %>%
  summarise(h2r_int = ssum(pdose),
            h2r_mean = mmean(dose),
            h2r_max = mmax(dose),
            h2r_maxdiff = mmax(sdose),
            h2r_tv = mmean(abs(sdose)))
print('summarize h2r')

gc()

# processing PPI

# threshold dates and get next fill
allmeds_ppi <- allmeds_ppi_keep %>%
  mutate(fill_thresh = pmax(start,Filldate),
         end_thresh = pmin(end,newenddate),
         next_ID = lead(ID),
         next_fill = lead(Filldate,default=Inf))
rm(allmeds_ppi_keep)

# Set skips to NA or Inf
skip_rows <- !(allmeds_ppi$ID == allmeds_ppi$next_ID)
allmeds_ppi$next_ID[skip_rows] <- NA
allmeds_ppi$next_fill[skip_rows] <- Inf

# create new tables
keep_enddate <- allmeds_ppi$end_thresh < allmeds_ppi$next_fill

allmeds_ppi_clean <- bind_rows(data.frame(ID=allmeds_ppi$ID,
                                          meddate=allmeds_ppi$fill_thresh,
                                          dose=allmeds_ppi$dd),
                               data.frame(ID=allmeds_ppi$ID[keep_enddate],
                                          meddate=allmeds_ppi$end_thresh[keep_enddate],
                                          dose=ifelse(allmeds_ppi$end_thresh==allmeds_ppi$end,allmeds_ppi$dd,0)[keep_enddate]))
rm(allmeds_ppi)

# arrange full table by ID then meddate
allmeds_ppi_clean <- arrange(allmeds_ppi_clean,ID,meddate)

# # lagged variables (ID, labdate and value fields)
allmeds_ppi_lead <- allmeds_ppi_clean %>%
  mutate(ID_next = lead(ID),
         meddate_next = lead(meddate),
         dose_next = lead(dose))
rm(allmeds_ppi_clean)
print('compute leading variables')

# find where IDs change
skip_rows <- !(allmeds_ppi_lead$ID == allmeds_ppi_lead$ID_next)
# make lagged values NAs
allmeds_ppi_lead$ID_next[skip_rows] <- NA
allmeds_ppi_lead$meddate_next[skip_rows] <- NA
allmeds_ppi_lead$dose_next[skip_rows] <- NA
print('remove ID changes')

# get diffs and slopes
allmeds_ppi_all <- allmeds_ppi_lead %>%
  mutate(dmeddate = pmax(1,meddate_next - meddate),
         ddose = dose_next - dose) %>%
  mutate(sdose = ddose / dmeddate,
         pdose = dose*dmeddate)
rm(allmeds_ppi_lead)
print('calculate slopes')

# group by ID
allmeds_ppi_group <- allmeds_ppi_all %>%
  group_by(ID)
print('group by ID')

# summary stats for dosage variable
summary_ppi <- allmeds_ppi_group %>%
  summarise(ppi_int = ssum(pdose),
            ppi_mean = mmean(dose),
            ppi_max = mmax(dose),
            ppi_maxdiff = mmax(sdose),
            ppi_tv = mmean(abs(sdose)))
print('summarize ppi')

# recombine (outer join) and save object
summary_allmeds <- full_join(summary_h2r,
                             summary_ppi,
                             by='ID')
print('join tables')

# check results against old code
# summary_allmeds_check <- readRDS('R_data/subsample/sub_allmeds_summary.rds')

# save result
#saveRDS(summary_allmeds,file='R_data/subsample/sub_allmeds_summary.rds')
saveRDS(summary_allmeds,file='R_data/allmeds_summary.rds')
print('save result')

# check results against subsample
#summary_allmeds_check <- readRDS('R_data/subsample/sub_allmeds_summary.rds')

# spot check that it worked for ID 1463679

