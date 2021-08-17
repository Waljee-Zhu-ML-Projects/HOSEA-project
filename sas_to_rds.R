# code to import files and save as rds objects

# sample
fname <- 'sample'
temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
temp <- arrange(temp,ID,NA)
saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
sample <- temp
rm(temp)
print(paste0('Saved ',fname,'.rds'))
timestamp()

# # colonoscopy
# fname <- 'colonoscopy'
# temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
# temp <- arrange(temp,ID,Procdate)
# saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
# rm(temp)
# print(paste0('Saved ',fname,'.rds'))
# timestamp()

# labs_a1c
fname <- 'labs_a1c'
temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
temp <- arrange(temp,ID,labdate)
saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
rm(temp)
print(paste0('Saved ',fname,'.rds'))
timestamp()

# labs_bmp
fname <- 'labs_bmp'
temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
temp <- arrange(temp,ID,labdate)
saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
rm(temp)
print(paste0('Saved ',fname,'.rds'))
timestamp()

# labs_cbc
fname <- 'labs_cbc'
temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
temp <- arrange(temp,ID,labdate)
saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
rm(temp)
print(paste0('Saved ',fname,'.rds'))
timestamp()

# labs_crp
fname <- 'labs_crp'
temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
temp <- arrange(temp,ID,labdate)
saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
rm(temp)
print(paste0('Saved ',fname,'.rds'))
timestamp()

# # labs_fobt
# fname <- 'labs_fobt'
# temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
# temp <- arrange(temp,ID,labdate)
# saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
# rm(temp)
# print(paste0('Saved ',fname,'.rds'))
# timestamp()

# labs_lft
fname <- 'labs_lft'
temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
temp <- arrange(temp,ID,labdate)
saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
rm(temp)
print(paste0('Saved ',fname,'.rds'))
timestamp()

# labs_lipid
fname <- 'labs_lipid'
temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
temp <- arrange(temp,ID,labdate)
saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
rm(temp)
print(paste0('Saved ',fname,'.rds'))
timestamp()

# # allmeds
# fname <- 'allmeds'
# temp <- read_sas(paste0('unzipped_data/',fname,'.sas7bdat'))
# #temp <- readRDS(paste0('R_data/',fname,'.rds'))
# temp <- arrange(temp,ID,Filldate)
# saveRDS(temp,file=paste0('R_data/',fname,'.rds'))
# rm(temp)
# print(paste0('Saved ',fname,'.rds'))
# timestamp()


