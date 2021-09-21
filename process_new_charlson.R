# new indicators for identifying Charlson scores

# example: COPD

# define the name of the disease as it appears in the full data
charl_name <- 'copd'
# define functions to identify ICD9 and 10 codes
check_code_icd9 <- function(code){
  as.integer(substr(code,1,2)=='49')
}

check_code_icd10 <- function(code){
  as.integer(substr(code,1,3)=='J44')
}

charl_chunk <- list()
length(charl_chunk) <- 6

for(cc in c(1:5,'cx')){
  print(paste0('For chunk ',chunk))
  # load data
  alldxs <- readRDS(paste0('R_data/alldxs',chunk,'_filter.rds'))
  print('Loaded')
  # new variables
  alldxs <- mutate(alldxs,charl9=check_code_icd9(icd9code))
  if(cc == 'cx'){
    alldxs <- mutate(alldxs,charl10=check_code_icd10(icd10code))
  }
  else{
    alldxs <- mutate(alldxs,charl10=0)
  }
  print('Checked codes')
  # summarize
  charl_chunk[[cc]] <- summarize(group_by(alldxs,ID),charl=max(charl9 + charl10))
  print('Summarized chunk')
}

# combine individual chunks
charl_chunk_rbound <- map_dfr(c(1:5,'cx'),function(ii){charl_chunk[[ii]]})
print('Bind individual chunks')

charl <- summarize(group_by(charl_chunk_rbound,ID),charl = max(charl))
colnames(charl)[2] <- charl_name
print('Combine over chunks')

# save as an rds file
saveRDS(charl,file=paste0('R_data/charl_',charl_name,'.rds'))
print('Saved')

