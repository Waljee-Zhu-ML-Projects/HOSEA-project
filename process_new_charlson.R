# new indicators for identifying Charlson scores

# helpers and code checking dictionary
source('R_code/hosea-project/utils_missingcharl.R')

# load data
alldxs <- list()
length(alldxs) <- 6
for(cc in c(1:5,'cx')){
  alldxs[[cc]] <- readRDS(paste0('R_data/alldxs',chunk,'_filter.rds'))
}

# loop over disease names
for(charl_name in charl_names){
  print(paste0('For disease ',charl_name,':'))
  timestamp()
  
  # define functions to identify ICD9 and 10 codes
  check_code_icd9 <- charl_checks[[charl_name]][['icd9']]
  
  check_code_icd10 <- charl_checks[[charl_name]][['icd10']]
  
  charl_chunk <- list()
  length(charl_chunk) <- 6
  
  for(cc in c(1:5,'cx')){
    print(paste0('For chunk ',chunk))
    # new variables
    alldxs_chunk <- mutate(alldxs[[cc]],charl9=check_code_icd9(icd9code))
    if(cc == 'cx'){
      alldxs_chunk <- mutate(alldxs_chunk,charl10=check_code_icd10(icd10code))
    }
    else{
      alldxs_chunk <- mutate(alldxs_chunk,charl10=0)
    }
    print('Checked codes')
    # summarize
    charl_chunk[[cc]] <- summarize(group_by(alldxs_chunk,ID),charl=max(charl9 + charl10))
    print('Summarized chunk')
    
    # remove chunk
    rm(alldxs_chunk)
    gc()
  }
  
  # combine individual chunks
  charl_chunk_rbound <- map_dfr(c(1:5,'cx'),function(cc){charl_chunk[[cc]]})
  print('Bind individual chunks')
  
  charl <- summarize(group_by(charl_chunk_rbound,ID),charl = max(charl))
  colnames(charl)[2] <- charl_name
  print('Combine over chunks')
  
  # save as an rds file
  saveRDS(charl,file=paste0('R_data/charlson_',charl_name,'.rds'))
  print('Saved')
}
