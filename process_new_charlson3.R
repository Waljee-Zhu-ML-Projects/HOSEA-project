# new indicators for identifying Charlson scores
# processing in blocks due to memory limits:
# block 1 1:4
# block 2 5:8
# block 3 9:12
# block 4 13:16

# helpers and code checking dictionary
source('R_code/hosea-project/utils_missingcharl.R')

# load data
alldxs <- list()
for(cc in as.character(c(1:5,'cx'))){
  alldxs[[cc]] <- readRDS(paste0('R_data/alldxs',cc,'_filter.rds'))
  print(paste0('Loaded chunk ',cc))
}

# loop over disease names
for(charl_name in charl_names[10:12]){
  print(paste0('For disease ',charl_name,':'))
  timestamp()
  
  # define functions to identify ICD9 and 10 codes
  check_code_icd9 <- charl_checks[[charl_name]][['icd9']]
  check_code_icd10 <- charl_checks[[charl_name]][['icd10']]
  
  charl_chunk <- list()

  for(cc in as.character(c(1:5,'cx'))){
    print(paste0('For chunk ',cc))
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
  }
  
  # combine individual chunks
  charl_chunk_rbound <- map_dfr(as.character(c(1:5,'cx')),function(cc){charl_chunk[[cc]]})
  print('Bind individual chunks')
  
  charl <- summarize(group_by(charl_chunk_rbound,ID),charl = max(charl))
  colnames(charl)[2] <- charl_name
  print('Combine over chunks')
  
  # save as an rds file
  saveRDS(charl,file=paste0('R_data/charlson_',charl_name,'.rds'))
  print('Saved')
}
