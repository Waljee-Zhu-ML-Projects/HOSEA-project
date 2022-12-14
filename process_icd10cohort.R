setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(HOSEA)
library(dplyr)
library(magrittr)

y0=3; y1=1

dir_out = "./R_data/processed_records/icd10cohort/"
dir_in = "./unzipped_data/icd10cohort/"
files_sample=c("sampleext.sas7bdat")
files_charlson=c(
  "alldxs1ext.sas7bdat",
  "alldxs2ext.sas7bdat",
  "alldxs3ext.sas7bdat",
  "alldxs4ext.sas7bdat",
  "alldxs5ext.sas7bdat"
  )
files_labs=c(
  'labs_a1cext.sas7bdat', 
  'labs_bmpext.sas7bdat', 
  'labs_cbcext.sas7bdat',
  'labs_crpext.sas7bdat', 
  'labs_lftext.sas7bdat', 
  'labs_lipidext.sas7bdat'
  )
files_meds=c('allmedsext.sas7bdat')
files = c(files_sample, files_charlson, files_labs, files_meds)

metadata = list(
  HOSEA.version=packageVersion("HOSEA"),
  years=c(-y0, -y1),
  datetime=NULL,
  dir=dir_in,
  files=files
)

df = load_process_data(
  dir=dir_in,
  files_sample=files_sample,
  files_charlson=files_charlson,
  files_labs=files_labs,
  files_meds=files_meds,
  start=-y0, end=-y1, 
  verbose=3
)


metadata$datetime = timestamp()
df$metadata = metadata
saveRDS(df, paste0(dir_out, y0, "-", y1, ".rds"))