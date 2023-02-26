setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(HOSEA)
library(dplyr)
library(magrittr)

dir_out = "./R_data/processed_records/" # few small fixes


y0s = c(5, 5, 5, 5, 4, 3, 2, 4)
y1s = c(1, 3, 2, 4, 1, 1, 1, 2)
i=1

for(i in c(1, 2, 6, 8)){
  y0 = y0s[i]; y1 = y1s[i]
  
  metadata = list(
    HOSEA.version=packageVersion("HOSEA"),
    # HOSEA.commit=system("cd ./R_code/hosea-project/hosea-package/; git rev-parse HEAD", intern=T),
    years=c(-y0, -y1),
    datetime=NULL,
    dir_cases="unzipped_data/cases/",
    dir_controls="unzipped_data/",
    files_cases=c("sample.sas7bdat", "alldxs.sas7bdat", 'labs_a1c.sas7bdat', 
                  'labs_bmp.sas7bdat', 'labs_cbccx.sas7bdat', 
                  'labs_crp.sas7bdat', 'labs_lftcx.sas7bdat', 'labs_lipid.sas7bdat',
                  'allmeds.sas7bdat'),
    files_controls=c("sample.sas7bdat", "alldxs1.sas7bdat", "alldxs2.sas7bdat", 
                     "alldxs3.sas7bdat",
                     "alldxs4.sas7bdat", "alldxs5.sas7bdat", "alldxscx.sas7bdat", 
                     'labs_a1c.sas7bdat', 
                     'labs_bmp.sas7bdat', 'labs_cbc.sas7bdat', 
                     'labs_crp.sas7bdat', 'labs_lft.sas7bdat', 'labs_lipid.sas7bdat',
                     'allmeds.sas7bdat')
  )
  
  
  
  # new cases
  df = load_process_data(
    dir="unzipped_data/cases/",
    files_sample=c("sample.sas7bdat"),
    files_charlson=c("alldxs.sas7bdat"),
    files_labs=c('labs_a1c.sas7bdat', 'labs_bmp.sas7bdat', 'labs_cbccx.sas7bdat',
                 'labs_crp.sas7bdat', 'labs_lftcx.sas7bdat', 'labs_lipid.sas7bdat'),
    # files_labs=c('labs_a1c.sas7bdat', 'labs_bmp.sas7bdat', 'labs_cbc.sas7bdat', 
    #              'labs_crp.sas7bdat', 'labs_lft.sas7bdat', 'labs_lipid.sas7bdat'),
    files_meds=c('allmeds.sas7bdat'),
    start=-y0, end=-y1, 
    verbose=3
  )
  saveRDS(df, paste0(dir_out, y0, "-", y1, "_newcases.rds"))
  # controls and old cases
  df = load_process_data(
    dir="unzipped_data/",
    files_sample=c("sample.sas7bdat"),
    files_charlson=c("alldxs1.sas7bdat", "alldxs2.sas7bdat", "alldxs3.sas7bdat",
                     "alldxs4.sas7bdat", "alldxs5.sas7bdat", "alldxscx.sas7bdat"),
    files_labs=c('labs_a1c.sas7bdat', 'labs_bmp.sas7bdat', 'labs_cbc.sas7bdat',
                 'labs_crp.sas7bdat', 'labs_lft.sas7bdat', 'labs_lipid.sas7bdat'),
    files_meds=c('allmeds.sas7bdat'),
    start=-y0, end=-y1,
    verbose=3
  )
  saveRDS(df, paste0(dir_out, y0, "-", y1, "_original.rds"))
  
  # patch new cases
  df_old = readRDS(paste0(dir_out, y0, "-", y1, "_original.rds"))
  df_cases = readRDS(paste0(dir_out, y0, "-", y1, "_newcases.rds"))
  df_old$df %<>% filter(casecontrol==0)
  df_old$master %<>% filter(casecontrol==0)
  df_old$df %<>% bind_rows(df_cases$df)
  df_old$master %<>% bind_rows(df_cases$master)
  
  metadata$datetime = timestamp()
  df_old$metadata = metadata
  
  saveRDS(df_old, paste0(dir_out, y0, "-", y1, "_merged.rds"))
}