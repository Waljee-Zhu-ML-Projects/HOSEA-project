setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(HOSEA)
library(dplyr)
library(magrittr)

dir_out = "./R_data/processed_records/" # few small fixes


y0s = c(5, 5, 5, 5, 4, 3, 2, 4)
y1s = c(1, 3, 2, 4, 1, 1, 1, 2)
i=1

for(i in c(2, 6, 8)){
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







d0 = master %>% select(id, casecontrol, end)
d1 = src_df %>% group_by(id) %>% summarize_all(min)




src_df %>% pull(id) %>% unique() %>% length()









# i = 4
# source("R_code/hosea-project/hosea-package/R/data_utils.R")
# source("R_code/hosea-project/hosea-package/R/read_sas.R")

for(i in seq(2, 8)){
  y0 = y0s[i]; y1 = y1s[i]
  df = load_process_data(dir=dir,start=-y0, end=-y1)
  saveRDS(df, paste0("R_data/processed_records/cases_", y0, "-", y1, ".rds"))
  # patch
  df_old = readRDS(paste0("R_data/processed_records/", y0, "-", y1, ".rds"))
  df_old$df %<>% filter(CaseControl==0)
  df_old$master %<>% filter(CaseControl==0)
  df_old$df %<>% bind_rows(df$df)
  df_old$master %<>% bind_rows(df$master)
  saveRDS(df_old, paste0("R_data/processed_records/", y0, "-", y1, ".rds"))
}

# ICD10-only 
df = load_process_data(dir=dir,start=-5, end=-1, icd10=T, icd10startdate=17229)
# patch
df_old = readRDS(paste0("R_data/processed_records/5-1_icd10only.rds"))
df_old$df %<>% filter(CaseControl==0)
df_old$master %<>% filter(CaseControl==0)
df_old$df %<>% bind_rows(df$df)
df_old$master %<>% bind_rows(df$master)
saveRDS(df_old, paste0("R_data/processed_records/5-1_icd10only.rds"))

# ICD9-only 
df = load_process_data(dir=dir,start=-5, end=-1, icd9=T, icd9enddate=17229-3*31)
# patch
df_old = readRDS(paste0("R_data/processed_records/5-1_icd9only.rds"))
df_old$df %<>% filter(CaseControl==0)
df_old$master %<>% filter(CaseControl==0)
df_old$df %<>% bind_rows(df$df)
df_old$master %<>% bind_rows(df$master)
saveRDS(df_old, paste0("R_data/processed_records/5-1_icd9only.rds"))

#
df_new = load_sas(paste0(dir, "sample.rds"), "sample")
df_old = load_sas(paste0(dir, "old/sample.rds"), "sample")

tab_new = table(
  df_new$CaseControl,
  df_new$BMI %>% cut(breaks=c(0, 20, 25, 30, 35, 40, 100)), 
  useNA="ifany"
)

tab_old = table(
  df_old$df$CaseControl,
  df_old$df$bmi %>% cut(breaks=c(0, 20, 25, 30, 35, 40, 100)), 
  useNA="ifany"
)

round(tab_new / sum(tab_new)*100, 3)
round(tab_old / sum(tab_old)*100, 3)

freq_new = round(tab_new / rbind(tab_new%>%colSums(), tab_new%>%colSums())*100, 2)
freq_old = round(tab_old / rbind(tab_old%>%colSums(), tab_old%>%colSums())*100, 2)

# subset to cases
cases_new = df_new %>% filter(CaseControl == 1) %>% select(IndexDate, facility, BMI)
cases_old = df_old %>% filter(CaseControl == 1) %>% select(IndexDate, facility, bmi)
cases = full_join(cases_new, cases_old, by="IndexDate")
cases %<>% filter(facility.x==facility.y)

plot(cases$bmi, cases$BMI, xlab="BMI(old)", ylab="BMI(new)")
abline(a=0, b=1)


# HOSEA Package test
library(HOSEA)
out = load_process_data(
  dir="unzipped_data/cases/",
  files_sample=c("sample.sas7bdat"),
  files_charlson=c("alldxs.sas7bdat"),
  files_labs=c('labs_a1c.sas7bdat', 'labs_bmp.sas7bdat', 'labs_cbc.sas7bdat', 
               'labs_crp.sas7bdat', 'labs_lft.sas7bdat', 'labs_lipid.sas7bdat'),
  files_meds=c("allmeds.sas7bdat"),
  start=-4, end=0, 
  verbose=3
)
pred = predict.HOSEA(out$df, n_imputations=10)
head(pred)
pred = predict.HOSEA(out$df, n_imputations=10, xgb_fits=list(ANY=XGB_ANY))
head(pred)
df=out$df
model=XGB_ANY
xgb_fit=model$xgb_fit
quantiles=model$quantiles
missing_prop=model$missing_prop

n_quantiles = nrow(quantiles)
qs = seq(n_quantiles)/(n_quantiles + 1)
df_quantiles = df %>%
  summarise_all(quantile, na.rm=T, probs=qs, type=1)
ql = which.max(qs > 0.005)
qu = which.min(qs < 0.995)

varname="totprot_mean"
l = (quantiles %>% pull(varname))[ql]
u = (quantiles %>% pull(varname))[qu]
plot(quantiles%>%pull(varname), df_quantiles%>%pull(varname))
abline(0,1)
abline(v=c(l, u), col="red")
