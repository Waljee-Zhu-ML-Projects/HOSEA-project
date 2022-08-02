setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(magrittr)
library(tidyr)
library(xgboost)
library(ggplot2)
theme_set(theme_minimal())
library(HOSEA)


# ==============================================================================
# load data & preprocessing
old_lft = "./unzipped_data/cases/labs_lft.rds"
new_lft = "./unzipped_data/cases/labs_lftcx.rds"
old_cbc = "./unzipped_data/cases/labs_cbc.rds"
new_cbc = "./unzipped_data/cases/labs_cbccx.rds"
samples = "./unzipped_data/cases/sample.rds"

old_lft %<>% HOSEA:::load_sas("labs")
new_lft %<>% HOSEA:::load_sas("labs")
old_cbc %<>% HOSEA:::load_sas("labs")
new_cbc %<>% HOSEA:::load_sas("labs")
samples %<>% HOSEA:::load_sas("sample")

colnames(old_lft) %<>% tolower()
colnames(new_lft) %<>% tolower()
colnames(old_cbc) %<>% tolower()
colnames(new_cbc) %<>% tolower()
colnames(samples) %<>% tolower()

old_lft %<>% left_join(samples %>% select(id, cancertype), by="id")
new_lft %<>% left_join(samples %>% select(id, cancertype), by="id")
old_cbc %<>% left_join(samples %>% select(id, cancertype), by="id")
new_cbc %<>% left_join(samples %>% select(id, cancertype), by="id")

cbc = old_cbc %>% full_join(
  y=new_cbc,
  suffix=c("_old", "_new"),
  by=c("id", "labdate", "cancertype")
)
lft = old_lft %>% full_join(
  y=new_lft,
  suffix=c("_old", "_new"),
  by=c("id", "labdate", "cancertype")
)
labs = full_join(
  x=lft,
  y=cbc,
  by=c("id", "labdate", "cancertype")
)


# ==============================================================================
# compare values and missing
vars = c("alkphos", "alt", "ast", "totprot", "baso", "eos", "hct", 
         "hgb", "lymph", "mch", 
         "mchc", "mcv", "mono", "mpv", "neut", "platelet", "rbc", 
         "rdw", "wbc")

summaries = lapply(vars, function(var){
  df = labs %>% select(id, cancertype, labdate, 
                       paste(var, "old", sep="_"), 
                       paste(var, "new", sep="_"))
  df %<>% rename(old=paste(var, "old", sep="_"), 
                 new=paste(var, "new", sep="_"))
  df %<>% filter(!is.na(cancertype))
  df %<>% mutate(
    diff=new-old,
    is_diff=ifelse((new-old)==0, 0, 1),
    both_na=is.na(old)*is.na(new),
    become_missing=(!is.na(old))*is.na(new),
    become_observed=is.na(old)*(!is.na(new))
  )
  
  tab = bind_rows(
    df %>% group_by(cancertype) %>% summarise_all(mean, na.rm=T) %>% select(-id, -labdate),
    df %>% summarise_all(mean, na.rm=T) %>% select(-id, -labdate) %>% mutate(cancertype="All")
  ) %>% mutate(variable=!!var)
}) %>% bind_rows()

summaries %<>% mutate(
  is_diff=100*is_diff,
  both_na=100*both_na,
  become_observed=100*become_observed,
  become_missing=100*become_missing
)

print(xtable::xtable(summaries %>% select(variable, cancertype, everything())), include.rownames=F)


# ==============================================================================
# compare models
