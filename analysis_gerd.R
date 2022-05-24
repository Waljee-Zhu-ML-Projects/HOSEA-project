setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(magrittr)

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"

# =========================================================
# read in data
df_new = readRDS(paste0(dir_path, "5-1_test_merged.rds"))$df
df_old = readRDS(paste0(dir_path, "5-1.rds"))$df

ids = intersect(df_new$id, df_old$ID)
df_both = data.frame(
  id = ids,
  gerd_new = df_new %>% filter(id %in%ids) %>% pull(gerd),
  gerd_old = df_old %>% filter(ID %in% ids) %>% pull(GerdAtIndex),
  casecontrol = df_new %>% filter(id %in%ids) %>% pull(casecontrol)
)
tab_both = with(df_both, table(gerd_new, gerd_old, casecontrol))
prop_both = 100 * tab_both[,,2] / (tab_both[,,1] + tab_both[,,2])
tab_old = tab_both[1,,] + tab_both[2,,]
tab_new = tab_both[,1,] + tab_both[,2,]

prop_new = apply(tab_new, 1, function(x) x/sum(x)) * 100
prop_old = apply(tab_old, 2, function(x) x/sum(x)) * 100
