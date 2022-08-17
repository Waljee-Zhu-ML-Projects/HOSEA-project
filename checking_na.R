setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(magrittr)
library(ggplot2)
theme_set(theme_minimal())

# =========================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir_figures = "R_code/hosea-project/figures/"

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1_merged.rds")
df = readRDS(file_path)
master = df$master
df = df$df

# =========================================================
# new controls
master %<>% mutate(new=ifelse(end>17230, "new", "old"))
df %<>% left_join(master%>%select(id, end, new), by="id")
table(df$new, is.na(df$a1c_mean))

ggplot(
  data=df %>% 
    mutate(endgr=round(end/10)*10) %>% 
    group_by(endgr) %>% summarize(prop_na=mean(is.na(k_mean))),
  mapping=aes(x=endgr, y=prop_na)
) + geom_line()
