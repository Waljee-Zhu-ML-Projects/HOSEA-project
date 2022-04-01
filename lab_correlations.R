setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(HOSEA)
library(ggcorrplot)

# ==============================================================================
# paths and parameters
dir_path = "R_data/processed_records/"
dir = "./unzipped_data/"
dir_figures = "R_code/hosea-project/figures/"
model_path = "R_data/results/models/XGB_nALL_typeANY.rds"

# =========================================================
# read in model
results = readRDS(model_path)
xgb_fit = results$xgb_fit
quantiles = results$quantiles
test_ids = results$test_ids
rm(results); gc()

# =========================================================
# read in data
file_path = paste0(dir_path, "5-1.rds")
df = readRDS(file_path)
master = df$master
# subset to test set
master %<>% filter(ID %in% test_ids)

df = df$df
df %<>% filter(ID %in% test_ids)

lab_vars = HOSEA:::lab_vars()

# =========================================================
# get lab data
which=HOSEA:::lab_types()

dfs = list()

for(file in which){
  cat(paste0("- ", file, " ...\n"))
  src_df = load_sas(paste0(dir, file, ".sas7bdat"), "file")
  src_df %<>% filter(ID %in% master$ID)
  dfs[[file]] = src_df
  # subtypes = tail(colnames(src_df), -2)
  # # restrict to prediction window
  # src_df %<>% left_join(master, by="ID")
  # src_df %<>% filter((labdate>=start)&(labdate<=end))
  # # ensure ordered
  # src_df %<>% arrange(ID, labdate)
  # 
  # cat("  ")
  # for(type in subtypes){
  #   cat(paste0(type, " "))
  #   tmp = src_df %>% select(ID, labdate, type)
  #   colnames(tmp) = c("ID", "labdate", "var")
  #   dfs[[paste(file, type)]] = tmp
  # }
  cat("\n  ...done.\n")
}
dff = dfs %>% purrr::reduce(full_join, by=c("ID", "labdate"))

dff2 = dff %>% select(-c(ID, labdate))

cormat = dff2 %>% cor(use="pairwise.complete.obs")

g = ggcorrplot(cormat, method="square")
ggsave(paste0(dir_figures, "lab_corr.pdf"), width=10, height=10)

cormat_df = cormat%>% data.frame()
cormat_df$var1 = rownames(cormat_df)
cormat_long = cormat_df %>% pivot_longer(cols=-var1)
cormat_long %<>% filter(value!=1)

rank = cormat_long$value %>% abs() %>% order()

xtable::xtable(cormat_long[rank, ] %>% tail(20), digits=3)
