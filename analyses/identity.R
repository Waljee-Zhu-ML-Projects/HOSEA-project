# ==============================================================================
# PERFORMANCE STRATIFIED BY IDENTITY GROUPS
# Author: Simon Fontaine (simfont@umich.edu)
# ------------------------------------------------------------------------------




# ==============================================================================
# REQUIRED PACKAGES
library(dplyr)
library(xgboost)
library(magrittr)
library(ggplot2)
library(HOSEA)
theme_set(theme_minimal())
# ------------------------------------------------------------------------------




# ==============================================================================
# PATHS
setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = "./R_code/hosea-project/figures/mice/identity/"
imputed_data = "5-1test_mice_any.rds"
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------




# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------



# ==============================================================================
# PARAMETERS
outcome = "ANY"
seed = 0
age_bins = c(0, 45, 55, 65, 75, 100)
# ------------------------------------------------------------------------------



# ==============================================================================
# GET SCORES
proba = predict.HOSEA(imputed_df, imputer=NULL)
# ------------------------------------------------------------------------------


for(outcome in c("ANY", "EAC", "EGJAC")){


# ==============================================================================
# PREPARE DF
df =  raw_df$df %>% 
    select(id, casecontrol, age, gender, black, asian, indianalaskan, hawaiianpacific)
df %<>% mutate(
  white=ifelse(pmax(black, asian, indianalaskan, hawaiianpacific)==0, 1, 0),
  nonwhite=ifelse(pmax(black, asian, indianalaskan, hawaiianpacific)>0, 1, 0)
)
df %<>% mutate(
  age_bin=cut(age, breaks=age_bins)
)
df %<>% patch_outcome(raw_df$master, outcome=outcome)
df = proba %>% left_join(df, by="id")
# ------------------------------------------------------------------------------



# ==============================================================================
# AGE
ages = df %>% pull(age_bin) %>% levels()
roc_age = lapply(ages, function(age){
  wdf = df %>% filter(age_bin==!!age) %>% select(id, !!outcome, casecontrol) %>%
    rename(!!age:=!!outcome, y=casecontrol)
  out = roc(wdf)[[age]]
  out$curve$Age = paste0(age, " (", out$display.ci, ")")
  return(out$curve)
})
names(roc_age) = ages
roc_age %<>% bind_rows()


filepath = paste0(dir_figures, outcome, "_age.pdf")
title = paste0("ROC Curve per age strutum\n",
               "Cancer type: ", outcome)
g = ggplot(data=roc_age, aes(x=fpr, y=recall, color=Age)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  ggtitle(title)
g
ggsave(filepath, g, width=8, height=6)
# ------------------------------------------------------------------------------



# ==============================================================================
# AGE
df %<>% mutate(Sex=as.factor(ifelse(gender==1, "Male", "Female")))
sexes = df %>% pull(Sex) %>% levels()
roc_sex = lapply(sexes, function(sex){
  wdf = df %>% filter(Sex==!!sex) %>% select(id, !!outcome, casecontrol) %>%
    rename(!!sex:=!!outcome, y=casecontrol)
  out = roc(wdf)[[sex]]
  out$curve$Sex = paste0(sex, " (", out$display.ci, ")")
  return(out$curve)
})
names(roc_sex) = sexes
roc_sex %<>% bind_rows()


filepath = paste0(dir_figures, outcome, "_sex.pdf")
title = paste0("ROC Curve per sex\n",
               "Cancer type: ", outcome)
g = ggplot(data=roc_sex, aes(x=fpr, y=recall, color=Sex)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  ggtitle(title)
g
ggsave(filepath, g, width=8, height=6)
# ------------------------------------------------------------------------------



# ==============================================================================
# RACE
races = c("white", "nonwhite", "black", "asian", "indianalaskan", "hawaiianpacific")
roc_race = lapply(races, function(race){
  wdf = df %>% filter(!!sym(race)==1) %>% select(id, !!outcome, casecontrol) %>%
    rename(!!race:=!!outcome, y=casecontrol)
  out = roc(wdf)[[race]]
  out$curve$Race = paste0(race, " (", out$display.ci, ")")
  return(out$curve)
})
names(roc_race) = races
roc_race %<>% bind_rows()


filepath = paste0(dir_figures, outcome, "_race.pdf")
title = paste0("ROC Curve per race\n",
               "Cancer type: ", outcome)
g = ggplot(data=roc_race, aes(x=fpr, y=recall, color=Race)) + geom_line() +
  theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  ggtitle(title)
g
ggsave(filepath, g, width=8, height=6)
# ------------------------------------------------------------------------------

}
      
      
      
      
