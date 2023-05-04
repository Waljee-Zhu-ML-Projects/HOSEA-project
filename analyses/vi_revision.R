# ==============================================================================
# VARIABLE IMPORTANCE, PDPs, SHAP
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
imputation = "srs"
center = T
setwd('/nfs/turbo/umms-awaljee-secure/umms-awaljee-HOSEA/Peter files')
dir_tables = paste0("./R_code/hosea-project/tables/", imputation, "/")
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/vi_revisions/")
# ------------------------------------------------------------------------------




# ==============================================================================
# LOAD VIs
gain = list(
  all=read.csv(paste0(dir_tables, "variable_importance/vi_group.csv")),
  above50=read.csv(paste0(dir_tables, "variable_importance50p/vi_group.csv")),
  below50=read.csv(paste0(dir_tables, "variable_importance50m/vi_group.csv"))
)
gain %<>% bind_rows(.id="cohort")
shap = list(
  all=read.csv(paste0(dir_tables, "variable_importance/shap_group", ifelse(center, "_centered", ""), ".csv")),
  above50=read.csv(paste0(dir_tables, "variable_importance50p/shap_group", ifelse(center, "_centered", ""), ".csv")),
  below50=read.csv(paste0(dir_tables, "variable_importance50m/shap_group", ifelse(center, "_centered", ""), ".csv"))
)
shap %<>% bind_rows(.id="cohort")
# ------------------------------------------------------------------------------




# ==============================================================================
# RENAME FEATURES
groups = data.frame(
  group=c(
    "gender",      "bmi_weight",  "race",        "agentorange", "age",         "smoke",      
    "gerd",        "chf",         "ctd",         "dem",         "diab_c",      "hiv",        
    "mld",         "msld",        "para",        "rd",          "cd",          "copd",       
    "diab_nc",     "mi",          "pud",         "pvd",         "h2r",         "ppi",        
    "a1c",         "bun",         "calc",        "chlor",       "co2",         "creat",      
    "gluc",        "k",           "na",          "baso",        "eos",         "hct",        
    "lymph",       "mch",         "mchc",        "mcv",         "mono",        "mpv",        
    "neut",        "platelet",    "wbc",         "crp",         "alkphos",     "alt",        
    "ast",         "totprot",     "hdl",         "ldl",         "trig"
  ),
  display=c(
    "Sex",      "BMI/Weight",  "Race",        
    "Agent orange", "Age",         "Smoking status",      
    "GERD",        "CHF",         "CTD",         
    "Dementia",         "DM with compl.",      "HIV",        
    "MLD",         "MLDS",        "Paraplesia/Hemiplegia",        
    "Renal Disease",          "Cerebrovascular disease",          "COPD",       
    "DM no compl.",     "Myocardial infarction",          "PUD",         
    "PVD",         "H2R",         "PPI",        
    "HgbA1c",         "BUN",         "Calcium",        
    "Chlorine",       "CO2",         "Creatinine",      
    "Glucose",        "Potassium",           "Sodium",          
    "Basophil count",        "Eosinophil count",         "Hematocrit",        
    "Lymphocyte count",       "MCH",         "MCHC",        
    "MCV",         "Monocyte count",        "MPV",        
    "Neutrophil count",        "Platelet",    "WBC",         
    "CRP",         "AlkPhos",     "Alanine transferase",        
    "Aspartate transferase",         "Protein",     "HDL",         
    "LDL",         "Triglycerides"
  )
)
gain %<>% left_join(groups, by="group")
shap %<>% left_join(groups, by="group")
cohorts = c(
  all="All",
  above50=">=50 yo",
  below50="<50 yo"
)
gain$cohort_display = cohorts[gain$cohort]
shap$cohort_display = cohorts[shap$cohort]
# ------------------------------------------------------------------------------




# ==============================================================================
# PLOT

# each cohort
for(cname in names(cohorts)){
  g = ggplot(
    data=gain %>% filter(cohort==!!cname),
    mapping=aes(
      y=reorder(display, ANY),
      x=ANY * 100
    )) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle(paste0("Gain Variable importance by group (", cohorts[cname], ")")) +
    ylab("Feature group") + xlab("Gain (%)")
  filename = paste0(dir_figures, "gain_", cname, ".pdf")
  ggsave(filename, g, width=6, height=8, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=8, bg="white")
  
  g = ggplot(
    data=shap %>% filter(cohort==!!cname),
    mapping=aes(
      y=reorder(display, ANY),
      x=ANY
    )) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle(paste0("SHAP Variable importance by group (", cohorts[cname], ")")) +
    ylab("Feature group") + xlab(ifelse(center, "mean|SHAP - median(SHAP)|", "mean|SHAP|"))
  filename = paste0(dir_figures, "shap_", cname, "", ifelse(center, "_centered", ""), ".pdf")
  ggsave(filename, g, width=6, height=8, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=8, bg="white")
}

# combined
g = ggplot(
  data=gain,
  mapping=aes(
    y=reorder(display, ANY),
    x= ANY * 100,
    fill=cohort_display
  )) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle(paste0("Gain Variable importance by group")) +
  ylab("Feature group") + xlab("Gain (%)")+ 
  theme(legend.position=c(0.8, 0.2)) + 
  labs(fill="Cohort")
filename = paste0(dir_figures, "gain.pdf")
ggsave(filename, g, width=6, height=8, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=8, bg="white")

g = ggplot(
  data=shap,
  mapping=aes(
    y=reorder(display, ANY),
    x=ANY,
    fill=cohort_display
  )) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle(paste0("SHAP Variable importance by group")) +
  ylab("Feature group") + xlab(ifelse(center, "mean|SHAP - median(SHAP)|", "mean|SHAP|")) + 
  theme(legend.position=c(0.8, 0.2)) + 
  labs(fill="Cohort")
filename = paste0(dir_figures, "shap", ifelse(center, "_centered", ""), ".pdf")
ggsave(filename, g, width=6, height=8, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=8, bg="white")
# ------------------------------------------------------------------------------


