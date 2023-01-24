# ==============================================================================
# MISC CALCULATIONS FOR THE ARTICLE
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
setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/comparison/")
dir_rocs = paste0("./R_code/hosea-project/figures/", imputation, "/roc/")
imputed_data = paste0("5-1test_", imputation, "_any.rds")
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------






# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------







# ==============================================================================
# RESULTS: MISSING PROPORTIONS

# overall 
raw_df$df %<>% mutate(
  complete = !(
    is.na(age) | 
    is.na(gender) | 
    is.na(smoke_current) | 
    is.na(bmi) | 
    is.na(gerd)
  )
)
with(raw_df$df, table(casecontrol, complete))

# ------------------------------------------------------------------------------











# ==============================================================================
# RESULTS: THRESHOLD & CLASSIFICATION METRICS

seed = 0 
representative=T

# representative
if(representative){
  set.seed(seed)
  imputed_df %<>% representative_sample()
}


models = load_models(
  files_meta=list(
    ANY=paste0("xgb_", imputation, "_any.meta")
  ),
  files_models=list(
    ANY=paste0("xgb_", imputation, "_any.model")
  )
)
proba = predict.HOSEA(imputed_df, imputer=NULL, models=models) %>% 
  select(id, !!outcome) %>% rename(HOSEA=!!outcome)
y = imputed_df %>% pull(casecontrol)
pred = proba$HOSEA

threshold = c(
  seq(0, 300, 5),
  seq(300, 500, 20), 
  seq(500, 1000, 50)
) / 100000
threshold = unique(threshold)
threshold = sort(threshold)

predmat = outer(pred, threshold, function(x, t) x>t)
ymat = outer(y, threshold, function(x, t) x)
out = data.frame(
  threshold=threshold*100000,
  N=colSums(!is.na(ymat)),
  p=colSums(predmat),
  tp=colSums(predmat*ymat),
  tn=colSums((1-predmat)*(1-ymat)),
  fp=colSums(predmat*(1-ymat)),
  fn=colSums(ymat*(1-predmat))
) %>% mutate(
  ppv=.data$tp/.data$p,
  tpr=.data$tp/(.data$tp+.data$fn),
  det_prev=.data$p/.data$N,
  spec=.data$tn/(.data$tn+.data$fp)
)
out %<>% mutate(tr_by_14=threshold/14)
out %<>% mutate(ppv=100*ppv, tpr=100*tpr, spec=100*spec, det_prev=100*det_prev)

write.csv(out, "./R_code/hosea-project/tables/srs/calibration/article.csv")

g_sensitivity = ggplot(
  data=out, 
  mapping=aes(x=threshold, y=tpr)
) + geom_line() + scale_x_continuous(trans="log10", breaks=c(5, 10, 20, 50, 100, 200, 500, 1000)) + 
  xlab("") + ylab("Sensitivity (%)") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept=235, linetype="dotted")+
  geom_hline(yintercept=out%>%filter(threshold==235)%>%pull(tpr), linetype="dotted")

g_specificity = ggplot(
  data=out, 
  mapping=aes(x=threshold, y=spec)
) + geom_line() + scale_x_continuous(trans="log10", breaks=c(5, 10, 20, 50, 100, 200, 500, 1000)) + 
  xlab("") + ylab("Specificity (%)") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_vline(xintercept=235, linetype="dotted")+
  geom_hline(yintercept=out%>%filter(threshold==235)%>%pull(spec), linetype="dotted")

g_ppv = ggplot(
  data=out, 
  mapping=aes(x=threshold, y=ppv)
) + geom_line() + scale_x_continuous(trans="log10", breaks=c(5, 10, 20, 50, 100, 200, 500, 1000)) + 
  xlab("") + ylab("PPV (%)") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(0, 1)+
  geom_vline(xintercept=235, linetype="dotted")+
  geom_hline(yintercept=out%>%filter(threshold==235)%>%pull(ppv), linetype="dotted")

g_det_prev = ggplot(
  data=out, 
  mapping=aes(x=threshold, y=det_prev)
) + geom_line() + scale_x_continuous(trans="log10", breaks=c(5, 10, 20, 50, 100, 200, 500, 1000)) + 
  xlab("Threshold (/100,000)") + ylab("Det. Prevalence") +
  geom_vline(xintercept=235, linetype="dotted")+
  geom_hline(yintercept=out%>%filter(threshold==235)%>%pull(det_prev), linetype="dotted")

g = cowplot::plot_grid(
  g_sensitivity,
  g_specificity,
  g_ppv,
  g_det_prev,
  ncol=1,
  align="v",
  axis="lr"
)

ggsave("./R_code/hosea-project/figures/srs/calibration/article.png", g, width=6, height=8, bg="white")

# ------------------------------------------------------------------------------







# ==============================================================================
# RESULTS: COHORT DATA
raw_df$master %<>% patch_staging()
raw_df$df %<>% left_join(raw_df$master %>% select(id, cancertype, nccn_stage_2017), by="id")

raw_df$df%<>%mutate(cancertype=ifelse(cancertype=="", "Control", cancertype))

out = raw_df$df %>% group_by(cancertype) %>% summarise(
  N=n(),
  
  mean_age=mean(age, na.rm=T),
  std_age=sd(age, na.rm=T),
  
  male_n=sum(gender, na.rm=T),
  male_prop=100*sum(gender, na.rm=T)/n(),
  
  bmi0.20_n=sum(bmi<20, na.rm=T),
  bmi0.20_prop=100*sum(bmi<20, na.rm=T)/n(),
  
  bmi20.25_n=sum(between(bmi, 20, 24.999), na.rm=T),
  bmi20.25_prop=100*sum(between(bmi, 20, 24.999), na.rm=T)/n(),
  
  bmi25.30_n=sum(between(bmi, 25, 29.999), na.rm=T),
  bmi25.30_prop=100*sum(between(bmi, 25, 29.999), na.rm=T)/n(),
  
  bmi30.35_n=sum(between(bmi, 30, 34.999), na.rm=T),
  bmi30.35_prop=100*sum(between(bmi, 30, 34.999), na.rm=T)/n(),
  
  bmi35.99_n=sum(bmi>=35, na.rm=T),
  bmi35.99_prop=100*sum(bmi>=35, na.rm=T)/n(),
  
  bmi_missing_n=sum(is.na(bmi)),
  bmi_missing_prop=100*sum(is.na(bmi))/n(),
  
  smoking_never_n=sum(pmax(smoke_current, smoke_former)==0, na.rm=T),
  smoking_never_prop=100*sum(pmax(smoke_current, smoke_former)==0, na.rm=T)/n(),
  
  smoking_former_n=sum(smoke_former, na.rm=T),
  smoking_former_prop=100*sum(smoke_former, na.rm=T)/n(),
  
  smoking_current_n=sum(smoke_current, na.rm=T),
  smoking_current_prop=100*sum(smoke_current, na.rm=T)/n(),
  
  smoking_missing_n=sum(is.na(smoke_current), na.rm=T),
  smoking_missing_prop=100*sum(is.na(smoke_current), na.rm=T)/n(),
  
  white_n=sum(pmax(black, asian, hawaiianpacific, indianalaskan)==0, na.rm=T),
  white_prop=100*sum(pmax(black, asian, hawaiianpacific, indianalaskan)==0, na.rm=T)/n(),
  
  black_n=sum(black, na.rm=T),
  black_prop=100*sum(black, na.rm=T)/n(),
  
  asian_n=sum(asian, na.rm=T),
  asian_prop=100*sum(asian, na.rm=T)/n(),
  
  other_n=sum(pmax(hawaiianpacific, indianalaskan), na.rm=T),
  other_prop=100*sum(pmax(hawaiianpacific, indianalaskan), na.rm=T)/n(),
  
  race_missing_n=sum(is.na(black), na.rm=T),
  race_missing_prop=100*sum(is.na(black), na.rm=T)/n(),
  
  gerd_n=sum(gerd, na.rm=T),
  gerd_prop=100*sum(gerd, na.rm=T)/n(),
  
  staging_0_n=sum(nccn_stage_2017=="0", na.rm=T),
  staging_0_prop=100*sum(nccn_stage_2017=="0", na.rm=T)/n(),
  
  staging_1_n=sum(nccn_stage_2017=="1", na.rm=T),
  staging_1_prop=100*sum(nccn_stage_2017=="1", na.rm=T)/n(),
  
  staging_2_n=sum(nccn_stage_2017=="2", na.rm=T),
  staging_2_prop=100*sum(nccn_stage_2017=="2", na.rm=T)/n(),
  
  staging_3_n=sum(nccn_stage_2017=="3", na.rm=T),
  staging_3_prop=100*sum(nccn_stage_2017=="3", na.rm=T)/n(),
  
  staging_34_n=sum(nccn_stage_2017=="3 or 4", na.rm=T),
  staging_34_prop=100*sum(nccn_stage_2017=="3 or 4", na.rm=T)/n(),
  
  staging_4_n=sum(nccn_stage_2017=="4", na.rm=T),
  staging_4_prop=100*sum(nccn_stage_2017=="4", na.rm=T)/n(),
  
  staging_u_n=sum(nccn_stage_2017=="u", na.rm=T),
  staging_u_prop=100*sum(nccn_stage_2017=="u", na.rm=T)/n(),
  
  staging_missing_n=sum(is.na(nccn_stage_2017), na.rm=T),
  staging_missing_prop=100*sum(is.na(nccn_stage_2017), na.rm=T)/n()
)

fout = out %>% summarise(
  cancertype=cancertype,
  N=format(N, trim=T),
  Age=paste0(format(mean_age, trim=T, digits=2, nsmall=2), " (", format(std_age, trim=T, digits=2, nsmall=2), ")"),
  Male=paste0(format(male_n, big.mark=","), " (", format(male_prop, trim=T, digits=1, nsmall=1), ")"),
  BMI="",
  `BMI<20`=paste0(format(bmi0.20_n, big.mark=","), " (", format(bmi0.20_prop, trim=T, digits=1, nsmall=1), ")"),
  `20<=BMI<25`=paste0(format(bmi20.25_n, big.mark=","), " (", format(bmi20.25_prop, trim=T, digits=1, nsmall=1), ")"),
  `25<=BMI<30`=paste0(format(bmi25.30_n, big.mark=","), " (", format(bmi25.30_prop, trim=T, digits=1, nsmall=1), ")"),
  `30<=BMI<35`=paste0(format(bmi30.35_n, big.mark=","), " (", format(bmi30.35_prop, trim=T, digits=1, nsmall=1), ")"),
  `35<=BMI`=paste0(format(bmi35.99_n, big.mark=","), " (", format(bmi35.99_prop, trim=T, digits=1, nsmall=1), ")"),
  `BMI Missing`=paste0(format(bmi_missing_n, big.mark=","), " (", format(bmi_missing_prop, trim=T, digits=1, nsmall=1), ")"),
  Smoking="",
  `Smoking Never`=paste0(format(smoking_never_n, big.mark=","), " (", format(smoking_never_prop, trim=T, digits=1, nsmall=1), ")"),
  `Smoking Former`=paste0(format(smoking_former_n, big.mark=","), " (", format(smoking_former_prop, trim=T, digits=1, nsmall=1), ")"),
  `Smoking Current`=paste0(format(smoking_current_n, big.mark=","), " (", format(smoking_current_prop, trim=T, digits=1, nsmall=1), ")"),
  `Smoking Missing`=paste0(format(smoking_missing_n, big.mark=","), " (", format(smoking_missing_prop, trim=T, digits=1, nsmall=1), ")"),
  Race="",
  `White`=paste0(format(white_n, big.mark=","), " (", format(white_prop, trim=T, digits=1, nsmall=1), ")"),
  `Black`=paste0(format(black_n, big.mark=","), " (", format(black_prop, trim=T, digits=1, nsmall=1), ")"),
  `Asian`=paste0(format(asian_n, big.mark=","), " (", format(asian_prop, trim=T, digits=1, nsmall=1), ")"),
  `Other`=paste0(format(other_n, big.mark=","), " (", format(other_prop, trim=T, digits=1, nsmall=1), ")"),
  `Race Missing`=paste0(format(race_missing_n, big.mark=","), " (", format(race_missing_prop, trim=T, digits=1, nsmall=1), ")") ,
  GERD=paste0(format(gerd_n, big.mark=","), " (", format(gerd_prop, trim=T, digits=1, nsmall=1), ")"),
  Stage="",
  `0`=paste0(format(staging_0_n, big.mark=","), " (", format(staging_0_prop, trim=T, digits=1, nsmall=1), ")"),
  I=paste0(format(staging_1_n, big.mark=","), " (", format(staging_1_prop, trim=T, digits=1, nsmall=1), ")"),
  II=paste0(format(staging_2_n, big.mark=","), " (", format(staging_2_prop, trim=T, digits=1, nsmall=1), ")"),
  III=paste0(format(staging_3_n, big.mark=","), " (", format(staging_3_prop, trim=T, digits=1, nsmall=1), ")"),
  `III/IV`=paste0(format(staging_34_n, big.mark=","), " (", format(staging_34_prop, trim=T, digits=1, nsmall=1), ")"),
  IV=paste0(format(staging_4_n, big.mark=","), " (", format(staging_4_prop, trim=T, digits=1, nsmall=1), ")"),
  U=paste0(format(staging_u_n, big.mark=","), " (", format(staging_u_prop, trim=T, digits=1, nsmall=1), ")"),
  `Stage Missing`=paste0(format(staging_missing_n, big.mark=","), " (", format(staging_missing_prop, trim=T, digits=1, nsmall=1), ")")
  
) %>% t()
write.csv(fout[, c(2, 3, 1)],
          "./R_code/hosea-project/tables/srs/cohort.csv")

# ------------------------------------------------------------------------------





