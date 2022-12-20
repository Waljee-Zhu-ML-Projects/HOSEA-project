# ==============================================================================
# LOGISTIC REGRESSION
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
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/labs/")
imputed_data = paste0("5-1test_", imputation, "_any.rds")
raw_data = "5-1_merged.rds"
# ------------------------------------------------------------------------------






# ==============================================================================
# READ IN DATA
imputed_df = readRDS(paste0(dir_imputed_data, imputed_data))
raw_df = readRDS(paste0(dir_raw_data, raw_data))
# ------------------------------------------------------------------------------







# ==============================================================================
# PREPARE DATA
# df = imputed_df
set.seed(0)
# subsample for SHAP/PDPs
# df = imputed_df %>% sample_n(1000000)
# subsample controls
df = imputed_df %>% filter(casecontrol==1) %>% 
  bind_rows(imputed_df%>%filter(casecontrol==0) %>% sample_n(100000))
# compute some variables
df %<>% mutate(
  white=(black+asian+hawaiianpacific+indianalaskan)==0,
  obesity=bmi>30,
  smoke_any=pmax(smoke_former, smoke_current),
  anion_gap_mean=na_mean + k_mean - co2_mean - chlor_mean
  )
df %<>% left_join(
  raw_df$df %>% select(id, ppi_max) %>%
    mutate(
      ppi=factor(
        ifelse(is.na(ppi_max), "NONE", ifelse(ppi_max<30, "LOW", "HIGH")), 
        levels=c("NONE", "LOW", "HIGH"))
      ),
  by="id"
)
# ------------------------------------------------------------------------------






# ==============================================================================
# MODELS
models = list(
  "age"=c(age=T),
  "sex"=c(gender=F),
  "gerd"=c(gerd=F),
  "white"=c(white=F),
  "obesity"=c(obesity=F),
  "smoking"=c(smoke_any=F),
  "ppi"=c(ppi=F),
  "known_predictors"=c(age=T, gender=F, gerd=F, white=F, obesity=F, smoke_any=F),
  "anion_gap"=c(anion_gap_mean=T),
  "known_predictors_anion_gap"=c(age=T, gender=F, gerd=F, white=F, obesity=F, 
                                 smoke_any=F, anion_gap_mean=T),
  "ppi_anion_gap"=c(ppi=F, anion_gap_mean=T),
  "known_predictors_ppi_anion_gap"=c(age=T, gender=F, gerd=F, white=F, obesity=F, 
                                     smoke_any=F, anion_gap_mean=T, ppi=F),
  "known_predictors_ppi"=c(age=T, gender=F, gerd=F, white=F, obesity=F, 
                                     smoke_any=F, ppi=F)
)
# ------------------------------------------------------------------------------






# ==============================================================================
# FORMULA BUILDER
to_formula = function(x){
  vnames = names(x)
  vterms = sapply(vnames, function(vname) ifelse(x[vname], paste0("s(", vname, ")"), vname))
  f = paste(vterms, collapse=" + ")
  f = paste0("casecontrol ~ ", f)
  return(f)
}
# ------------------------------------------------------------------------------






# ==============================================================================
# FIT ALL MODELS
fits = lapply(names(models), function(mname){
  model = models[[mname]]
  marginal = length(model)==1
  fit = mgcv::gam(formula=to_formula(model) %>% formula, family=binomial, data=df)
  fit$name = ifelse(marginal, "marginal", mname)
  return(fit)
})
names(fits) = names(models)
# ------------------------------------------------------------------------------





# ==============================================================================
# COEFFICIENT DF
coef_df = lapply(fits, function(fit){
  coefs = summary(fit)$p.table %>% as.data.frame()
  coefs %<>% mutate(
    L=Estimate-1.96*`Std. Error`,
    U=Estimate+1.96*`Std. Error`
  ) %>% mutate(
    Odds=exp(Estimate),
    LOdds=exp(L),
    UOdds=exp(U)
  ) %>% mutate(
    model=fit$name
  ) 
  coefs$feature = rownames(coefs)
  return(coefs %>% filter(feature!="(Intercept)"))
}) %>% bind_rows()
# ------------------------------------------------------------------------------






# ==============================================================================
# PLOT COEFFICIENTS
filepath = paste0(dir_figures, "logreg_coefs.pdf")

g = ggplot(data=coef_df) + 
  geom_point(aes(x=model, y=Odds, color=model)) +
  geom_errorbar(mapping=aes(x=model, ymin=LOdds, ymax=UOdds, color=model)) + 
  facet_wrap(~feature, scales="free") + 
  ggtitle("Logistic regression: estimate comparison") + 
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position="bottom"
    )
ggsave(filepath, g, width=8, height=8, bg="white")
ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=8, height=8, bg="white")
# ------------------------------------------------------------------------------











# ==============================================================================
# COEFFICIENT DF
spline_df = lapply(fits, function(fit){
  if(length(fit$smooth)>0){
    plots = plot(fit) 
    coefs = lapply(plots, function(plt){
      feature = plt$xlab
      x = plt$x
      estimate = plt$fit
      se = plt$se
      df = data.frame(
        feature=feature,
        x=x,
        estimate=estimate,
        L=estimate-1.96*se,
        U=estimate+1.96*se
      ) %>% mutate(
        Odds=exp(estimate),
        LOdds=exp(L),
        UOdds=exp(U)
      )
      return(df)
    }) %>% bind_rows() %>% mutate(model=fit$name)
    return(coefs)
  }
}) %>% bind_rows()
# ------------------------------------------------------------------------------








# ==============================================================================
# PLOT SPLINES
filepath = paste0(dir_figures, "logreg_splines.pdf")
xlims = list(age=c(18, 90), anion_gap_mean=c(5, 20))

spline_df = lapply(names(xlims), function(vname){
  df = spline_df %>% filter(feature==!!vname) %>%
    filter(x>=xlims[[vname]][1]) %>%
    filter(x<=xlims[[vname]][2])
  return(df)
}) %>% bind_rows()

g = ggplot(data=spline_df) + 
  geom_line(mapping=aes(x=x, y=Odds, color=model)) +
  # geom_ribbon(mapping=aes(x=x, ymin=LOdds, ymax=UOdds, color=model), alpha=0.2) + 
  facet_wrap(~feature, scales="free", ncol=1) + 
  ggtitle("Logistic regression: splines comparison") + 
  theme(legend.position="bottom")
ggsave(filepath, g, width=8, height=8, bg="white")
ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=8, height=8, bg="white")
# ------------------------------------------------------------------------------






# ==============================================================================
# PPI x ANION GAP
filepath = paste0(dir_figures, "logreg_splines_ppi_x_anion_gap.pdf")
fit = mgcv::gam(formula("casecontrol ~ s(age) + white + gerd + smoke_any + 
                        gender + obesity + 
                        s(anion_gap_mean, by=ppi)"), family=binomial, data=df)
plots = plot(fit) 
plots = list(
  NONE=plots[[2]],
  LOW=plots[[3]],
  HIGH=plots[[4]]
)
coefs = lapply(names(plots), function(ppilvl){
  plt = plots[[ppilvl]]
  feature = ppilvl
  x = plt$x
  estimate = plt$fit
  se = plt$se
  df = data.frame(
    feature=feature,
    x=x,
    estimate=estimate,
    L=estimate-1.96*se,
    U=estimate+1.96*se
  ) %>% mutate(
    Odds=exp(estimate),
    LOdds=exp(L),
    UOdds=exp(U)
  )
  return(df)
}) %>% bind_rows()




g = ggplot(data=coefs) + 
  geom_line(mapping=aes(x=x, y=Odds, color=feature)) +
  # geom_ribbon(mapping=aes(x=x, ymin=LOdds, ymax=UOdds, color=model), alpha=0.2) + 
  ggtitle("Logistic regression: splines by PPI level") + 
  theme(legend.position="bottom") + 
  xlim(0, 30) + ylim(0,1.5)
ggsave(filepath, g, width=8, height=6, bg="white")
ggsave(stringr::str_replace(filepath, "pdf", "png"), g, width=8, height=6, bg="white")
# ------------------------------------------------------------------------------






