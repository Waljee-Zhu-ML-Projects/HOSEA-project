# ==============================================================================
# INVESTIGATING LAB FEATURE ASSOCIATION
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
# PARAMETERS
seed = 0
models = load_models(
  files_meta=list(
    ANY=paste0("xgb_", imputation, "_any.meta"), 
    EAC=paste0("xgb_", imputation, "_eac.meta"), 
    EGJAC=paste0("xgb_", imputation, "_egjac.meta")
  ),
  files_models=list(
    ANY=paste0("xgb_", imputation, "_any.model"), 
    EAC=paste0("xgb_", imputation, "_eac.model"), 
    EGJAC=paste0("xgb_", imputation, "_egjac.model")
  )
)
features = feature_groups()
# ------------------------------------------------------------------------------




# ==============================================================================
# PREPARE DATA
df = imputed_df
set.seed(seed)
# subsample for SHAP/PDPs
df = imputed_df %>% sample_n(1000000)
# subsample controls
df = imputed_df %>% filter(casecontrol==1) %>% bind_rows(imputed_df%>%filter(casecontrol==0) %>% sample_n(100000))
# ------------------------------------------------------------------------------




# ==============================================================================
# Get SHAP & prepare
mname = "ANY"
model = models[[mname]]
xgb_df = xgb.DMatrix(as.matrix(df %>% select(model$feature_names)),
                     label=df$casecontrol)
proba = predict(model, newdata=xgb_df, predcontrib=TRUE, approxcontrib=F)

#by group
groups = lapply(unique(features$group), 
                function(gr) unique(features$name[features$group==gr]))
names(groups) = unique(features$group)
shap_groups = lapply(names(groups), function(gr){
  pr = proba[, groups[[gr]]]
  if(!is.vector(pr)) pr %<>% rowSums()
  pr
})
names(shap_groups) = names(groups)
shap_groups = bind_cols(shap_groups)

shap_groups %<>% mutate(
  anion_gap_eu = na + k + co2 + chlor,
  anion_gap_us = na + co2 + chlor,
)

# merge
colnames(shap_groups) = paste0("shap_", colnames(shap_groups))
shap_groups$id = df$id
var_shap = shap_groups %>% left_join(raw_df$df, by="id")

# Compute anion gap (from imputed data ...)
anion = df %>% mutate(
  anion_gap_eu_mean = na_mean + k_mean - co2_mean - chlor_mean,
  anion_gap_eu_min = na_min + k_min - co2_max - chlor_max,
  anion_gap_eu_max = na_max + k_max - co2_min - chlor_min,
  anion_gap_us_mean = na_mean - co2_mean - chlor_mean,
  anion_gap_us_min = na_min - co2_max - chlor_max,
  anion_gap_us_max = na_max - co2_min - chlor_min
)
var_shap %<>% left_join(anion %>% select(id, anion_gap_eu_mean, anion_gap_eu_min, anion_gap_eu_max,
                                         anion_gap_us_mean, anion_gap_us_min, anion_gap_us_max),
                        by="id")

# binary as factors
for(cname in colnames(var_shap)){
  if(stringr::str_sub(cname, 1L, 5L) == "shap_") next
  if(length(unique(var_shap[[cname]]))<=3){
    var_shap[[cname]] = as.factor(var_shap[[cname]])
  }
}

# obesity and ppi
var_shap %<>% mutate(
  ppi=factor(
    ifelse(is.na(ppi_int), "0", ifelse(ppi_mean<20, "<20", ">=20")),
    levels=c("0", "<20", ">=20")
  )
)
# ------------------------------------------------------------------------------







# ==============================================================================
# Graphical LASSO
covmat = var_shap %>% 
  select(starts_with("shap_")) %>% 
  select(-shap_anion_gap_us, -shap_anion_gap_eu) %>% 
  cov
grlasso_fit = glasso::glasso(covmat, rho=0.001)
prec = grlasso_fit$wi
rownames(prec) = rownames(covmat)
colnames(prec) = colnames(covmat)
prec = -prec/abs(prec)
diag(prec) = 1
g = ggcorrplot::ggcorrplot(prec, type="full")
g = g + ggtitle("GLasso on SHAP values aggregated by groups of features")
filename = paste0(dir_figures, "glasso_group.pdf")
ggsave(filename, g, width=10, height=10, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=10, height=10, bg="white")

regcoef_list = lapply(rownames(prec), function(vname){
  regcoef = -prec[vname, ] / prec[vname, vname]
  regcoef[vname] = 0
  return(regcoef[regcoef!=0])
})
names(regcoef_list) = rownames(prec)
for(vname in names(regcoef_list)){
  cat(vname, "\n")
  print(round(regcoef_list[[vname]], 3))
}
# ------------------------------------------------------------------------------







# ==============================================================================
# anion_gap vs SHAP anion gap
xvar = "wbc_mean"
yvar = "shap_anion_gap_eu"
xname = "WBC (mean)"
yname = "SHAP Anion Gap (log-odds)"
file = paste0(xvar, "_vs_", yvar)

wdf = var_shap %>% select(casecontrol, xvar, yvar) %>% tidyr::drop_na()
wdf[[yvar]] = wdf[[yvar]] - mean(wdf[[yvar]])

xrange = quantile(wdf %>% pull(get(xvar)), c(0.01, 0.99), na.rm=T)
yrange = quantile(wdf %>% pull(get(yvar)), c(0.1, 0.9), na.rm=T)
g_density = ggplot() +
  geom_histogram(
    data= wdf %>% group_by(casecontrol) %>% mutate(weight=1/n()),
    mapping=aes(x=get(xvar), color=casecontrol, weight=weight, fill=casecontrol),
    alpha=0.2,
    position="identity",
    binwidth=1
  ) + xlab(xname) + xlim(xrange) + ylab("Frequency")
g_shap = ggplot() +
  geom_smooth(
    data=wdf,
    mapping=aes(x=get(xvar), y=get(yvar)),
    method="gam"
  ) +
  geom_hline(yintercept=0) + xlab("") +
  theme(axis.text.x = element_blank()) + xlim(xrange) + 
  ylab(yname) +
  guides(color="none")
g = cowplot::plot_grid(g_shap, g_density, nrow=2, rel_heights=c(2, 1), align="v", axis="lr")

filename = paste0(dir_figures, file, ".pdf")
ggsave(filename, g, width=6, height=6, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=6, bg="white")

# ------------------------------------------------------------------------------





wdf = var_shap %>% filter(casecontrol==1) %>% bind_rows(var_shap%>%filter(casecontrol==0) %>% sample_n(2500))

g = GGally::ggpairs(
  wdf,
  columns=c("na_mean", "co2_mean", "chlor_mean", "bun_mean", "wbc_mean", "alt_mean", "age", "anion_gap_eu_mean"),
  mapping=aes(color=ppi),
  lower=list(continuous="smooth")
)
filename = paste0(dir_figures, "pairs_byppi.pdf")
ggsave(filename, g, width=16, height=16, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=16, height=16, bg="white")


g = GGally::ggpairs(
  wdf,
  columns=c("shap_na", "shap_co2", "shap_chlor", "shap_bun", "shap_wbc", "shap_alt", "shap_age", "shap_anion_gap_eu"),
  mapping=aes(color=ppi),
  lower=list(continuous="smooth")
)
filename = paste0(dir_figures, "shap_pairs_nbyppi.pdf")
ggsave(filename, g, width=16, height=16, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=16, height=16, bg="white")


wdf %<>% mutate(
  co2_bin = cut(co2_mean, c(0, 25, 29, 50))
)

wdf %<>% filter(!is.na(co2_bin))



g = GGally::ggpairs(
  wdf,
  columns=c("shap_na", "shap_co2", "shap_bun", "shap_wbc", "shap_age", "shap_chlor", "shap_creat"),
  mapping=aes(color=co2_bin)
)
filename = paste0(dir_figures, "shap_pairs_by_co2.pdf")
ggsave(filename, g, width=14, height=14, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=14, height=14, bg="white")

g = GGally::ggpairs(
  wdf,
  columns=c("na_mean", "co2_mean", "bun_mean", "wbc_mean", "age", "chlor_mean", "creat_mean"),
  mapping=aes(color=co2_bin)
)
filename = paste0(dir_figures, "pairs_by_co2.pdf")
ggsave(filename, g, width=14, height=14, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=14, height=14, bg="white")



# ==============================================================================
# co2
wdf = var_shap %>% select(c(
  "casecontrol", 
  "shap_bmi_weight", "shap_gerd", "shap_ppi", "shap_anion_gap_eu",
  "shap_bun", "shap_co2", "shap_na", "shap_chlor", "shap_k", "shap_alt",
  "bmi", "weight", "gerd", "bun_mean", "co2_mean", "na_mean", "k_mean", "chlor_mean", "alt_mean",
  "anion_gap_eu_mean"
))
wdf %<>% mutate(
  co2_bin = cut(co2_mean, c(0, 25, 29, 50))
)


ggplot() + 
  geom_smooth(
    data=wdf,
    mapping=aes(x=anion_gap_eu_mean, y=exp(shap_anion_gap_eu), color=co2_bin)
  ) + xlim(wdf$anion_gap_eu_mean %>% quantile(c(0.01, 0.99), na.rm=T)) + 
  geom_hline(yintercept=1)

# ------------------------------------------------------------------------------


