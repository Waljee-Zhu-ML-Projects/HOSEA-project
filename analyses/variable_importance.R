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
setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = "./R_code/hosea-project/figures/variable_importance/"
dir_shap = "./R_code/hosea-project/figures/shap/"
dir_pdp = "./R_code/hosea-project/figures/pdp/"
dir_tables = "./R_code/hosea-project/tables/variable_importance/"
dir_models = "./R_data/results/models/"
imputed_data = "test_mice_any.rds"
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
models = load_models()
features = feature_groups()
# ------------------------------------------------------------------------------




# ==============================================================================
# PREPARE DATA
set.seed(seed)
# subsample for SHAP/PDPs
df = imputed_df %>% sample_n(1000000)
# ------------------------------------------------------------------------------




# ==============================================================================
# VARIABLE IMPORTANCE

vis = list(
  feature=list(),
  group=list(),
  category=list()
)

for(mname in names(models)){
  
  model = models[[mname]]
  vi_raw = xgboost::xgb.importance(model$feature_names, model)
  vi = vi_raw %>% left_join(features, by=c("Feature"="name")) %>%
    select(-c(Cover, Frequency))
  vi %<>% replace(is.na(.), 0)
  vi_cat = vi %>% group_by(category) %>% summarize(Gain=sum(Gain))
  vi_group = vi %>% group_by(group) %>% summarize(Gain=sum(Gain))
  vi_group %<>% left_join(
    features %>% group_by(group, category) %>%
      summarize(group=first(group), category=first(category)), 
    by="group")
  vi_group %<>% arrange(category)
  vi_group = vi_group[, c("category", "group", "Gain")]
  vi_cat = vi_cat[, c("category", "Gain")]
  
  vi %<>% rename(!!mname:=Gain)
  vi_group %<>% rename(!!mname:=Gain)
  vi_cat %<>% rename(!!mname:=Gain)
  
  vis$feature[[mname]] = vi
  vis$group[[mname]] = vi_group
  vis$category[[mname]] = vi_cat
  
}

vi = vis$feature %>% purrr::reduce(full_join, by=c("Feature", "group", "category")) %>% replace(is.na(.), 0)
vi_group = vis$group %>% purrr::reduce(full_join, by=c("group", "category")) %>% replace(is.na(.), 0)
vi_cat = vis$category %>% purrr::reduce(full_join, by=c("category")) %>% replace(is.na(.), 0)

vi %<>% select((any_of(c("category", "group", "Feature", names(models)))))
vi_group %<>% select((any_of(c("category", "group", "Feature", names(models)))))
vi_cat %<>% select((any_of(c("category", "group", "Feature", names(models)))))

write.csv(vi, paste0(dir_tables, "vi.csv"))
write.csv(vi_group, paste0(dir_tables, "vi_group.csv"))
write.csv(vi_cat, paste0(dir_tables, "vi_cat.csv"))

vi_cat_long = vi_cat %>% tidyr::pivot_longer(names(models), names_to="Model", values_to="Gain")
category_order = vi_cat_long %>% filter(Model=="ANY") %>% arrange(Gain) %>% pull(category)
vi_cat_long$category %<>% factor(levels=category_order, ordered=T)

g = ggplot(data=vi_cat_long, aes(y=reorder(category, Gain), x=Gain*100, fill=Model)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle("Variable importance by category") +
  xlim(0, 100) + ylab("Category") + xlab("Gain (%)")
filename = paste0(dir_figures, "vi_cat.pdf")
ggsave(filename, g, width=6, height=4)

vi_group_long = vi_group %>% tidyr::pivot_longer(names(models), names_to="Model", values_to="Gain")
group_order = vi_group_long %>% filter(Model=="ANY") %>% arrange(Gain) %>% pull(group)
vi_group_long$group %<>% factor(levels=group_order, ordered=T)

g = ggplot(data=vi_group_long, aes(y=reorder(group, Gain), x=Gain*100, fill=Model)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle("Variable importance by group") +
  ylab("Feature group") + xlab("Gain (%)")
filename = paste0(dir_figures, "vi_group.pdf")
ggsave(filename, g, width=6, height=8)



# ------------------------------------------------------------------------------




# ==============================================================================
# SHAP

vis = list(
  feature=list(),
  group=list(),
  category=list(),
  raw=list()
)

for(mname in names(models)){
  
  model = models[[mname]]
  
  xgb_df = xgb.DMatrix(as.matrix(df %>% select(model$feature_names)),
                        label=df$casecontrol)
  proba = predict(model, newdata=xgb_df, predcontrib=TRUE, approxcontrib=F)
  
  vi = colSums(abs(proba))
  vi_all = data.frame(feature=names(vi), SHAP=vi)
  
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
  shap_group_agg = abs(shap_groups) %>% colMeans()
  vi_group = data.frame(group=names(shap_group_agg), SHAP=shap_group_agg)
  
  #by category
  categories = lapply(unique(features$category), 
                      function(cat) unique(features$name[features$category==cat]))
  names(categories) = unique(features$category)
  shap_categories = lapply(names(categories), function(cat){
    pr = proba[, categories[[cat]]]
    if(!is.vector(pr)) pr %<>% rowSums()
    pr
  })
  names(shap_categories) = names(categories)
  shap_categories = bind_cols(shap_categories)
  shap_categories_agg = abs(shap_categories) %>% colMeans()
  vi_cat = data.frame(category=names(shap_categories_agg), SHAP=shap_categories_agg)
  
  vi_group %<>% left_join(
    features %>% group_by(group, category) %>%
      summarize(group=first(group), category=first(category)), 
    by="group")
  
  vi_all %<>% left_join(features, by=c("feature"="name")) %>% rename(Feature=feature)
  
  vi_all %<>% rename(!!mname:=SHAP)
  vi_group %<>% rename(!!mname:=SHAP)
  vi_cat %<>% rename(!!mname:=SHAP)
  
  
  proba %<>% bind_cols(df %>% select(id, casecontrol))
  vis$feature[[mname]] = vi_all
  vis$group[[mname]] = vi_group
  vis$category[[mname]] = vi_cat
  vis$raw[[mname]] = proba
  
}

vi = vis$feature %>% purrr::reduce(full_join, by=c("Feature", "group", "category")) %>% replace(is.na(.), 0)
vi_group = vis$group %>% purrr::reduce(full_join, by=c("group", "category")) %>% replace(is.na(.), 0)
vi_cat = vis$category %>% purrr::reduce(full_join, by=c("category")) %>% replace(is.na(.), 0)

vi %<>% select((any_of(c("category", "group", "Feature", names(models)))))
vi_group %<>% select((any_of(c("category", "group", "Feature", names(models)))))
vi_cat %<>% select((any_of(c("category", "group", "Feature", names(models)))))

write.csv(vi, paste0(dir_tables, "shap.csv"))
write.csv(vi_group, paste0(dir_tables, "shap_group.csv"))
write.csv(vi_cat, paste0(dir_tables, "shap_cat.csv"))

vi_cat_long = vi_cat %>% tidyr::pivot_longer(names(models), names_to="Model", values_to="Gain")
category_order = vi_cat_long %>% filter(Model=="ANY") %>% arrange(Gain) %>% pull(category)
vi_cat_long$category %<>% factor(levels=category_order, ordered=T)

g = ggplot(data=vi_cat_long, aes(y=reorder(category, Gain), x=Gain, fill=Model)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle("Variable importance by category") +
  ylab("Category") + xlab("mean|SHAP|")
filename = paste0(dir_figures, "shap_cat.pdf")
ggsave(filename, g, width=6, height=4)

vi_group_long = vi_group %>% tidyr::pivot_longer(names(models), names_to="Model", values_to="Gain")
group_order = vi_group_long %>% filter(Model=="ANY") %>% arrange(Gain) %>% pull(group)
vi_group_long$group %<>% factor(levels=group_order, ordered=T)

g = ggplot(data=vi_group_long, aes(y=reorder(group, Gain), x=Gain, fill=Model)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle("Variable importance by group") +
  ylab("Feature group") + xlab("mean|SHAP|")
filename = paste0(dir_figures, "shap_group.pdf")
ggsave(filename, g, width=6, height=8)

# ------------------------------------------------------------------------------




# ==============================================================================
# SHAP Marginals

shaps = list()

for(mname in names(models)){
  model = models[[mname]]
  xgb_df = xgb.DMatrix(as.matrix(df %>% select(model$feature_names)),
                       label=df$casecontrol)
  proba = predict(model, newdata=xgb_df, predcontrib=TRUE, approxcontrib=F)
  proba %<>% bind_cols(df %>% select(id))
  proba$Model = mname
  shaps[[mname]] = proba
}

shap = shaps %>% bind_rows()

variables = shap %>% select(-c(id, Model, BIAS)) %>% colnames()

for(var in variables){
  wdf = shap %>% select(id, Model, !!var)
  wdf %<>% left_join(df %>% select(id, !!var), by="id")
  wdf %<>% left_join(raw_df$master %>% select(id, cancertype), by="id")
  colnames(wdf) = c("id", "Model", "SHAP", "Var", "CancerType")
  bin = (wdf %>% pull(Var)%>% unique() %>% length()) <3
  xrange = quantile(wdf %>% pull(Var), c(0.001, 0.999))
  g_shap = ggplot() + 
    geom_smooth(
      data=wdf,
      mapping=aes(x=Var, y=exp(SHAP), color=Model),
      method=ifelse(bin, "lm", "gam")
    ) + 
    geom_hline(yintercept=1) + xlab("") +
    theme(axis.text.x = element_blank()) + xlim(xrange)
  g_density = ggplot() + 
    geom_density(
      data=wdf,
      mapping=aes(Var, fill=CancerType, color=CancerType),
      alpha=0.2
    ) +
    xlab(var) + xlim(xrange)
  g = cowplot::plot_grid(g_shap, g_density, nrow=2, rel_heights=c(2, 1), align="v")
  filename = paste0(dir_shap, var, ".pdf")
  ggsave(filename, g, width=6, height=6)
}

# ------------------------------------------------------------------------------







# ==============================================================================
# PDPs

models_rds = lapply(names(models), function(mname) readRDS(paste0(dir_models, "XGB_all_", mname, ".rds")))
names(models_rds) = names(models)

set.seed(seed)
df_ = df %>% sample_n(100000)


for(var in models_rds[["ANY"]]$xgb_fit$feature_names){
  cat(paste0(var, "\n"))
  wdf = df %>% select(id, !!var) %>% left_join(raw_df$master %>% select(id, cancertype), by="id")
  bin = (wdf %>% pull(var)%>% unique() %>% length()) < 3
  xrange = quantile(wdf %>% pull(var), c(0.001, 0.999), na.rm=T)
  colnames(wdf) = c("id", var, "CancerType")
  
  g_density = ggplot() + 
    geom_density(
      data=wdf,
      mapping=aes(get(var), fill=CancerType, color=CancerType),
      alpha=0.2
    ) +
    xlab(var) + xlim(xrange)
  
  
  pdps = list()
  
  for(mname in names(models)){
    model = models_rds[[mname]]
    out = pdp::partial(model$xgb_fit, pred.var=var, 
                       train=df_ %>% select(model$xgb_fit$feature_names), 
                       type="classification", prob=T, which.class=1,
                       plot=F, progress="text")
    out$yhat = out$yhat * 100000
    out$odds = out$yhat / mean(out$yhat)
    
    pdps[[mname]] = out
  }
  
  pdps %<>% bind_rows(.id="Model")
  
  g_pdp = ggplot(
    data=pdps,
    mapping=aes(x=get(var), y=odds, color=Model)
  ) + geom_line() + 
    ylab("Odds") + xlab("") +
    theme(axis.text.x = element_blank()) + xlim(xrange)
  
  g = cowplot::plot_grid(g_pdp, g_density, nrow=2, rel_heights=c(2, 1), align="v")
  filename = paste0(dir_pdp, var, ".pdf")
  ggsave(filename, g, width=6, height=6)
  
}


# ------------------------------------------------------------------------------
      