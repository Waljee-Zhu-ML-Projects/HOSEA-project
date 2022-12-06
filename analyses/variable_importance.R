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
setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
dir_imputed_data = "./R_data/imputed_records/"
dir_raw_data = "./R_data/processed_records/"
dir_figures = paste0("./R_code/hosea-project/figures/", imputation, "/variable_importance/")
dir_tables = paste0("./R_code/hosea-project/tables/", imputation, "/variable_importance/")
dir_shap = paste0("./R_code/hosea-project/figures/", imputation, "/shap/")
dir_pdp = paste0("./R_code/hosea-project/figures/", imputation, "/pdp/")
dir_models = "./R_data/results/models/"
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
  ggtitle("Gain Variable importance by category") +
  xlim(0, 100) + ylab("Category") + xlab("Gain (%)")
filename = paste0(dir_figures, "vi_cat.pdf")
ggsave(filename, g, width=6, height=4, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=4, bg="white")

vi_group_long = vi_group %>% tidyr::pivot_longer(names(models), names_to="Model", values_to="Gain")
group_order = vi_group_long %>% filter(Model=="ANY") %>% arrange(Gain) %>% pull(group)
vi_group_long$group %<>% factor(levels=group_order, ordered=T)

g = ggplot(data=vi_group_long, aes(y=reorder(group, Gain), x=Gain*100, fill=Model)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle("Gain Variable importance by group") +
  ylab("Feature group") + xlab("Gain (%)")
filename = paste0(dir_figures, "vi_group.pdf")
ggsave(filename, g, width=6, height=8, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=8, bg="white")

for(mname in names(models)){
  g = ggplot(data=vi_group_long %>% filter(Model==!!mname), aes(y=reorder(group, Gain), x=Gain*100)) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle(paste0("Gain Variable importance by group (", mname, ")")) +
    ylab("Feature group") + xlab("Gain (%)")
  filename = paste0(dir_figures, "vi_group_", mname, ".pdf")
  ggsave(filename, g, width=6, height=8, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=8, bg="white")
}


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
  ggtitle("SHAP Variable importance by category") +
  ylab("Category") + xlab("mean|SHAP|")
filename = paste0(dir_figures, "shap_cat.pdf")
ggsave(filename, g, width=6, height=4, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g,width=6, height=4, bg="white")

vi_group_long = vi_group %>% tidyr::pivot_longer(names(models), names_to="Model", values_to="Gain")
group_order = vi_group_long %>% filter(Model=="ANY") %>% arrange(Gain) %>% pull(group)
vi_group_long$group %<>% factor(levels=group_order, ordered=T)

g = ggplot(data=vi_group_long, aes(y=reorder(group, Gain), x=Gain, fill=Model)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle("SHAP Variable importance by group") +
  ylab("Feature group") + xlab("mean|SHAP|")
filename = paste0(dir_figures, "shap_group.pdf")
ggsave(filename, g, width=6, height=8, bg="white")
ggsave(stringr::str_replace(filename, "pdf", "png"), g,width=6, height=8, bg="white")

for(mname in names(models)){
  g = ggplot(data=vi_group_long %>% filter(Model==!!mname), aes(y=reorder(group, Gain), x=Gain)) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle(paste0("SHAP Variable importance by group (", mname, ")")) +
    ylab("Feature group") + xlab("mean|SHAP|")
  filename = paste0(dir_figures, "shap_group_", mname, ".pdf")
  ggsave(filename, g, width=6, height=8, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g,width=6, height=8, bg="white")
}

# ------------------------------------------------------------------------------




# ==============================================================================
# SHAP Marginals
set.seed(seed)
# subsample for SHAP/PDPs
df = imputed_df %>% sample_n(100000)

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
  wdf %<>% mutate(cancertype=ifelse(cancertype=="", "Control", cancertype))
  colnames(wdf) = c("id", "Model", "SHAP", "Var", "CancerType")
  bin = (wdf %>% pull(Var)%>% unique() %>% length()) <3
  xrange = quantile(wdf %>% pull(Var), c(0.001, 0.999))
  if(bin){
    g_density = ggplot() +
      geom_histogram(
        data=wdf %>% group_by(CancerType) %>% mutate(weight=1/n()),
        mapping=aes(Var, fill=CancerType, weight=weight),
        position="dodge",
        stat="count",
        width=0.1
      ) +
      xlab(var) + ylab("Frequency") + scale_x_continuous(breaks=0:1)
    xrange=c(-0.05, 1.05)
  }else{
    r = wdf %>% pull(Var) %>% range
    g_density = ggplot() +
      geom_histogram(
        data=wdf %>% group_by(CancerType) %>% mutate(weight=1/n()),
        mapping=aes(x=Var, color=CancerType, weight=weight, fill=CancerType),
        alpha=0.2,
        position="identity",
        binwidth=ifelse(r[2]-r[1] < 200, 1, (r[2]-r[1])/20)
      ) + xlab(var) + xlim(xrange) + ylab("Frequency")
  }
  g_shap = ggplot() +
    geom_smooth(
      data=wdf,
      mapping=aes(x=Var, y=exp(SHAP), color=Model),
      method=ifelse(bin, "lm", "gam"),
      se=F
    ) +
    geom_hline(yintercept=1) + xlab("") +
    theme(axis.text.x = element_blank()) + xlim(xrange)
  g = cowplot::plot_grid(g_shap, g_density, nrow=2, rel_heights=c(2, 1), align="v")
  filename = paste0(dir_shap, var, ".pdf")
  ggsave(filename, g, width=6, height=6, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=6, bg="white")
}

# ------------------------------------------------------------------------------







# ==============================================================================
# SHAP Correlation

# aggregate by group
groups = lapply(unique(features$group), 
                function(gr) unique(features$name[features$group==gr]))
names(groups) = unique(features$group)
shap_groups = lapply(names(groups), function(gr){
  pr = shap[, groups[[gr]]]
  if(!is.vector(pr)) pr %<>% rowSums()
  pr
})
names(shap_groups) = names(groups)
shap_groups = bind_cols(shap_groups)
shap_groups$Model = shap$Model

# aggregate by category
categories = lapply(unique(features$category), 
                    function(cat) unique(features$name[features$category==cat]))
names(categories) = unique(features$category)
shap_categories = lapply(names(categories), function(cat){
  pr = shap[, categories[[cat]]]
  if(!is.vector(pr)) pr %<>% rowSums()
  pr
})
names(shap_categories) = names(categories)
shap_categories = bind_cols(shap_categories)
shap_categories$Model = shap$Model

for(mname in names(models)){
  corr_groups = cor(shap_groups %>% filter(Model==!!mname) %>% select(-Model))
  corr_groups_categories = cor(
    shap_groups %>% filter(Model==!!mname) %>% select(-Model),
    shap_categories %>% filter(Model==!!mname) %>% select(-Model)
    )
  corr_lab_comorbibidities = cor(
    shap_groups %>% filter(Model==!!mname) %>% select(one_of(features %>% filter(category=="Lab") %>% pull(group))),
    shap_groups %>% filter(Model==!!mname) %>% select(one_of(features %>% filter(category=="Comorbidities") %>% pull(group)))
  )
  
  filename = paste0(dir_figures, "shap_corr_lab_comorbidities_", mname, ".pdf")
  g = ggcorrplot::ggcorrplot(
    corr_lab_comorbibidities,
    method="square"
  ) + ggtitle(
    paste0("Correlation between Lab and Comorbidities group SHAP (", mname, ")")
  )
  ggsave(filename, g, width=8, height=4, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=8, height=4, bg="white")
  
  filename = paste0(dir_figures, "shap_corr_groups_", mname, ".pdf")
  g = ggcorrplot::ggcorrplot(
    corr_groups,
    method="square"
  ) + ggtitle(
    paste0("Correlation between group SHAP (", mname, ")")
  )
  ggsave(filename, g, width=8, height=8, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=8, height=8, bg="white")
  
  filename = paste0(dir_figures, "shap_corr_groups_categories_", mname, ".pdf")
  g = ggcorrplot::ggcorrplot(
    corr_groups_categories,
    method="square"
  ) + ggtitle(
    paste0("Correlation between group and category SHAP (", mname, ")")
  )
  ggsave(filename, g, width=8, height=3, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=8, height=3, bg="white")
}

# ------------------------------------------------------------------------------







# ==============================================================================
# PDPs
df = imputed_df

models_rds = lapply(
  names(models),
  function(mname) readRDS(paste0(system.file('extdata', package = 'HOSEA'), "/xgb_", 
                                 imputation, "_", tolower(mname), ".meta"))
)
names(models_rds) = names(models)


for(var in features$name[111:174]){
  cat(paste0(var, "\n"))
  wdf = df %>% select(id, !!var) %>% left_join(raw_df$master %>% select(id, cancertype), by="id")
  bin = (wdf %>% pull(var)%>% unique() %>% length()) < 3
  wdf %<>% mutate(cancertype=ifelse(cancertype=="", "Control", cancertype))
  xrange = quantile(wdf %>% pull(var), c(0.001, 0.999), na.rm=T)
  colnames(wdf) = c("id", var, "CancerType")
  
  if(bin){
    g_density = ggplot() + 
      geom_histogram(
        data=wdf %>% group_by(CancerType) %>% mutate(weight=1/n()),
        mapping=aes(get(var), fill=CancerType, weight=weight),
        position="dodge",
        stat="count",
        width=0.1
      ) +
      xlab(var) + ylab("Frequency") + scale_x_continuous(breaks=0:1)
    xrange=c(-0.05, 1.05)
  }else{
    r = xrange
    binwidth = NULL
    if(r[2]-r[1] < 100 & r[2]-r[1] >10) binwidth=1 
    g_density = ggplot() +
      geom_histogram(
        data=wdf %>% group_by(CancerType) %>% mutate(weight=1/n()),
        mapping=aes(x=get(var), color=CancerType, weight=weight, fill=CancerType),
        alpha=0.2,
        position="identity",
        binwidth=binwidth
      ) + xlab(var) + xlim(xrange) + ylab("Frequency")
  }
  
  set.seed(seed)
  wdf %<>% sample_n(100000)
  pdps = list()
  
  for(mname in names(models)){
    model = models_rds[[mname]]$xgb_fit
    out = pdp::partial(model, pred.var=var, 
                       train=df %>% select(model$feature_names), 
                       type="classification", prob=T, which.class=1,
                       plot=F, progress=T, approx=F, 
                       quantiles=!bin)
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
    theme(axis.text.x = element_blank()) + xlim(xrange) + 
    geom_hline(yintercept=1)
  if(bin) g_pdp = g_pdp + scale_x_continuous(breaks=0:1)
  
  g = cowplot::plot_grid(g_pdp, g_density, nrow=2, rel_heights=c(2, 1), align="v")
  filename = paste0(dir_pdp, var, ".pdf")
  ggsave(filename, g, width=6, height=6, bg="white")
  ggsave(stringr::str_replace(filename, "pdf", "png"), g, width=6, height=6, bg="white")
  
}


# ------------------------------------------------------------------------------
