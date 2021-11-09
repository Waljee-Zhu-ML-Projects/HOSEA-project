source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')

# import data
complete_data = readRDS('R_data/complete_data_raw.rds')
complete_data$n_visits = NULL 
cc_test <- readRDS('R_data/cc_complete_data.rds')

# replciability
set.seed(0)
n_quantiles = 10000

# logging
log = function(df) cat(paste("Full data set: ", nrow(df), "observations,", 
                             df$CaseControl%>%sum, "cases", 
                             (df$CaseControl==0)%>%sum, "controls"), fill=T)

# subsampling of controls
n_controls = 2e6
sub_complete_data = subsample_controls(complete_data, n_controls)
log(sub_complete_data)
rm(complete_data)
gc()

# train-test split
out = train_test_split(df=sub_complete_data, weights=c(3, 1))
train = out[[1]]
test = out[[2]]

log(train)
log(test)

# loop over n from here
n = 1e5

# base XGBoost parameters
param_xg_base = list(
  max_depth=5,
  subsample=0.5,
  eta=.05,
  objective='binary:logistic',
  eval_metric='auc',
  scale_pos_weight = 1.
)


# grids
values = list(
  max_depth = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
  subsample = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.),
  eta = c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.)
)
for(param in names(values)){
  aucs = data.frame(matrix(0., 1, 5))
  colnames(aucs) = c("value", "train", "valid", "test", "cc")
  for(value in values[[param]]){
    param_xg = param_xg_base
    param_xg[[param]] = value
    # prepare output
    sub_train = subsample_controls(train, n)
    log(sub_train)
    out = train_test_split(df=sub_train, weights=c(2, 1))
    train_n = out[[1]]
    valid_n = out[[2]]
    log(train_n)
    log(valid_n)
    n0 = (train_n$CaseControl==0)%>%sum
    n1 = (train_n$CaseControl==1)%>%sum
    
    # resampling
    train_resample = balanced_resample(train_n)
    log(train_resample)
    
    quantiles_resample = compute_quantiles(train_resample, n_quantiles)
    set.seed(0)
    train_resample = impute_srs(train_resample, quantiles_resample)
    test_resample = impute_srs(valid_n, quantiles_resample)
    valid_resample = impute_srs(test, quantiles_resample)
    
    dwatchlist_resample = xgb_prep(train_resample,
                                   test_resample,
                                   valid_resample,
                                   cc=cc_test)
    
    # fit
    set.seed(1)
    xgb_fit_resample = xgb.train(param_xg,
                                 dwatchlist_resample$train,
                                 nrounds=5000,
                                 dwatchlist_resample,
                                 verbose=1,print_every_n=10,
                                 early_stopping_rounds=50)
    
    
    aucs = rbind(aucs, c(value, best_auc(xgb_fit_resample)))
    
    write.csv(aucs, paste0("R_data/results/best_aucs_tuning_", param,".csv"))
  }
}

# for(param in names(values)){
for(param in c("max_depth", "subsample")){
  aucs = read.csv(paste0("R_data/results/best_aucs_tuning_", param,".csv"))
  aucs$X = NULL
  aucs = aucs[-1, ]
  # to long format
  aucs_long = tidyr::pivot_longer(aucs, c("train", "valid", "test", "cc"), 
                                  names_to="dataset", values_to="auc")
  value_base = param_xg_base[[param]]
  
  pdf(paste0("R_code/hosea-project/figures/best_aucs_tuning_", param, ".pdf"), width=8, height=5)
  ggplot(aucs_long, aes(value, auc, group=dataset, colour=dataset)) + 
    geom_line() + geom_vline(xintercept=value_base, linetype="dashed") +
    xlab(param) + ylab("AUC")
  dev.off()
}

# aucs = read.csv("R_data/results/best_aucs_large.csv")
# aucs$X = NULL
# aucs$tuning = "old"
# aucs = aucs[-1, ]
# 
# for(n in c(1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6)){
#   aucs_n = read.csv(paste0("R_data/results/best_aucs_large_", n, ".csv"))
#   aucs_n$X = NULL
#   aucs_n$tuning = "new"
#   aucs_n = aucs_n[-1, ]
#   aucs = rbind(aucs, aucs_n)
# }
# 
# aucs = aucs[order(aucs$n), ]
# aucs = aucs[order(aucs$method), ]
# aucs = aucs[order(aucs$type), ]
# xtable::xtable(aucs, digits=3, include.rownames=F)
# 
# library(ggplot2)
# #drop sq
# aucs = aucs %>% filter(method != "weighted_sq")
# 
# ggplot(aucs, aes(n, test, group=interaction(method, tuning), colour=method, linetype=tuning)) +
#   geom_line() + scale_x_continuous(trans="log10") +
#   xlab("Nb. controls") + ylab("Test AUC")
# 
# 
# ggplot(aucs, aes(n, cc, group=interaction(method, type), colour=method, linetype=type)) +
#   geom_line() + scale_x_continuous(trans="log10") +
#   xlab("Nb. controls") + ylab("CC AUC")
