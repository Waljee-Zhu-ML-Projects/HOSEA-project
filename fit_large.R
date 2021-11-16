setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
source('R_code/hosea-project/setup.R')
source('R_code/hosea-project/utils_subsample.R')
source('R_code/hosea-project/impute_missing.R')
source('R_code/hosea-project/utils_xgb.R')
source('R_code/hosea-project/compute_quantiles.R')
source('R_code/hosea-project/classification_metrics.R')
source('R_code/hosea-project/evaluation_split.R')

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

# prepare stuff
ns = c(1e6)
# ns = c(1e4, 2e4, 5e4, 1e5, 2e5, 5e5)
methods = c(
  # "resample",
  # "unweighted",
  # "weighted",
  "weighted_sqrt"
)
param_xg = list(
  max_depth = 5,
  subsample = 0.5,
  eta = .05,
  objective = 'binary:logistic',
  eval_metric = 'auc'
)

# loop
for(n in ns){
  sub_train = subsample_controls(train, n)
  log(sub_train)
  set.seed(0)
  out = train_test_split(df=sub_train, weights=c(2, 1))
  train_n = out[[1]]
  valid_n = out[[2]]
  n0 = (train_n$CaseControl==0)%>%sum
  n1 = (train_n$CaseControl==1)%>%sum
  for(method in methods){
    cat("===================================================================", fill=T)
    cat(paste0("===== method = ", method, ", n = ", n), fill=T)
    cat("===================================================================", fill=T)
    # produce training set
    train_ = switch(method,
                   "downsample" = subsample_controls(train_n, n1),
                   "resample" = balanced_resample(train_n),
                   "unweighted" = train_n,
                   "weighted" = balanced_resample(train_n),
                   "weighted_sqrt" = balanced_resample(train_n),
                   "weighted_sq" = balanced_resample(train_n)
    )
    # imputation
    quantiles = compute_quantiles(train_n, n_quantiles)
    set.seed(0)
    train_ = impute_srs(train_, quantiles)
    test_ = impute_srs(test, quantiles)
    valid_ = impute_srs(valid_n, quantiles)
    dwatchlist = xgb_prep(train=train_,
                          test=test_,
                          valid=valid_,
                          cc=cc_test)
    dwatchlist$test = NULL
    dwatchlist$cc = NULL
    # fit
    set.seed(0)
    param_xg$scale_pos_weight = switch(
      method,
      "downsample" = 1.,
      "resample" = 1.,
      "unweighted" = 1.,
      "weighted" = n0/n1,
      "weighted_sqrt" = sqrt(n0/n1),
      "weighted_sq" = (n0/n1)^2
    )
    xgb_fit = xgb.train(param_xg,
                        dwatchlist$train,
                        nrounds=5000,
                        dwatchlist,
                        verbose=1,print_every_n=50,
                        early_stopping_rounds=50)
    # prepare testing sets
    test_sets = evaluation_split(test, test_)
    # evaluate
    calibration_metrics = lapply(test_sets, function(df) calibration(xgb_fit, df))
    threshold_metrics = lapply(test_sets, function(df) 
      classification_metrics(xgb_fit, df, seq(0.01, 0.99, 0.01)))
    # save
    out = list(
      xgb_fit=xgb_fit,
      metrics=list(calibration=calibration_metrics, classification=threshold_metrics)
    )
    filepath = paste0("R_data/results/models/", method, "_n", n, ".rds")
    saveRDS(out, filepath)
  }
}

# plot results

library(ggplot2)
ns = c(1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6)
methods = c(
  "resample",
  "unweighted",
  "weighted",
  "weighted_sqrt"
)

# AUCs
results_to_df = function(method, n) {
  filepath = paste0("R_data/results/models/", method, "_n", n, ".rds")
  results = readRDS(filepath)
  auc_roc = sapply(results$metrics$calibration, function(x) x$roc$auc)
  auc_pr = sapply(results$metrics$calibration, function(x) x$pr$auc.integral)
  auc_df = data.frame(t(cbind(auc_roc, auc_pr)))
  colnames(auc_df) = c("all", "cc", "0_5", "5_10", "10_30", "30_100")
  auc_df$curve = c("roc", "pr")
  auc_df$n = n
  auc_df$method = method
  tr_metrics = results$metrics$classification
  tr_metrics_df = data.frame(tr_metrics[[1]])
  tr_metrics_df$test_df = names(tr_metrics)[1]
  for(df in names(tr_metrics)[-1]){
    tr_metrics[[df]]$test_df = df
    tr_metrics_df = rbind(tr_metrics_df, tr_metrics[[df]])
  }
  tr_metrics_df$n = n
  tr_metrics_df$method = method
  return(list(tr=tr_metrics_df, auc=auc_df))
}

auc_df = data.frame()
tr_df = data.frame()
for(n in ns){
  for(method in methods){
    dfs = results_to_df(method, n)
    auc_df = bind_rows(auc_df, dfs$auc)
    tr_df = bind_rows(tr_df, dfs$tr)
  }
}

# auc plots
#drop sq
aucs = aucs %>% filter(method != "weighted_sq")

ggplot(aucs, aes(n, test, group=interaction(method, tuning), colour=method, linetype=tuning)) +
  geom_line() + scale_x_continuous(trans="log10") +
  xlab("Nb. controls") + ylab("Test AUC")


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
