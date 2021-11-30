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

n0all = (complete_data$CaseControl==0)%>%sum
n1all = (complete_data$CaseControl==1)%>%sum

# replciability
set.seed(0)
n_quantiles = 10000

# logging
log = function(df) cat(paste("Full data set: ", nrow(df), "observations,", 
                             df$CaseControl%>%sum, "cases", 
                             (df$CaseControl==0)%>%sum, "controls"), fill=T)


# train-test split
n_controls = 2e6
complete_data = subsample_controls(complete_data, n_controls)
out = train_test_split(df=complete_data, weights=c(3, 1))
train = out[[1]]
test = out[[2]]
rm(complete_data)
gc()

log(train)
log(test)

# prepare stuff
methods = c(
  "resample",
  "unweighted"
  # "weighted",
  # "weighted_sqrt"
)
param_xg = list(
  max_depth = 5,
  subsample = 0.2,
  eta = .05,
  objective = 'binary:logistic',
  eval_metric = 'auc'
)

# loop
n = "2M"
set.seed(0)
out = train_test_split(df=train, weights=c(2, 1))
rm(train);gc()
train_n = out[[1]]
valid_n = out[[2]]
log(train_n)
log(valid_n)
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
  log(train_)
  # imputation
  quantiles = compute_quantiles(train_n, n_quantiles)
  set.seed(0)
  train_ = impute_srs(train_, quantiles)
  test_ = impute_srs(test, quantiles)
  valid_ = impute_srs(valid_n, quantiles)
  log(train_);log(valid_);log(test_)
  dwatchlist = xgb_prep(train=train_,
                        test=test_,
                        valid=valid_,
                        cc=cc_test)
  dwatchlist$test = NULL
  dwatchlist$cc = NULL
  gc()
  # fit
  set.seed(0)
  param_xg$scale_pos_weight = switch(
    method,
    "downsample" = 1.,
    "resample" = (n1all/n0all)*(n0/n0),
    "unweighted" = (n1all/n0all)*(n0/n1),
    "weighted" = n0/n1,
    "weighted_sqrt" = sqrt(n0/n1),
    "weighted_sq" = (n0/n1)^2
  )
  xgb_fit = xgb.train(param_xg,
                      dwatchlist$train,
                      nrounds=10000,
                      dwatchlist,
                      verbose=1,print_every_n=100,
                      early_stopping_rounds=100)
  # prepare testing sets
  test_sets = evaluation_split(test, test_)
  test_sets_groups = split_by_vargoups(test, test_)
  test_sets = append(test_sets, test_sets_groups)
  rm(test_sets_groups); gc()
  test_set_summaries = t(rbind(
    N=sapply(test_sets, nrow),
    Ncases=sapply(test_sets, function(df) sum(getinfo(df, "label"))),
    propcases=sapply(test_sets, function(df) mean(getinfo(df, "label")))
  ))
  # evaluate
  calibration_metrics = lapply(test_sets, function(df) calibration(xgb_fit, df))
  thresholds = c(seq(0.0001, 0.00098, 0.00002), 
                 seq(0.001, 0.0098, 0.0002), 
                 seq(0.01, 0.098, 0.002), 
                 seq(0.10, 0.99, 0.01))
  threshold_metrics = lapply(test_sets, function(df) 
    classification_metrics(xgb_fit, df, thresholds))
  calibration_curves = lapply(test_sets, function(df) calibration_curve(xgb_fit, df))
  # save
  out = list(
    xgb_fit=xgb_fit,
    metrics=list(calibration=calibration_metrics, 
                 classification=threshold_metrics,
                 calibration_curves=calibration_curves),
    test_sets=test_set_summaries
  )
  filepath = paste0("R_data/results/models/", method, "_n", n, ".rds")
  saveRDS(out, filepath)
}


out=readRDS(filepath)
xgb_fit=out$xgb_fit

curve = calibration_curves[[1]]
plot(curve$mid, curve$propcase, log="xy", ylim=c(1e-5, 1.), xlim=c(1e-5, 1.))
abline(0, 1)

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
  tr_metrics_df$tr = as.numeric(rownames(tr_metrics_df))
  tr_metrics_df$test_df = names(tr_metrics)[1]
  for(df in names(tr_metrics)[-1]){
    tr_metrics[[df]]$tr = as.numeric(rownames(tr_metrics[[df]]))
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
tmethod = "resample"
tcurve = "pr"
# filter & melt
aucs = auc_df %>% filter(method == tmethod)
aucs = aucs %>% filter(curve == tcurve)
dfs = c("all", "cc", "0_5", "5_10", "10_30", "30_100")

aucs = as.data.frame(aucs %>% tidyr::pivot_longer(cols=dfs))
colnames(aucs) = c("curve", "n", "method", "test_df", "auc")

pdf(paste0("R_code/hosea-project/figures/size_", tmethod, "_", tcurve, ".pdf"), width=8, height=5)
ggplot(aucs, aes(n, auc, group=test_df, colour=test_df)) +
  geom_line() + scale_x_continuous(trans="log10") +
  xlab("Nb. controls") + ylab(paste("auc", tcurve)) +
  ggtitle(paste0("Method: ", tmethod))
dev.off()

# threshold metrics
colnames(tr_df)
metric = "detection_prevalance"
tmethod = "resample"
ttest_df = "all"
# filter
trm = tr_df %>% filter(method==tmethod)
trm = trm %>% filter(test_df==ttest_df)

pdf(paste0("R_code/hosea-project/figures/size_", tmethod, "_", ttest_df, "_", metric, ".pdf"), 
    width=8, height=5)
ggplot(trm, aes(tr, detection_prevalance, group=n, colour=n)) +
  geom_line() +
  xlab("Threshold") + ylab(metric) +
  ggtitle(paste0("Method: ", tmethod, ", Test df: ", ttest_df))
dev.off()
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
