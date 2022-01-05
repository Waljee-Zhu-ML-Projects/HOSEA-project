setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(ggplot2)

n = "1e+05"
results = list()
resample = readRDS(paste0("R_data/results/models45/resample_n", n, ".rds"))
unweighted = readRDS(paste0("R_data/results/models/unweighted_n", n, ".rds"))


# merge results into common dfs
all = lapply(resample$metrics$classification, function(df){
  df$method = "resample"
  return(df)
})
for(nm in names(all)){
  # df = unweighted$metrics$classification[[nm]]
  # df$method = "unweighted"
  all[[nm]]$tr = as.numeric(rownames(all[[nm]]))
  # df$tr = as.numeric(rownames(df))
  # all[[nm]] = rbind(all[[nm]], df)
}
alldf = data.frame()
for(nm in names(all)){
  df = all[[nm]]
  df$df = nm
  alldf = bind_rows(alldf, df)
}


pdf("R_code/hosea-project/figures/all_ppv_allcurves.pdf", 8, 5)
plotdf = alldf
ggplot(data=plotdf, aes(x=tr, y=ppv, group=interaction(method, df), 
                       color=df, linetype=method)) + 
  geom_line() + scale_x_continuous(trans="log10")
dev.off()


pdf("R_code/hosea-project/figures/all_ppv_missing.pdf", 8, 5)
plotdf = alldf %>% filter(df  %in% c("all", "cc", "all_0_5", "all_5_10", "all_10_30",
                                     "all_30_100"))
ggplot(data=plotdf, aes(x=tr, y=ppv, group=interaction(method, df), 
                        color=df, linetype=method)) + 
  geom_line() + scale_x_continuous(trans="log10")
dev.off()

pdf("R_code/hosea-project/figures/all_ppv_vars.pdf", 8, 5)
plotdf = alldf %>% filter(!(df %in% c("all", "cc", "all_0_5", "all_5_10", "all_10_30",
                                     "all_30_100")))
ggplot(data=plotdf, aes(x=tr, y=ppv, group=interaction(method, df), 
                        color=df, linetype=method)) + 
  geom_line() + scale_x_continuous(trans="log10")
dev.off()


df = alldf %>% filter(df == "all")
df = df %>% select(one_of("tpr", "ppv", "detection_prevalance"))
df = df[46:225, ]
rownames(df) = format(round(as.numeric(rownames(df)), 5)*100000, digits=5)
df = df*100
cat(print(xtable::xtable(df)), file="R_code/hosea-project/figures/all_calibration.tex")


# calibration (AUC)
aucs = data.frame(matrix(0, 1, 3))
colnames(aucs) = c("df", "method", "auc")
for(df in names(resample$metrics$calibration)){
  aucs = rbind(aucs, 
               c(df=df, method="resample", 
                 auc=resample$metrics$calibration[[df]]$roc$auc))
}
for(df in names(unweighted$metrics$calibration)){
  aucs = rbind(aucs, 
               c(df=df, method="unweighted", 
                 auc=unweighted$metrics$calibration[[df]]$roc$auc))
}
aucs = aucs[-1, ]
aucs$auc = as.numeric(aucs$auc)

pdf("R_code/hosea-project/figures/25_all_aucs.pdf", 8, 5)
ggplot(data=aucs, aes(x=df, y=auc, fill=method)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  coord_cartesian(ylim=c(0.6, 0.88))
dev.off()

# ROC curves
rocs = data.frame(matrix(0, 1, 5))
colnames(rocs) = c("fpr", "tpr", "tr", "df", "method")

for(df in names(resample$metrics$calibration)){
  roc = data.frame(resample$metrics$calibration[[df]]$roc$curve)
  roc = roc[seq(1, nrow(roc), length.out=100), ]
  colnames(roc) = c("fpr", "tpr", "tr")
  roc$df = df
  roc$method = "resample"
  rocs = rbind(rocs, roc)
}

for(df in names(unweighted$metrics$calibration)){
  roc = data.frame(unweighted$metrics$calibration[[df]]$roc$curve)
  roc = roc[seq(1, nrow(roc), length.out=100), ]
  colnames(roc) = c("fpr", "tpr", "tr")
  roc$df = df
  roc$method = "unweighted"
  rocs = rbind(rocs, roc)
}
rocs = rocs[-1, ]

tdf = "all"
tmethod = "resample"
auroc = aucs[aucs$df==tdf & aucs$method==tmethod, "auc"]
filepath = paste0("R_code/hosea-project/figures/45_all_roc_", tmethod, "_", tdf, ".pdf")
plotdf = rocs %>% filter(df == tdf & method == tmethod)
colnames(plotdf) = c("fpr", "tpr", "Threshold", "df", "method")
g = ggplot(data=plotdf, aes(x=fpr, y=tpr, color=Threshold)) + geom_line() +
  scale_color_gradientn(colors=rainbow(10)) + theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted")
# dev.off()
ggsave(filepath, g, width=8, height=7)

for(tdf in names(resample$metrics$calibration)){
  for(tmethod in c("resample", "unweighted")){
    auroc = aucs[aucs$df==tdf & aucs$method==tmethod, "auc"]
    # pdf(paste0("R_code/hosea-project/figures/roc/", tmethod, "_", tdf, ".pdf"), 6, 5)
    plotdf = rocs %>% filter(df == tdf & method == tmethod)
    g = ggplot(data=plotdf, aes(x=fpr, y=tpr, color=tr)) + geom_line() +
      scale_color_gradientn(colors=rainbow(20)) + theme(aspect.ratio=1) +
      ggtitle(paste0("Df: ", tdf, ", method: ", tmethod, ", AUROC=", round(auroc, 4)))
    # dev.off()
    ggsave(paste0("R_code/hosea-project/figures/roc/45_all_", tmethod, "_", tdf, ".pdf"), g,
           width=6, height=5)
  }
}

# calibration plot

all = lapply(resample$metrics$calibration_curves, function(df){
  df$method = "resample"
  return(df)
})
for(nm in names(all)){
  # df = unweighted$metrics$calibration_curves[[nm]]
  # df$method = "unweighted"
  all[[nm]]$tr = as.numeric(rownames(all[[nm]]))
  # df$tr = as.numeric(rownames(df))
  # all[[nm]] = rbind(all[[nm]], df)
}
alldf = data.frame()
for(nm in names(all)){
  df = all[[nm]]
  df$df = nm
  alldf = bind_rows(alldf, df)
}


tdf = "all"
tmethod = "resample"
log = F

for(tdf in names(unweighted$metrics$calibration)){
  for(tmethod in c("resample", "unweighted")){
    for(log in c(T, F)){
      filename = paste0("R_code/hosea-project/figures/calibration_curves/25_all_", tmethod, "_", 
                 tdf,  ifelse(log, "_log", ""), ".pdf")
      plotdf = alldf %>% filter(df == tdf & method == tmethod)
      g = ggplot(data=plotdf, aes(x=mid, y=propcase)) + theme(aspect.ratio=1) + 
        geom_point()  +
        ggtitle(paste0("Df: ", tdf, ", method: ", tmethod)) +
        geom_abline(slope=1, intercept=0, linetype="dashed") +
        ylab("Observed") + xlab("Predicted")
      if(log){
        g = g + scale_x_log10(limits=c(1e-5, 1.)) + scale_y_log10(limits=c(1e-5, 1.))
      }else{
        g = g + xlim(0., 0.01) + ylim(0., 0.01)
      }
      ggsave(filename, g, width=5, height=5)
    }
  }
}

# Calibration metrics
df = resample$metrics$calibration_curves$all
pred = df$mid
obs = df$propcase
n = df$N
n1 = df$Ncase

cor.test(pred[1:50], obs[1:50], method="spearman")
cor.test(pred[1:50], obs[1:50], method="pearson")

o1 = n1
o0 = n-n1
e1 = pred*n
e0 = n-e1

H51 = sum(((o1-e1)^2/e1 + (o0-e0)^2/e0))
df51 = 49
p51 = pchisq(H51, df51, lower.tail=F)
H50 = sum(((o1-e1)^2/e1 + (o0-e0)^2/e0)[1:50])
df50 = 48
p50 = pchisq(H50, df50, lower.tail=F)

# variable importance
imp = xgboost::xgb.importance(model=resample$xgb_fit)
which.max(imp$Feature == "smoke_current")
which.max(imp$Feature == "smoke_former")
which.max(imp$Feature == "HIV")
