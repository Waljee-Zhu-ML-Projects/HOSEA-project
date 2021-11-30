setwd('/nfs/turbo/umms-awaljee/umms-awaljee-HOSEA/Peter files')
library(dplyr)
library(ggplot2)

n = "2M"
results = list()
resample = readRDS(paste0("R_data/results/models/resample_n", n, ".rds"))
unweighted = readRDS(paste0("R_data/results/models/unweighted_n", n, ".rds"))


# merge results into common dfs
all = lapply(resample$metrics$classification, function(df){
  df$method = "resample"
  return(df)
})
for(nm in names(all)){
  df = unweighted$metrics$classification[[nm]]
  df$method = "unweighted"
  all[[nm]]$tr = as.numeric(rownames(all[[nm]]))
  df$tr = as.numeric(rownames(df))
  all[[nm]] = rbind(all[[nm]], df)
}
alldf = data.frame()
for(nm in names(all)){
  df = all[[nm]]
  df$df = nm
  alldf = bind_rows(alldf, df)
}


pdf("R_code/hosea-project/figures/dp_allcurves.pdf", 8, 5)
plotdf = alldf
ggplot(data=plotdf, aes(x=tr, y=detection_prevalance, group=interaction(method, df), 
                       color=df, linetype=method)) + 
  geom_line() + scale_x_continuous(trans="log10")
dev.off()


pdf("R_code/hosea-project/figures/dp_missing.pdf", 8, 5)
plotdf = alldf %>% filter(df  %in% c("all", "cc", "all_0_5", "all_5_10", "all_10_30",
                                     "all_30_100"))
ggplot(data=plotdf, aes(x=tr, y=detection_prevalance, group=interaction(method, df), 
                        color=df, linetype=method)) + 
  geom_line() + scale_x_continuous(trans="log10")
dev.off()

pdf("R_code/hosea-project/figures/dp_vars.pdf", 8, 5)
plotdf = alldf %>% filter(!(df %in% c("all", "cc", "all_0_5", "all_5_10", "all_10_30",
                                     "all_30_100")))
ggplot(data=plotdf, aes(x=tr, y=detection_prevalance, group=interaction(method, df), 
                        color=df, linetype=method)) + 
  geom_line() + scale_x_continuous(trans="log10")
dev.off()



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

pdf("R_code/hosea-project/figures/aucs.pdf", 8, 5)
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
filepath = paste0("R_code/hosea-project/figures/roc_", tmethod, "_", tdf, ".pdf")
plotdf = rocs %>% filter(df == tdf & method == tmethod)
colnames(plotdf) = c("fpr", "tpr", "Threshold", "df", "method")
g = ggplot(data=plotdf, aes(x=fpr, y=tpr, color=Threshold)) + geom_line() +
  scale_color_gradientn(colors=rainbow(10)) + theme(aspect.ratio=1) +
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  geom_abline(intercept=0, slope=1, linetype="dotted")
# dev.off()
ggsave(filepath, g, width=8, height=7)

for(tdf in names(unweighted$metrics$calibration)){
  for(tmethod in c("resample", "unweighted")){
    auroc = aucs[aucs$df==tdf & aucs$method==tmethod, "auc"]
    # pdf(paste0("R_code/hosea-project/figures/roc/", tmethod, "_", tdf, ".pdf"), 6, 5)
    plotdf = rocs %>% filter(df == tdf & method == tmethod)
    g = ggplot(data=plotdf, aes(x=fpr, y=tpr, color=tr)) + geom_line() +
      scale_color_gradientn(colors=rainbow(20)) + theme(aspect.ratio=1) +
      ggtitle(paste0("Df: ", tdf, ", method: ", tmethod, ", AUROC=", round(auroc, 4)))
    # dev.off()
    ggsave(paste0("R_code/hosea-project/figures/roc/", tmethod, "_", tdf, ".pdf"), g,
           width=6, height=5)
  }
}

# calibration plot

all = lapply(resample$metrics$calibration_curves, function(df){
  df$method = "resample"
  return(df)
})
for(nm in names(all)){
  df = unweighted$metrics$calibration_curves[[nm]]
  df$method = "unweighted"
  all[[nm]]$tr = as.numeric(rownames(all[[nm]]))
  df$tr = as.numeric(rownames(df))
  all[[nm]] = rbind(all[[nm]], df)
}
alldf = data.frame()
for(nm in names(all)){
  df = all[[nm]]
  df$df = nm
  alldf = bind_rows(alldf, df)
}


tdf = "all"
tmethod = "resample"
log = T

for(tdf in names(unweighted$metrics$calibration)){
  for(tmethod in c("resample", "unweighted")){
    for(log in c(T, F)){
      filename = paste0("R_code/hosea-project/figures/calibration_curves/", tmethod, "_", 
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
        g = g + xlim(0., 1.) + ylim(0., 1.)
      }
      ggsave(filename, g, width=5, height=5)
    }
  }
}

# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrr}
# \hline
# & N & Ncases & propcases \\ 
# \hline
# all & 502849.00 & 2849.00 & 0.01 \\ 
# cc & 1571.00 & 12.00 & 0.01 \\ 
# all\_0\_5 & 83640.00 & 286.00 & 0.00 \\ 
# all\_5\_10 & 133193.00 & 369.00 & 0.00 \\ 
# all\_10\_30 & 163373.00 & 1122.00 & 0.01 \\ 
# all\_30\_100 & 121072.00 & 1060.00 & 0.01 \\ 
# lab\_vars\_0\_6 & 89372.00 & 304.00 & 0.00 \\ 
# lab\_vars\_6\_14 & 165220.00 & 593.00 & 0.00 \\ 
# lab\_vars\_14\_33 & 120373.00 & 843.00 & 0.01 \\ 
# lab\_vars\_33\_100 & 88946.00 & 444.00 & 0.00 \\ 
# lab\_vars\_100\_101 & 38938.00 & 665.00 & 0.02 \\ 
# charlson\_vars\_0\_1 & 495273.00 & 2537.00 & 0.01 \\ 
# charlson\_vars\_1\_101 & 7576.00 & 312.00 & 0.04 \\ 
# demo\_vars\_0\_1 & 416583.00 & 2244.00 & 0.01 \\ 
# demo\_vars\_1\_101 & 86266.00 & 605.00 & 0.01 \\ 
# other\_vars\_0\_101 & 502849.00 & 2849.00 & 0.01 \\ 
# smoke\_vars\_0\_50 & 284455.00 & 2044.00 & 0.01 \\ 
# smoke\_vars\_50\_101 & 218394.00 & 805.00 & 0.00 \\ 
# \hline
# \end{tabular}
# \end{table}