library(dplyr)
library(hmeasure)
library(pbapply)

# all taxonomic groups -----------------------------------------
oc_pred <- read.csv2(file = paste0(pathway, "osoc_pred.csv"))

#-1 to remove Eukaryota
unique_ocs <- unique(unlist(strsplit(as.character(os_pred[["oc"]]), "; ")))[-1]

full_oc_metrics <- pblapply(unique_ocs, function(single_oc) {
  dat <- os_pred[grepl(single_oc, os_pred[["oc"]]), c("signalHsmm", "signalP", "real")]
  both <- length(unique(dat[, "real"])) == 2
  if(both) {
    HMeasure(dat[["real"]], dat[, c("signalHsmm", "signalP")])[["metrics"]]
  } else {
    NULL
  }
})

oc_counts <- sapply(unique_ocs, function(single_oc)
  length(grep(single_oc, oc_pred[["oc"]])))[!sapply(full_oc_metrics, is.null)]
  
#proper oc metrics
oc_pmetrics <- full_oc_metrics[!sapply(full_oc_metrics, is.null)]

oc_AUC <- data.frame(names(oc_counts), oc_counts, t(sapply(oc_pmetrics, function(i) i[, "AUC"])))
colnames(oc_AUC) <- c("OC", "Count", "signalHsmm", "signalP")
oc_AUC <- cbind(oc_AUC, diff = oc_AUC[, "signalHsmm"] - oc_AUC[, "signalP"])
oc_AUC <- oc_AUC[order(oc_AUC[, "diff"], decreasing = TRUE), ]
write.csv2(oc_AUC[oc_AUC[, "Count"] > 50, ], file = "oc_AUC50.csv")
write.csv2(oc_AUC, file = "oc_AUC.csv")