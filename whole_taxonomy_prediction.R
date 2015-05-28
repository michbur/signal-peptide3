library(dplyr)
library(hmeasure)
library(pbapply)


if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "D:/michal/doktorat/grant_data/signal_peptides/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/Qnap/Publikacje/signal_peptides/"

# all taxonomic groups -----------------------------------------
oc_pred <- read.csv2(file = paste0(pathway, "osoc_pred.csv"))

#-1 to remove Eukaryota
unique_ocs <- unique(unlist(strsplit(as.character(oc_pred[["oc"]]), "; ")))[-1]

full_oc_metrics <- pblapply(unique_ocs, function(single_oc) {
  dat <- oc_pred[grepl(single_oc, oc_pred[["oc"]]), c("signalHsmm", "signalP", "real")]
  both <- length(unique(dat[, "real"])) == 2
  list(metrics = if(both) {
    HMeasure(dat[["real"]], dat[, c("signalHsmm", "signalP")])[["metrics"]]
  } else {
    NULL
  }, count = nrow(dat),
  what = single_oc)
})

save(full_oc_metrics, unique_ocs, oc_counts, file = "/home/michal/Dropbox/signal-peptide2_data/whole_taxonomy_snap.RData")

oc_counts <- pbsapply(unique_ocs, function(single_oc)
  sum(grepl(single_oc, oc_pred[["oc"]])))
  
# #proper oc metrics
# oc_pmetrics <- full_oc_metrics[!sapply(full_oc_metrics, is.null)]
# 
# oc_AUC <- data.frame(names(oc_counts), oc_counts, t(sapply(oc_pmetrics, function(i) i[, "AUC"])))
# colnames(oc_AUC) <- c("OC", "Count", "signalHsmm", "signalP")
# oc_AUC <- cbind(oc_AUC, diff = oc_AUC[, "signalHsmm"] - oc_AUC[, "signalP"])
# oc_AUC <- oc_AUC[order(oc_AUC[, "diff"], decreasing = TRUE), ]
# write.csv2(oc_AUC[oc_AUC[, "Count"] > 50, ], file = paste0(pathway, "oc_AUC50.csv"))
# write.csv2(oc_AUC, file = paste0(pathway, "oc_AUC.csv"))