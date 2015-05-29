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

unique_ocs <- sub(".", "", unique_ocs, fixed = TRUE)

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

#load(file = paste0(pathway, "whole_taxonomy_snap.RData"))

informative <- !sapply(full_oc_metrics, function(i) is.null(i[["metrics"]]))

oc_AUC <- data.frame(unique_ocs[informative], t(sapply(full_oc_metrics[informative], function(i)
  c(i[["count"]], i[["metrics"]][, "AUC"])
)))

colnames(oc_AUC) <- c("OC", "Count", "signalHsmm", "signalP")
oc_AUC <- cbind(oc_AUC, diff = oc_AUC[, "signalHsmm"] - oc_AUC[, "signalP"])
oc_AUC <- oc_AUC[order(oc_AUC[, "diff"], decreasing = TRUE), ]

#no necessary now, I'm filtering unique_ocs
#levels(oc_AUC[["OC"]]) <- sub(".", "", levels(oc_AUC[["OC"]]), fixed = TRUE)
#oc_AUC <- oc_AUC[!duplicated(oc_AUC), ]


write.csv2(oc_AUC[oc_AUC[, "Count"] > 50, ], file = paste0(pathway, "oc_AUC50.csv"))
write.csv2(oc_AUC, file = paste0(pathway, "oc_AUC.csv"))


