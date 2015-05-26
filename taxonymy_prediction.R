library(pbapply)
library(slam)

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "D:/michal/doktorat/grant_data/signal_peptides/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/Qnap/Publikacje/signal_peptides/"

read_signalp41 <- function(connection) {
  all_lines <- readLines(connection)[-c(1L:2)]
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(sp.probability = line[10] == "Y",
                      sp.start = ifelse(line[10] == "Y", 1, NA),
                      sp.end = ifelse(line[10] == "Y", as.numeric(line[5]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}


#

file_names <- c(paste0("nonsignal_peptides", 1L:20), "signal_peptides")

# all_oc <- unique(unlist(pblapply(file_names, function(single_file) {
#   load(paste0(pathway, single_file, ".RData"))
#   #just to remove Eukaryota from each oc
#   oc <- unique(unlist(sapply(strsplit(sapply(seqs, function(i) attr(i, "OC")), "; "), 
#                              function(i) i[-1])))
#   sub(".", "", oc, fixed = TRUE)
# })))
# 
# save(all_oc, file = paste0(pathway, "all_oc.RData"))

load(paste0(pathway, "all_oc.RData"))

# cool code below allocates too much of my memory
# all_oc_seqs <- pblapply(file_names, function(single_file) {
#   load(paste0(pathway, single_file, ".RData"))
#   write.csv2(vapply(seqs, function(single_seq) {
#     oc <- strsplit(attr(single_seq, "OC"), "; ")[[1]][-1]
#     oc[length(oc)] <- sub(".", "", oc[length(oc)], fixed = TRUE)
#     #c(attr(single_seq, "OS"))
#     all_oc %in% oc
#   }, rep(TRUE, length(all_oc))), 
#   file = paste0(pathway, single_file, "oc.csv"))
#   0
# })

signalHsmm_preds <- pblapply(file_names, function(single_file) {
  read.csv2(paste0(pathway, single_file, "signalHsmm.csv"))[, "sp.probability"]
})

signalP_preds <- pblapply(file_names, function(single_file) {
  read_signalp41(paste0(pathway, single_file, ".short_out"))[["sp.probability"]]
})

real_labels_len <- pbsapply(file_names, function(single_file) {
  nrow(read.csv2(paste0(pathway, single_file, "signalHsmm.csv")))
})

real_labels <- c(rep(FALSE, sum(real_labels_len[1L:20])), 
                 rep(TRUE, real_labels_len[21]))

os <- pblapply(file_names, function(single_file) {
  load(paste0(pathway, single_file, ".RData"))
  sapply(seqs, function(i) attr(i, "OS"))
})

os_pred <- data.frame(signalHsmm = unlist(signalHsmm_preds),
                      signalP = unlist(signalP_preds),
                      real = unlist(real_labels),
                      os = unlist(os))

library(dplyr)
library(hmeasure)



os_metrics <- lapply(levels(os_pred[["os"]]), function(i) {
  dat <- os_pred[os_pred[["os"]] == i, ]
  both <- length(unique(dat[, "real"])) == 2
  if(both) {
    HMeasure(dat[["real"]], dat[, c("signalHsmm", "signalP")])[["metrics"]]
  } else {
    NULL
  }
})
  
os_counts <- as.data.frame(table(os_pred[["os"]]))[!sapply(os_metrics, is.null), ]
#proper os metrics
os_pmetrics <- os_metrics[!sapply(os_metrics, is.null)]

os_AUC <- data.frame(os_counts, t(sapply(os_pmetrics, function(i) i[, "AUC"])))
colnames(os_AUC) <- c("OS", "Count", "signalHsmm", "signalP")
os_AUC <- cbind(os_AUC, diff = os_AUC[, "signalHsmm"] - os_AUC[, "signalP"])
os_AUC <- os_AUC[order(os_AUC[, "diff"], decreasing = TRUE), ]
write.csv2(os_AUC[os_AUC[, "Count"] > 50, ], file = "os_AUC50.csv")
write.csv2(os_AUC, file = "os_AUC.csv")

