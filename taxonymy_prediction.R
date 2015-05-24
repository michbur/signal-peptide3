library(pbapply)
library(slam)

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "D:/michal/doktorat/grant_data/signal_peptides/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/Qnap/Publikacje/signal_peptides/"

#os <- sapply(seqs, function(i) attr(i, "OS"))

file_names <- c(paste0("nonsignal_peptides", 1L:9), "signal_peptides")

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

all_oc_seqs <- pblapply(file_names, function(single_file) {
  load(paste0(pathway, single_file, ".RData"))
  write.csv2(vapply(seqs, function(single_seq) {
    oc <- strsplit(attr(single_seq, "OC"), "; ")[[1]][-1]
    oc[length(oc)] <- sub(".", "", oc[length(oc)], fixed = TRUE)
    #c(attr(single_seq, "OS"))
    all_oc %in% oc
  }, rep(TRUE, length(all_oc))), 
  file = paste0(pathway, single_file, "oc.csv"))
  0
})

