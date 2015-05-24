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

file_names <- c(paste0("nonsignal_peptides", 1L:9), "signal_peptides")

all_oc <- unique(unlist(pblapply(file_names, function(single_file) {
  load(paste0(pathway, single_file, ".RData"))
  #just to remove Eukaryota from each oc
  oc <- unique(unlist(sapply(strsplit(sapply(seqs, function(i) attr(i, "OC")), "; "), 
                             function(i) i[-1])))
  sub(".", "", oc, fixed = TRUE)
})))

save(all_oc, file = paste0(pathway, "all_oc.RData"))

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

real_labels <- c(rep(FALSE, sum(real_labels_len[1L:9])), rep(TRUE, real_labels_len[10]))

os <- unlist(pblapply(file_names, function(single_file) {
  load(paste0(pathway, single_file, ".RData"))
  sapply(seqs, function(i) attr(i, "OS"))
}))

data.frame(signalHsmm = signalHsmm_preds,
           signalP = signalP_preds,
           real = real_labels,
           os = os)
