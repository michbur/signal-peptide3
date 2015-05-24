#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:manual) AND reviewed:yes
#read and filter data --------------------------------------------------

#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:manual) AND reviewed:yes
#positive data set, 26,737 proteins
#QNAP file signal_peptides.txt

#taxonomy:"Eukaryota [2759]" NOT annotation:(type:signal evidence:any) AND reviewed:yes
#negative data set, 152,185 proteins
#QNAP file nonsignal_peptides.txt

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/Qnap/Publikacje/signal_peptides/"

# all_lines <- readLines(paste0(pathway, "nonsignal_peptides.txt"))
# prot_ids <- grep("\\<ID   ", all_lines)
# border_ids <- round(seq(1L, 152185, length.out = 10), 0)
# #stupid idea, seq sets will be so uneven, but sigh
# sapply(1L:(length(border_ids) - 1), function(i)
#   writeLines(all_lines[prot_ids[border_ids[i]]:(prot_ids[border_ids[i + 1]] - 1)], 
#              con = paste0(pathway, "nonsignal_peptides", i, ".txt")))

library(signalHsmm)
library(seqinr)

source("signalp_benchmark_fun.R")

file_names <- c("nonsignal_peptides10", "nonsignal_peptides1", "nonsignal_peptides2", 
  "nonsignal_peptides3", "nonsignal_peptides4", "nonsignal_peptides5", 
  "nonsignal_peptides6", "nonsignal_peptides7", "nonsignal_peptides8", 
  "nonsignal_peptides9", "signal_peptides")

#sapply(file_names, benchmark_class)
benchmark_class("signal_peptides")
