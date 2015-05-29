#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:manual) AND reviewed:yes
#read and filter data --------------------------------------------------

#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:manual) AND reviewed:yes
#positive data set, 26,737 proteins
#QNAP file signal_peptides.txt

#taxonomy:"Eukaryota [2759]" NOT annotation:(type:signal evidence:any) AND reviewed:yes
#negative data set, 152,185 proteins
#QNAP file nonsignal_peptides.txt


if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "D:/michal/doktorat/grant_data/signal_peptides/"


if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/Qnap/Publikacje/signal_peptides/"

all_lines <- readLines(paste0(pathway, "nonsignal_peptides.txt"))
prot_ids <- grep("\\<ID   ", all_lines)


library(cvTools)
n_splits <- 20
splits <- c(cvFolds(length(prot_ids), K = n_splits, type = "consecutive")[["which"]], n_splits + 1)
prot_ids <- c(prot_ids, length(all_lines))

sapply(1L:n_splits, function(i)
  writeLines(all_lines[(prot_ids[splits == i][1]):((prot_ids[splits == (i + 1)][1]) - 1)], 
             con = paste0(pathway, "nonsignal_peptides", i, ".txt")))

library(signalHsmm)
library(seqinr)

source("signalp_benchmark_fun.R")

file_names <- paste0("nonsignal_peptides", 1L:n_splits)

sapply(file_names, benchmark_class)
#benchmark_class("signal_peptides")
