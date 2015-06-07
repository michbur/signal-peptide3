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

read_predsi <- function(connection) {
  dat <- read.table(connection, sep = "\t")
  data.frame(sp.probability = dat[, 4] == "Y",
             sig.start = ifelse(dat[, 4] == "Y", 1, NA),
             sig.end = ifelse(dat[, 4] == "Y", as.numeric(dat[, 3]), NA),
             row.names = dat[, 1])
}

read_phobius <- function(connection) {
  all_lines <- readLines(connection)
  all_lines <- all_lines[-1]
  splited <- strsplit(all_lines, " ")
  #remove "" characters
  purged <- t(sapply(splited, function(i) i[i != ""]))
  cl_sites <- sapply(which(purged[, 3] == "Y"), function(i)
    as.numeric(strsplit(strsplit(purged[i,4], "/")[[1]][1], "c")[[1]][[2]]))
  res <- data.frame(sp.probability = purged[, 3] == "Y",
                    sig.start = ifelse(purged[, 3] == "Y", 1, NA),
                    sig.end = rep(NA, nrow(purged)), 
                    row.names = purged[, 1])
  res[purged[, 3] == "Y", "sig.end"] <- cl_sites
  res
}

read_philius <- function(connection) {
  require(XML)
  all_dat <- xmlToList(xmlTreeParse(connection, asTree = TRUE))
  seq_dat_id <- 1L:(length(all_dat)/2)*2
  #data for table
  table_dat <- sapply(seq_dat_id, function(i) 
    unlist(all_dat[i][[1]][[1]][c(24, 22)]))
  cleaved <- sapply(table_dat, function(i)
    !(is.null(i[1]) || is.na(i[1])))
  res <- data.frame(sp.probability = cleaved,
                    sig.start = ifelse(cleaved, 1, NA),
                    sig.end = rep(NA, length(seq_dat_id)),
                    row.names = unlist(all_dat[1L:(length(all_dat)/2)*2 - 1]))
  res[cleaved, "sig.end"] <- as.numeric(sapply(table_dat[cleaved], function(i) i[2])) - 1
  res
}

library(signalHsmm)
library(seqinr)
library(hmeasure)






#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500601 TO 20100101] AND reviewed:yes
pos_seqs2010 <- read_uniprot("sp1986_2010.txt", what = "signal", euk = TRUE)
pos_seqs2010 <- pos_seqs2010[sapply(pos_seqs2010, length) > 80]
#2169 seqs

#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[20110101 TO 20150601] AND reviewed:yes
pos_seqs2015 <- read.fasta("sp2010_2015.fasta", seqtype = "AA")
#218 seqs

#taxonomy:"Eukaryota [2759]" NOT annotation:(type:signal evidence:any) created:[20100101 TO 20150601] AND reviewed:yes
neg_seqs2015 <- read.fasta("nsp2010_2015.fasta", seqtype = "AA")
#18610 seqs

#benchmark_dat <- c(pos_seqs2015, sample(neg_seqs2015, length(pos_seqs2015)))
#write.fasta(benchmark_dat, names(benchmark_dat), file = "benchmark_data.fasta")

real_labels <- c(rep(1, length(pos_seqs2015)), rep(0, length(pos_seqs2015)))

bench_metrics <- HMeasure(real_labels, 
                          data.frame(signalPnotm = read_signalp41("sp4_notm.txt")[["sp.probability"]], 
                                     signalPtm = read_signalp41("sp4_tm.txt")[["sp.probability"]], 
                                     predsi = read_predsi("predisi.txt")[["sp.probability"]],
                                     phobius = read_phobius("phobius.txt")[["sp.probability"]],
                                     philius = read_philius("philius.xml")[["sp.probability"]],
                                     signalHsmm2010 = pred2df(predict(train_hsmm(read_uniprot("sp1950_2010.txt", euk = TRUE, what = "signal"), aa_group = all_groups[[22]]),
                                                                      benchmark_dat))[["sp.probability"]],
                                     signalHsmm1989 = pred2df(predict(train_hsmm(read_uniprot("sp1950_1989.txt", euk = TRUE, what = "signal"), aa_group = all_groups[[22]]),
                                                                      benchmark_dat))[["sp.probability"]]),
                          threshold = c(rep(0.5, 5), 
                                        0.05, 0.05))[["metrics"]]


TP <- as.numeric(bench_metrics[, "TP"])
TN <- bench_metrics[, "TN"]
FP <- bench_metrics[, "FP"]
FN <- bench_metrics[, "FN"]
bench_metrics <- cbind(bench_metrics, 
                       MCC = unname((TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))
bench_metrics[, c("AUC", "MCC")]