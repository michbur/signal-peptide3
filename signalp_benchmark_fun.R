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


benchmark_class <- function(filename) {
  seqs <- read_uniprot(paste0(pathway, filename, ".txt"), euk = TRUE, what = NULL)
  
  lengths <- sapply(seqs, length)
  
  seqs <- seqs[lengths > 80]
  
  save(seqs, file = paste0(pathway, filename, ".RData"))
  
  write.fasta(seqs, names(seqs), paste0(pathway, filename, ".fasta"))
  
  system(paste0("/home/michal/signalp-4.1/signalp -t euk -f short ", 
                paste0(pathway, filename, ".fasta"), " > ",
                paste0(pathway, filename, ".short_out")))
  
  write.csv2(pred2df(run_signalHsmm(seqs)), file = paste0(pathway, filename, "signalHsmm.csv"))
  
  file.create(paste0("/home/michal/Dropbox/signal-peptide2_data/", filename, ".sucess"))
  
  0
}
