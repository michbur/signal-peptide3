#counts of signal peptide for species and taxonomic groups

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "D:/michal/doktorat/grant_data/signal_peptides/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/Qnap/Publikacje/signal_peptides/"

osoc <- read.csv2(paste0(pathway, "osoc_pred.csv"))

sp_counts <-data.frame(table(osoc[osoc[["real"]], "os"]))

os_AUC50 <- read.csv2("os_AUC.csv")[, "OS"]

sp_counts <- sp_counts[sp_counts[, 1] %in% os_AUC50, ]

write.csv2(sp_counts, file = "os_AUC50_counts.csv")
