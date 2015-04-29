library(signalHsmm)
library(biogram)
library(seqinr)
library(ggplot2)
library(reshape2)
library(dplyr)

# working directory ------------------------------

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"