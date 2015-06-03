#cross-validation of signalHsmm with amino acid aggregation

library(signalHsmm)
library(seqinr)
library(pbapply)
library(cvTools)

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/Qnap/Publikacje/signalHsmm_cv/"

# data sets -----------------------------------------------

#UniProt release 2015_06

#/Qnap/Publikacje/signalHsmm_cv/signal_peptides.txt
#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) AND reviewed:yes
#2589 proteins in UniProt 

pos_seqs <- read_uniprot(paste0(pathway, "/data/signal_peptides.txt"), what = "signal", euk = TRUE)
pos_seqs <- pos_seqs[sapply(pos_seqs, length) > 80]
#2354 proteins after purification

#/Qnap/Publikacje/signalHsmm_cv/nonsignal_peptides.fasta
#taxonomy:"Eukaryota [2759]" NOT annotation:(type:signal evidence:any) AND reviewed:yes
#152272 proteins in UniProt 
neg_seqs <- read.fasta(paste0(pathway, "/data/nonsignal_peptides.fasta"), seqtype = "AA")
neg_seqs <- neg_seqs[sapply(neg_seqs, length) > 80]
neg_seqs <- neg_seqs[-which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B", "U"))))]
#138648 after filtering

# aggregation groups -------------------------------------
source("aggregation_groups.R")


# cross-validation -------------------------------------------
#number of repeats:
#qbinom(0.05, 155:165, prob=2589/138648)

fold_res <- pblapply(1L:50, function(dummy) {
  pos_ids <- cvFolds(length(pos_seqs), K = 5)
  cv_neg <- neg_seqs[sample(1L:length(neg_seqs), length(pos_seqs))]
  lapply(all_groups, function(agg_group) {
    lapply(1L:5, function(fold) {
      model_cv <- train_hsmm(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] != fold]], agg_group)
      test_dat <- c(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] == fold]],
                    cv_neg[pos_ids[[4]][,][pos_ids[[5]] == fold]])
      preds <- cbind(t(sapply(predict(model_cv, test_dat), function(single_pred)
        c(prob = single_pred[["sp_probability"]], cs_pred = single_pred[["sp_end"]]))),
        cs_real = sapply(test_dat, function(i) 
          ifelse(is.null(attr(i, "sig")[2]), NA, attr(i, "sig")[2])))
      preds
    })
  })
})

save(fold_res, file = paste0(pathway, "fold_res.RData"))

if(exists("fold_res")) {
  file.create("/home/michal/Dropbox/cv_success.txt")
} else {
  file.create("/home/michal/Dropbox/cv_fail.txt")
}