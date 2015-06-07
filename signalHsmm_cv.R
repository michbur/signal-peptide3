#cross-validation of signalHsmm with amino acid aggregation

library(signalHsmm)
library(seqinr)
library(pbapply)
library(cvTools)
library(hmeasure)
library(dplyr)
library(reshape2)

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "D:/michal/doktorat/grant_data/signalHsmm_cv/"

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

load(paste0(pathway, "fold_res.RData"))


total_AUC <- data.frame(repetition = sort(rep(1L:50, length(all_groups)*5)), 
                        do.call(rbind, lapply(fold_res, function(single_repeat) {
                          data.frame(agg = sort(rep(1L:length(single_repeat), 5)),
                                     do.call(rbind, lapply(single_repeat, function(single_group)
                                       data.frame(fold = 1L:5, t(sapply(single_group, function(single_fold) {
                                         unlist(HMeasure(as.numeric(!is.na(single_fold[, "cs_real"])), 
                                                  single_fold[, "prob"])[["metrics"]][c("H", "AUC")])
                                       }))
                                       )
                                     )))
                        })))

#total_AUC[, "H"] <- unlist(total_AUC[, "H"])
#total_AUC[, "AUC"] <- unlist(total_AUC[, "AUC"])
total_AUC %>% 
  group_by(agg) %>% 
  summarize(m_AUC = mean(AUC), sd_AUC = sd(AUC)) %>%
  arrange(desc(m_AUC))
#Source: local data frame [27 x 3]
#
#agg     m_AUC      sd_AUC
#1   22 0.9600330 0.006160980
#2    7 0.9593754 0.006076867

thresh_vals <- c(0.005, 0.01, 1:12/20)
mcc_optimization <- lapply(thresh_vals, function(thres)
  sapply(fold_res, function(i) rowMeans(sapply(i[[22]], function(single_fold) {
    res <- unlist(HMeasure(as.numeric(!is.na(single_fold[, "cs_real"])), 
                           single_fold[, "prob"],
                           threshold = thres)[["metrics"]][c("TP", "FP", "TN", "FN", "Sens", "Spec")])
    TP <- as.numeric(res["TP"])
    TN <- res["TN"]
    FP <- res["FP"]
    FN <- res["FN"]
    c(mcc = unname((TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))),
      sens = res[["Sens"]], spec = res[["Spec"]])
  }))))


mmcc <- t(sapply(mcc_optimization, rowMeans)) %>% data.frame %>% mutate(sensspec = (sens + spec)/2) %>% melt
mmcc <- cbind(rep(thresh_vals, 4), mmcc)
colnames(mmcc) <- c("thresh", "perf", "value")



library(ggplot2)
ggplot(mmcc, aes(x = thresh, y = value, colour = perf)) +
  geom_line()


#final values of MCC
final_metrics <- sapply(fold_res, function(i) rowMeans(sapply(i[[22]], function(single_fold) {
    res <- unlist(HMeasure(as.numeric(!is.na(single_fold[, "cs_real"])), 
                           single_fold[, "prob"],
                           threshold = 0.01)[["metrics"]])
    TP <- as.numeric(res["TP"])
    TN <- res["TN"]
    FP <- res["FP"]
    FN <- res["FN"]
    c(res, MCC = unname((TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))
  })))
  

nice_finals <- format(t(apply(final_metrics, 1, function(i)
  c(mean(i), sd(i), min(i), max(i)))), trim = TRUE, digits = 1,  small.interval = 0, scientific = FALSE, justify = "none")

nice_finals <- nice_finals[c(3, 1, 23, 2, 5, 11L:14, 19L:22), ]
rownames(nice_finals) <- c("Area Under the Curve",
                           "H-measure",
                           "Matthew's Correlation Coefficient",
                           "Gini index",
                           "Kolmogorov-Smirnoff statistic",
                           "Sensitivity",
                           "Specificity",
                           "Precision",
                           "Recall",
                           "TP",
                           "FP",
                           "TN",
                           "FN")

write.table(nice_finals, file = "clipboard", sep = "\t")
