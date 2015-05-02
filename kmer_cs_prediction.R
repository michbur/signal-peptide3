#adding n-gram to hsmm to predict cleavage site more accurately

source("start.R")

signal.hsmm2010 <- train_hsmm(read_uniprot(paste0(pathway, "pub_pos_train.txt"), euk = TRUE),
                              aaaggregation)


train_hsmm(train_data, aa_group, max_length = 32)

add_k_mer_state(c("111", "112"), signal.hsmm2010[["pipar"]], signal.hsmm2010[["tpmpar"]],
                signal.hsmm2010[["od"]], signal.hsmm2010[["params"]], 3, 4, 0.25)

