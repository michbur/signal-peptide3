#ngrams in cleavage sites.

source("start.R")

#sequences with cleavage sites
seqs_pos <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)

#positive sequences with at least 12 residues long SP
seqs_posl <- seqs_pos[sapply(seqs_pos, function(i) attr(i, "sig")[2]) - 12 > 1]

#cleavage site data
cs_data <- t(sapply(seqs_posl, function(i) {
  #cs is the posion of the FIRST amino acid after cleavage site
  cs <- attr(i, "sig")[2]
  #positions: -3, -1, +1
  i[(cs - 4):(cs + 3)]
}))

#did I choose appropriate positions?
# apply(cs_data, 2, table) %>% melt %>% ggplot(aes(x = Var2, y = value)) +
#   geom_bar(stat = "identity") + facet_wrap(~ Var1)

before_cs <- t(sapply(seqs_posl, function(i) {
  #cs is the posion of the FIRST amino acid after cleavage site
  cs <- attr(i, "sig")[2]
  #positions: -3, -1, +1
  i[(cs - 12):(cs - 5)]
}))

after_cs <- t(sapply(seqs_posl, function(i) {
  #cs is the posion of the FIRST amino acid after cleavage site
  cs <- attr(i, "sig")[2]
  #positions: -3, -1, +1
  i[(cs + 4):(cs + 11)]
}))

learning_dat <- rbind(cs_data, before_cs, after_cs) %>%
  apply(2, tolower) %>% degenerate(aaaggregation)
  
targets <- c(rep(1, nrow(cs_data)), rep(0, nrow(cs_data)*2))

#takes too long - start it on Thursday
#constructed_ngrams <- construct_ngrams(targets, learning_dat, as.character(1L:4), 4)
#save(constructed_ngrams, file = paste0(pathway, "constructed_ngrams.RData"))

distances1 <- 0
distances2 <- 0L:4
distances3 <- unlist(apply(expand.grid(0L:2, 0L:2), 1, list), recursive = FALSE)

bitrigrams <- count_multigrams(seq = learning_dat, 
                               ns = c(1, rep(2, length(distances2)), rep(3, length(distances3))), 
                               u = as.character(1L:4), 
                               ds = c(list(distances1), as.list(distances2), distances3),
                               pos = TRUE)

test_res <- test_features(targets, bitrigrams, adjust = NULL)
save(test_res, file = "cs_ngrams.RData")
