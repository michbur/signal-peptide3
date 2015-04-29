source("start.R")

#to predict cleavage sites, I need to analyze their properties.

#sequences with cleavage sites
pos_seqs <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)

#index number of the cleavage site
id_cs <- sapply(pos_seqs, function(i) attr(i, "sig")[2])

#conensus after Heijne 1983 and Hiller 2004

consensus <- t(sapply(pos_seqs, function(i) {
  #cs is the posion of the FIRST amino acid after cleavage site
  cs <- attr(i, "sig")[2]
  #positions: -3, -1, +1
  i[(cs - 6):(cs+4)]
}))

consensus <- data.frame(consensus)
colnames(consensus) <- paste0("P", c(-5:5))

mconsensus <- melt(data.frame(a = rownames(apply(consensus, 2, table)), apply(consensus, 2, table)))

levels(mconsensus[["variable"]]) <- sapply(as.character(c(-5:5)), function(i) 
  ifelse(substr(i, 0, 1) == "-", i, paste0("+", i)))

#posible way how to reorder amino acids to make them appear in the order as they are in our aaagregation
#mconsensus[["a"]] <- factor(as.character(mconsensus[["a"]]), toupper(unlist(aaaggregation)))

aa_groups <- list(`1` = c("k", "r", "h"), 
                  `2` = c("v", "i", "l", "m", "f", "w", "c"), 
                  `3` = c("s", "t", "n", "q"), 
                  `4` = c("d", "e", "a", "p", "y", "g"))

mconsensus <- cbind(group = sapply(names(unlist(aa_groups)), substr, 0, 1), mconsensus)

ggplot(mconsensus, aes(x = variable, y = value, fill = a)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete("Position") +
  scale_y_continuous("Frequency") +
  scale_fill_discrete("Amino acid") +
  facet_wrap(~ group)





ggplot(mconsensus, aes(x = variable, y = value)) +
  geom_bar()