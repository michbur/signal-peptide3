#here I can plot easily unigrams against groups of amino acids. 
#I will use this code to create better groups for cleavage site prediction.

source("start.R")

#to predict cleavage sites, I need to analyze their properties.

#sequences with cleavage sites
seqs_pos <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)

#aminoacids in background
bg_data <- data.frame(t(sapply(seqs_pos, function(i) {
  i[2L:51]
}))) %>% na.omit #too short proteins will leave NAs in sequences

#cleavage site data
cs_data <- t(sapply(seqs_pos, function(i) {
  #cs is the posion of the FIRST amino acid after cleavage site
  cs <- attr(i, "sig")[2]
  #positions: -3, -1, +1
  i[(cs - 6):(cs+4)]
}))


##
# calculate background -----------------------------------
##

bg_counts <- bg_data %>% unlist %>% table %>% data.frame %>% select(Freq) %>% unlist

##
# calculate amino acids frequency in cs -----------------------------------
##

mcs_counts <- melt(data.frame(a = rownames(apply(cs_data, 2, table)), apply(cs_data, 2, table)))

levels(mcs_counts[["variable"]]) <- sapply(as.character(c(-5:5)), function(i) 
  ifelse(substr(i, 0, 1) == "-", i, paste0("+", i)))

#posible way how to reorder amino acids to make them appear in the order as they are in our aaagregation
#mconsensus[["a"]] <- factor(as.character(mconsensus[["a"]]), toupper(unlist(aaaggregation)))

aa_groups <- list(`1` = c("k", "r", "h"), 
                  `2` = c("v", "i", "l", "m", "f", "w", "c"), 
                  `3` = c("s", "t", "n", "q"), 
                  `4` = c("d", "e", "a", "p", "y", "g"))

#add group information to amino acids
suppressWarnings(mcs_counts <- cbind(group = sapply(names(unlist(aa_groups)), substr, 0, 1), mcs_counts))
#suppress expected warning 'row names were found from a short variable and have been discarded'

mcs_counts[["value"]] <- mcs_counts[["value"]]/bg_counts


ggplot(mcs_counts, aes(x = variable, y = value, fill = a)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete("Position") +
  scale_y_continuous("Frequency") +
  scale_fill_discrete("Amino acid") +
  geom_text(aes(label=a), position=position_dodge(width=0.9), vjust=-0.25) +
  facet_wrap(~ group)