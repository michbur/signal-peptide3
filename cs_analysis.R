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


#roi - region of the interest
#important - use full d argument, ngrams_freqs doesn't recycle d as other function do
ngrams_freqs <- function(roi, background, n = 1, d = 0) {
  
  ##
  # calculate background -----------------------------------
  ##
  
  bg_table <- background %>% as.matrix %>% seq2ngrams(n, a()[-1], d, pos = FALSE) %>%
    unlist %>% table %>% data.frame 
  
  bg_counts <- bg_table[["Freq"]]
  names(bg_counts) <- bg_table[["."]]

  
  ##
  # calculate amino acids frequency in cs -----------------------------------
  ##
  
  roi_ngrams <- roi %>% as.matrix %>% seq2ngrams(n, a()[-1], d, pos = FALSE)
  
  roi_counts <- cbind(bg_table[["."]], do.call(cbind, lapply(1L:ncol(roi_ngrams), function(single_position) {
    factor(roi_ngrams[, single_position], levels = bg_table[["."]]) %>% 
      table %>% data.frame %>% select(Freq)
  })))
  #solve it later by smart select call with .dots
  
  colnames(roi_counts) <- c("ngrams", sapply(as.character(c(-5:(5 - n - sum(d) + 1))), function(i) 
    ifelse(substr(i, 0, 1) == "-", i, paste0("+", i))))
    
  mroi_counts <- melt(roi_counts)
    
  mroi_counts[["freq"]] <- mroi_counts[["value"]]/bg_counts
  mroi_counts
}

freqs2 <- ngrams_freqs(cs_data, bg_data, 3, d = c(1, 0))

library(entropy)
#calculate entropy to assess which sites ae the most informative
group_by(freqs2, variable) %>% summarise(entropy = entropy.empirical(value))

group_by(freqs2, variable) %>% filter(freq > quantile(freq, 0.99))


ggplot(freqs2, aes(x = variable, y = value, fill = ngrams)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete("Position") +
  scale_y_continuous("Frequency") +
  scale_fill_discrete("Amino acid") +
  geom_text(aes(label=ngrams), position=position_dodge(width=0.9), vjust=-0.25) 

