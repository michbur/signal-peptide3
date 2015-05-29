library(ggvis)
library(dplyr)

os50 <- read.csv2("os_AUC50.csv")[, -1]

levels(os50[["OS"]]) <- sub("\\.[^\\.]*$", "", levels(os50[["OS"]]))

#how odd
os50 %>% ggvis(~signalHsmm, ~signalP, size = ~Count) %>% layer_points()  

os50 <- cbind(cCount = cut(os50[["Count"]], c(50, 100, 150, 300, 500, 1000, 70000)), os50)

better_hsmm <- filter(os50, signalP < 0.8, diff > 0) %>% droplevels
levels(better_hsmm[["OS"]]) <- c("Aplysia californica", 
                                 "Brassica oleracea", 
                                 "Citrus sinensis", 
                                 "Coffea arabica", 
                                 "Crithidia fasciculata", 
                                 "Cucumis sativus", 
                                 "Culex quinquefasciatus", 
                                 "Daucus carota", 
                                 "Leishmania donovani", 
                                 "Leishmania major", 
                                 "Phaeodactylum tricornutum", 
                                 "Plasmodium berghei (strain Anka)", 
                                 "Plasmodium falciparum", 
                                 "Plasmodium falciparum (isolate 3D7)", 
                                 "Plasmodium falciparum (isolate K1 / Thailand)", 
                                 "Plasmodium vivax (strain Salvador I)", 
                                 "Theileria annulata", 
                                 "Trichoplax adhaerens", 
                                 "Trypanosoma cruzi", 
                                 "Vigna radiata var. radiata")

library(ggplot2)
os50 %>% filter(diff > 0) %>% ggplot(aes(x = signalHsmm, y = signalP, colour = cCount)) +
  geom_point(size = 5) +
  scale_color_discrete(h = c(60, 210) + 15, direction = -1) +
  geom_text(aes(x = signalHsmm, y = signalP, label = OS), better_hsmm, hjust=0.3, vjust=2, size = 4)


