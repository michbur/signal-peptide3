#creation of aggregation groups for cv
load(paste0(pathway, "data/aa_nprop.RData"))

traits <- list(size = c(30, 36, 54),
               polarity = 202,
               pi = 203,
               hydroph = c(1, 26, 57),
               alpha = c(12, 23, 232))

grouping_properties <- t(aa_nprop[unlist(traits), ])

#traits without polarity
all_traits_combn_nonpol <- cbind(expand.grid(traits[["size"]], traits[["hydroph"]], 
                                             traits[["alpha"]]), traits[["pi"]])

colnames(all_traits_combn_nonpol) <- c("size", "hydroph", "alpha", "pi")


all_groups_nonpol <- lapply(1L:nrow(all_traits_combn_nonpol), function(single_trait_combn) {
  cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn_nonpol[single_trait_combn, ]), ])))
  gr <- cutree(cl, k = 4)
  names(gr) <- tolower(names(gr))
  agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
  names(agg_gr) <- 1L:length(agg_gr)
  agg_gr
})


all_groups <- all_groups_nonpol
