library(tidyverse)

dens <- read_csv("output/Table S4 Imputed density.csv")

islands <- read_tsv("data/Islands.txt")

anti_join(islands, dens, by = c("Species_Name" = "Binomial.1.2"))

islands$Species_Name %>% duplicated() %>% any

dens <- dens %>%
  filter(Binomial.1.2 %in% islands$Species_Name) %>% 
  transmute(Species_Name = Binomial.1.2,
            log10_Density = log10.density.mean,
            log10_Density_CI_lwr = log10.lower.95hpd,
            log10_Density_CI_upr = log10.upper.95hpd)

islands <- islands %>% left_join(dens)


write_tsv(islands, "output/Islands.txt")


# dens.post <- read_csv("builds/333_densities_post.pred.10k.sample.csv")
# We shuold switch to this in next update
dens.post <- read_csv("builds/3_densities_post.pred.csv")

# Build likelyhood for 500 animals
islands <- read_tsv("data/Islands.txt")

# Lonta felina not included in out dataset
islands$Species_Name[-(1:2)][!islands$Species_Name[-(1:2)] %in% names(dens.post)]

# Backup clean version
islands.bak <- islands

# Species numbers:
# 3:132, 1:2 are extra data
species.index <- 3:132
# Island numbers:
# 7:184, 1:6 are extra data
island.index <- 7:184


for (i in species.index) {
  for(j in island.index) {
    island.size <- islands[[2, j]] 
    species <- islands$Species_Name[i]
    if(species == "Lontra_felina") {
      islands[i, j] <- NA
      next
    } 
    pop.size <- island.size * 10^dens.post[species]
    result <- mean(pop.size >= 500)
    islands[i, j] <- result
  }
}
write_tsv(islands, "output/Islands_500_likelihood_all_islands.txt")


# Build likelyhood for 1000 animals
islands <- islands.bak

for (i in species.index) {
  for(j in island.index) {
    island.size <- islands[[2, j]] 
    species <- islands$Species_Name[i]
    if(species == "Lontra_felina") {
      islands[i, j] <- NA
      next
    } 
    pop.size <- island.size * 10^dens.post[species]
    result <- mean(pop.size >= 1000)
    islands[i, j] <- result
  }
}
write_tsv(islands, "output/Islands_1000_likelihood_all_islands.txt")
