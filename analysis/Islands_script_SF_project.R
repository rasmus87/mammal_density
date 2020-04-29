library(tidyverse)

dens <- read_csv("builds/imputed.density_3.csv")

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


dens.post <- read_csv("builds/333_densities_post.pred.10k.sample.csv")

# small.pop.sp <- "Prionailurus_bengalensis"
# island.Iriomote <- 290.68 # km2
# dens %>% filter(Binomial.1.2 == small.pop.sp) %>% select(log10.density.median, log10.lower.95hpd, log10.upper.95hpd) %>% 10^.*290.68

islands <- read_tsv("data/Islands.txt")

islands$Species_Name[-(1:2)][!islands$Species_Name[-(1:2)] %in% names(dens.post)]

i = 19 # Dingo
for (i in 3:132) {
  # home <- which(islands[i, 7:184] == TRUE)
  home <- (7:184)-6 # Fix calc for all islands
  j <- home[1]
  for(j in home) {
    island.size <- islands[[2, j + 6]] 
    species <- islands$Species_Name[i]
    if(species == "Lontra_felina") next
    pop.size <- island.size * 10^dens.post[species]
    result <- mean(pop.size >= 500)
    islands[i, j + 6] <- result
  }
}

write_tsv(islands, "output/Islands_500_likelihood_all_islands.txt")



islands <- read_tsv("data/Islands.txt")

islands$Species_Name[-(1:2)][!islands$Species_Name[-(1:2)] %in% names(dens.post)]

for(i in 3:132) {
  # home <- which(islands[i, 7:184] == TRUE)
  home <- (7:184)-6 # Fix calc for all islands
  j <- home[1]
  for(j in home) {
    island.size <- islands[[2, j + 6]] 
    species <- islands$Species_Name[i]
    if(species == "Lontra_felina") next
    pop.size <- island.size * 10^dens.post[species]
    result <- mean(pop.size >= 1000)
    islands[i, j + 6] <- result
  }
}

write_tsv(islands, "output/Islands_1000_likelihood_all_islands.txt")
