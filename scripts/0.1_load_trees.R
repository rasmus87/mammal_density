# Loads the big forest and stores it as an R-object for faster handling
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(ape)

# Load PHYLACINE 1.2.1 posterior phylogeny distribution
forest <- read.nexus("../PHYLACINE_1.2/Data/Phylogenies/Complete_phylogeny.nex")

# Non terrestrials and humans
bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")
humans <- "Homo"

terrestrial <- mam %>% 
  filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
         !Family.1.2 %in% c(whale.families, seal.families),
         !Binomial.1.2 %in% marine.carnivores,
         !Genus.1.2 %in% humans) %>% 
  pull(Binomial.1.2)


# Reduce forest to the designated species list
species <- forest[[1]]$tip.label
drop.species <- species[!species %in% terrestrial]
forest <- lapply(forest, drop.tip, tip = drop.species)

# Store for fast loading
write_rds(forest, "builds/forest.rds")
