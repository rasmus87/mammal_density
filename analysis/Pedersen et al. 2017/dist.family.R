# Calculate mean pairwise distances

library(tidyverse)
library(ape)
library(tictoc)

mam <- read_csv("data/mam.9014.csv", col_types = cols())
mam <-  mam %>%
  # Order Cetacea (Whales s.l.)
  filter(!Family.1.2 %in% c("Balaenidae", "Balaenopteridae", "Delphinidae",
                            "Eschrichtiidae", "Iniidae", "Monodontidae",
                            "Neobalaenidae", "Phocoenidae", "Physeteridae",
                            "Ziphiidae")) %>%
  # Families in the clade Pinnipedia (Seal s.l.):
  # Odobenidae (walruses)
  # Otariidae (fur seals and sea lions)
  # Phocidae (true seals)
  filter(!Family.1.2 %in% c("Odobenidae", "Otariidae", "Phocidae")) %>%
  # Order Sirenia (Sea cows s.l.):
  filter(Order.1.2 != "Sirenia") %>%
  # Order Chiroptera (Bats)
  filter(Order.1.2 != "Chiroptera") %>%
  select(Family.1.2, Genus.1.2) %>% 
  filter(!duplicated(Genus.1.2))

phyl <- read.nexus("../large_datasets/Phylogenies/Complete_phylogeny.nex")
i <- 1
tree <- phyl[[i]]
tree$tip.label <- str_replace(tree$tip.label, "_.*", "")
tree$tip.label[which(tree$tip.label == "Pachyarmaterium")] <- "Pachyarmatherium"
matches <- match(tree$tip.label, mam$Genus.1.2)
tree$tip.label <- mam$Family.1.2[matches]
dup <- which(duplicated(tree$tip.label) | is.na(tree$tip.label))
tree <- drop.tip(tree, dup)
ord <- order(tree$tip.label)

rm(dist.all)
dist.all <- array(NA, c(length(ord), length(ord), 1000))
row.names(dist.all) <- tree$tip.label[ord]
colnames(dist.all) <- tree$tip.label[ord]
i <- 1
tic()
for(i in 1:1000) {
  tree <- phyl[[i]]
  tree$tip.label <- str_replace(tree$tip.label, "_.*", "")
  tree$tip.label[which(tree$tip.label == "Pachyarmaterium")] <- "Pachyarmatherium"
  matches <- match(tree$tip.label, mam$Genus.1.2)
  tree$tip.label <- mam$Family.1.2[matches]
  dup <- which(duplicated(tree$tip.label) | is.na(tree$tip.label))
  tree <- drop.tip(tree, dup)
  ord <- order(tree$tip.label)
  
  dist <- cophenetic.phylo(tree)
  dist <- dist[ord, ord]
  dist.all[, , i] <- dist
}
toc()


tic()
median.family.distance <- apply(dist.all, c(1,2), median)
toc()

saveRDS(median.family.distance, "builds/median.family.distance.rds")
