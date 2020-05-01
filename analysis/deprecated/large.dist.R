# Calculate mean pairwise distances

library(tidyverse)
library(ape)
library(tictoc)
library(parallel)
cl <- makeCluster(8)
n.trees <- 1
test.size = 1000

mam <- read_csv("data/mam.9014.csv", col_types = cols())
mam <-  mam %>%
  select(Family.1.2, Binomial.1.2)

phyl <- read.nexus("../large_datasets/Phylogenies/Complete_phylogeny.nex")
i <- 1
tree <- phyl[[i]]
ord <- order(tree$tip.label)

rm(dist.all)
dist.all <- array(as.numeric(NA), c(length(ord), length(ord), n.trees))
row.names(dist.all) <- tree$tip.label[ord]
colnames(dist.all) <- tree$tip.label[ord]
i <- 1
tic()
for(i in 1:n.trees) {
  tree <- phyl[[i]]
  ord <- order(tree$tip.label)
  
  dist <- cophenetic.phylo(tree)
  dist <- dist[ord, ord]
  dist.all[, , i] <- dist
}
toc()

tic()
median.species.distance <- parApply(cl, dist.all[1:test.size, 1:test.size , ], c(1,2), median)
toc()

stopCluster(cl)

saveRDS(median.species.distance, "builds/median.species.distance.rds")
