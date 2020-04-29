# Calculate mean pairwise distances

library(tidyverse)
library(ape)
library(tictoc)
n.trees <- 10

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())
mam <-  mam %>%
  select(Family.1.2, Binomial.1.2)

phyl <- read.nexus("../PHYLACINE_1.1/Data/Phylogenies/Complete_phylogeny.nex")
i <- 1
tree <- phyl[[i]]
ord <- order(tree$tip.label)
median.species.distance <- cophenetic.phylo(tree)
median.species.distance <- median.species.distance[ord, ord]
median.species.distance[] <- NA
storage.mode(median.species.distance) <- "integer"

rm(dist.all)
dist.all <- matrix(as.integer(NA), length(ord)^2, n.trees)
i <- 1
tic()
for(i in 1:n.trees) {
  tree <- phyl[[i]]
  ord <- order(tree$tip.label)
  
  dist <- cophenetic.phylo(tree)
  dist <- as.integer(dist[ord, ord] * 1000)
  dist.all[, i] <- dist
}
toc()

write_csv(as.data.frame(dist.all), "builds/dist.all.csv", col_names = FALSE)
rm(dist.all)
gc()

timestamp()
tic()
size <- length(ord)^2
chunk.size <- 10000
connection <- file("builds/dist.all.csv", "r")
result <- rep(as.integer(NA), size)
for(i in 1:ceiling(size/chunk.size)) {
  print(paste0(i, "/", size/chunk.size))
  data <- readLines(connection, chunk.size)
  data <- read.csv(textConnection(data), header = FALSE)
  result[(1:chunk.size)+((i-1)*chunk.size)] <- apply(data, 1, median)
}
result <- result[1:size]
close(connection)
toc()

median.species.distance[] <- result

saveRDS(median.species.distance, "builds/median.species.distance.rds")

median.species.distance <- as_tibble(median.species.distance, rownames = "Species")

write_csv(median.species.distance, "builds/median.species.distance.csv")
