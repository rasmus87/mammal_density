library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)

## Set options:
# Set parralell cluster size
cluster.size <- 2
cluster.size <- 6
#cluster.size <- 20
# How many trees do you want to run this for? 2-1000?
n.trees <- 2
n.trees <- 6
# n.trees <- 6*10
#n.trees <- 1000

pantheria <- read_tsv("../metabolic_rate/data/PanTHERIA_1-0_WR05_Aug2008.txt", col_types = cols())
names(pantheria) <- make.names(names(pantheria))
pantheria <- pantheria %>% 
  filter(!is.na(X21.1_PopulationDensity_n.km2)) %>% 
  transmute(Binomial.Pantheria = paste(MSW05_Genus, MSW05_Species, sep = "_"),
            density = X21.1_PopulationDensity_n.km2,
            source = "PanTHERIA.2008")
# Load taxonomy solver
syn <- read_csv("../PHYLACINE_1.1/Data/Taxonomy/Synonymy_table_with_unaccepted_species.csv", col_types = cols(), guess_max = 5000)
syn <- syn %>% 
  transmute(Binomial.1.2,
            Binomial.EltonTraits.1.0 = paste(EltonTraits.1.0.Genus, EltonTraits.1.0.Species, sep = "_"),
            Binomial.1.0 = paste(Genus.1.0, Species.1.0, sep = "_"),
            Binomial.1.1 = paste(Genus.1.1, Species.1.0, sep = "_")) %>% 
  filter(!str_detect(Binomial.1.2, "000"))

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

syn <- syn %>% 
  gather(System, Binomial, -Binomial.1.2) %>% 
  filter(Binomial.1.2 != Binomial) %>% 
  filter(!duplicated(Binomial)) %>% 
  bind_rows(mam %>% transmute(Binomial.1.2,
                              System = "Binomial.1.2",
                              Binomial = Binomial.1.2))
anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial"))

# Cercopithecus_pogonias
# Possibly: Cercopithecus denti
# http://www.iucnredlist.org/details/136885/0

# Cebus_olivaceus is split into two species:
# Cebus brunneus and Cebus kaapori - pick a random

extra <- data_frame(Binomial.1.2 = c("Cercopithecus_denti", "Cebus_brunneus"), 
                    System = "Personal.judgement", 
                    Binomial = c("Cercopithecus_pogonias", "Cebus_olivaceus"))
syn <- syn %>% bind_rows(extra)
stopifnot(!anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial")) %>% nrow)

pantheria <- left_join(pantheria, syn %>% select(-System), by = c("Binomial.Pantheria" = "Binomial"))

pantheria <- pantheria %>% left_join(mam)

dup <- pantheria %>% filter(duplicated(Binomial.1.2)) %>% pull(Binomial.1.2)
pantheria %>% filter(Binomial.1.2 %in% dup) %>% select(Binomial.1.2, Binomial.Pantheria)

# Remove duplicated species and keep nominal
pantheria <- pantheria %>% 
  filter(!Binomial.Pantheria %in% c("Damaliscus_korrigum",
                                    "Felis_catus",
                                    "Alcelaphus_caama",
                                    "Alcelaphus_lichtensteinii"))
stopifnot(!anyDuplicated(pantheria$Binomial.1.2))

pantheria <- pantheria %>% 
  mutate(log10density = log10(density),
         log10BM = log10(Mass.g)) %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10density, log10BM, source)

terrestrial <- mam %>% filter(Terrestrial == 1) %>% pull(Binomial.1.2)

# Linear model:
pantheria <- pantheria %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "density")
pantheria <- as.data.frame(pantheria)

mam <- mam %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  mutate(log10BM = log10(Mass.g), log10density = NA)
mam <- mam %>% mutate(dataset = "mam")
mam <- mam %>% select(names(pantheria))
mam <- as.data.frame(mam)

n.mam <- nrow(mam)

df <- rbind(mam, pantheria)
df <- as.data.frame(df)

forest <- readRDS("../metabolic_rate/builds/forest.rds")
species <- forest[[1]]$tip.label
drop.species <- species[!species %in% terrestrial]
forest <- lapply(forest, drop.tip, tip = drop.species)

prior <- list(G = list(G1 = list(V = 1, nu = 0.02)), 
              R = list(V = 1, nu = 0.02))
thin <- 75
burnin <- thin * 10
nitt <- 333 * thin + burnin
i = 1
mcmc.regression <- function(i) {
  tree <- forest[[i]]
  inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
  chain.1 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = pantheria, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE)
  chain.2 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = pantheria, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE)
  chain.3 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = pantheria, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE)
  if(i == 1) {
    saveRDS(chain.1, paste0("builds/mcmcglmms/tree", i, ".chain1.rds"), compress = FALSE)
    saveRDS(chain.2, paste0("builds/mcmcglmms/tree", i, ".chain2.rds"), compress = FALSE)
    saveRDS(chain.3, paste0("builds/mcmcglmms/tree", i, ".chain3.rds"), compress = FALSE)
  }
  
  gc()
  pred1 <- predict(chain.1, df, marginal = NULL, interval = "confidence")[1:n.mam, ]
  pred2 <- predict(chain.2, df, marginal = NULL, interval = "confidence")[1:n.mam, ]
  pred3 <- predict(chain.3, df, marginal = NULL, interval = "confidence")[1:n.mam, ]
  
  pred1 <- as.data.frame(pred1)
  pred2 <- as.data.frame(pred2)
  pred3 <- as.data.frame(pred3)

  names(pred1) <- c("fit", "lwrCI", "uprCI")
  names(pred2) <- c("fit", "lwrCI", "uprCI")
  names(pred3) <- c("fit", "lwrCI", "uprCI")
  
  pred1["Binomial.1.2"] <- mam$Binomial.1.2
  pred2["Binomial.1.2"] <- mam$Binomial.1.2
  pred3["Binomial.1.2"] <- mam$Binomial.1.2

  pred <- rbind(pred1, pred2, pred3)
  pred["tree"] <- i
  
  solution <- rbind(chain.1$Sol[, 1:2],
                    chain.2$Sol[, 1:2],
                    chain.3$Sol[, 1:2])
  solution <- as.data.frame(solution)
  random.effects <- 3:ncol(chain.1$Sol)
  solution["random.effect"] <- c(rowMeans(chain.1$Sol[, random.effects]),
                                 rowMeans(chain.1$Sol[, random.effects]),
                                 rowMeans(chain.1$Sol[, random.effects]))
  solution["tree"] <- i
  solution["chain"] <- as.numeric(gl(3, 333))
  
  return(list(pred, solution))
}

comb <- function(p, q) {
  list(rbind(p[[1]], q[[1]]),
       rbind(p[[2]], q[[2]]))
}

cl <- makeCluster(cluster.size)
registerDoSNOW(cl)
timestamp()
tic()
pb <- txtProgressBar(max = n.trees, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
imputed <- foreach(i = 1:n.trees,
                   .packages = c('MCMCglmm'),
                   .inorder = FALSE,
                   .options.snow = opts,
                   .combine = comb,
                   .multicombine = FALSE) %dopar% mcmc.regression(i)
toc()
stopCluster(cl)
gc()

write_csv(imputed[[1]], "builds/test_imputed.densities_all.samples.csv")
write_csv(as_data_frame(imputed[[2]]), "builds/test_densities_fit.solution.csv")
# 
# write_csv(imputed[[1]], "builds/imputed.metabolic.rate_all.samples.csv")
# write_csv(as_data_frame(imputed[[2]]), "builds/metabolic.rate_fit.solution.csv")
