library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)

# 36 hours runtime for 1000 trees on 20 cores
# ~ 80 GB RAM required
# ~ 80 GB disk space

## Set options:
# Set parralell cluster size
cluster.size <- 2
# cluster.size <- 6
cluster.size <- 20
# How many trees do you want to run this for? 2-1000?
n.trees <- 2
# n.trees <- 6
# n.trees <- 6*10
n.trees <- 1000

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

bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")

terrestrial <- mam %>% filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
                              !Family.1.2 %in% c(whale.families, seal.families),
                              !Binomial.1.2 %in% marine.carnivores) %>% pull(Binomial.1.2)

pantheria <- pantheria %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "density")
write_csv(pantheria, "builds/imputation_dataset.csv")
pantheria <- as.data.frame(pantheria)

# Linear model:
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
samples <- 33
thin <- 75
burnin <- thin * 10
nitt <- samples * thin + burnin
i = 1
mcmc.regression <- function(i) {
  tree <- forest[[i]]
  inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
  chain.1 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = pantheria, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.2 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = pantheria, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.3 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = pantheria, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  if(i == 1) {
    saveRDS(chain.1, paste0("builds/mcmcglmms/tree", i, ".chain1.rds"), compress = FALSE)
    saveRDS(chain.2, paste0("builds/mcmcglmms/tree", i, ".chain2.rds"), compress = FALSE)
    saveRDS(chain.3, paste0("builds/mcmcglmms/tree", i, ".chain3.rds"), compress = FALSE)
  }
  
  gc()
  pred1 <- MCMC.predict.v2(chain.1, df)
  pred2 <- MCMC.predict.v2(chain.2, df)
  pred3 <- MCMC.predict.v2(chain.3, df)

  pred <- rbind(pred1[[2]], pred2[[2]], pred3[[2]])
  pred["tree"] <- i
  
  post.pred <- rbind(pred1[[1]], pred2[[1]], pred3[[1]])
  
  solution <- rbind(chain.1$Sol[, 1:2],
                    chain.2$Sol[, 1:2],
                    chain.3$Sol[, 1:2])
  solution <- as.data.frame(solution)
  random.effects <- 3:ncol(chain.1$Sol)
  solution["random.effect"] <- c(rowMeans(chain.1$Sol[, random.effects]),
                                 rowMeans(chain.2$Sol[, random.effects]),
                                 rowMeans(chain.3$Sol[, random.effects]))
  solution["tree"] <- i
  solution["chain"] <- as.numeric(gl(3, samples))
  
  return(list(pred, solution, post.pred))
}

MCMC.predict.v2 <- function(object, newdata) {
  object2 <- MCMCglmm(fixed=object$Fixed$formula, 
                      random=object$Random$formula, 
                      rcov=object$Residual$formula, 
                      family=object$Residual$original.family,
                      data=newdata, 
                      nitt=1, 
                      thin=1,
                      burnin=0, 
                      ginverse=object$ginverse, 
                      verbose=FALSE, 
                      pr=TRUE)
  
  W <- cbind(object2$X, object2$Z)
  post.pred <- t(apply(object$Sol, 1, function(x){(W %*% x)@x}))[, 1:n.mam]
  
  se <- apply(post.pred, 2, sd)
  pred <- colMeans(post.pred)
  pred <- cbind(mam$Binomial.1.2, pred, coda::HPDinterval(mcmc(post.pred), prob=0.95), se)
  colnames(pred) <- c("Binomial.1.2", "fit", "lwrCI", "uprCI", "se")
  pred <- as_tibble(pred)

  colnames(post.pred) <- mam$Binomial.1.2

  return(list(post.pred, pred))
}

comb <- function(p, q) {
  list(rbind(p[[1]], q[[1]]),
       rbind(p[[2]], q[[2]]),
       rbind(p[[3]], q[[3]]))
}

cl <- makeCluster(cluster.size)
registerDoSNOW(cl)
timestamp()
tic()
pb <- txtProgressBar(max = n.trees, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
imputed <- foreach(i = 1:n.trees,
                   .packages = c('MCMCglmm', 'tibble'),
                   .inorder = FALSE,
                   .options.snow = opts,
                   .combine = comb,
                   .multicombine = FALSE) %dopar% mcmc.regression(i)
toc()
stopCluster(cl)
gc()

write_csv(as_data_frame(imputed[[1]]), "builds/test_imputed.densities_all.samples.csv")
write_csv(as_data_frame(imputed[[2]]), "builds/test_densities_fit.solution.csv")
write_csv(as_data_frame(imputed[[3]]), "builds/test_densities_post.pred.csv")
write_csv(as_data_frame(imputed[[3]]) %>% sample_n(10000), "builds/test_densities_post.pred.10k.sample.csv")

b <- as_data_frame(imputed[[1]])
b <- b %>% filter(Binomial.1.2 == "Abditomys_latidens")
a <- as_data_frame(imputed[[3]]) %>% sample_n(5000)
c <- as_data_frame(imputed[[3]]) %>% sample_n(5000)
d <- as_data_frame(imputed[[3]]) %>% sample_n(10000)
e <- as_data_frame(imputed[[3]]) %>% sample_n(10000)
ggplot(a[,1], aes(x = Abditomys_latidens)) +
  geom_density(aes(col = "post.pred.5k.1")) +
  geom_density(data = c[,1], aes(col = "post.pred.5k.2")) +
  geom_density(data = d[,1], aes(col = "post.pred.10k.1")) +
  geom_density(data = e[,1], aes(col = "post.pred.10k.2")) +
  geom_density(data = data_frame(Abditomys_latidens = rnorm(5000, mean(as.numeric(b$fit)), mean(as.numeric(b$se)))),
               aes(col = "sample.pred"))
