# Run GLMM model
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)

# Switch between PanTHERIA or alternative dataset -------------------------

dataset = "" # Normal dataset
# dataset = ".alt" # Alternative dataset
if(dataset == "") {
  # Normal PanTHERIA dataset
  # Turn pantheria into data.frame for MCMCglmm
  density.dataset <- read.csv("builds/imputation_dataset_PanTHERIA.csv")
  
} else {
  # Alternative mixed dataset
  # Turn pantheria into data.frame for MCMCglmm
  density.dataset <- read.csv("builds/imputation_dataset_PanTHERIA_TetraDENSITY.csv")
}

# Load clean dataset of all species
mam <- read.csv("builds/imputation_dataset_mam.csv")
# Number of species
n.mam <- nrow(mam)

# Combine datasets
df <- rbind(mam, pantheria)

# Load forest
forest <- read_rds("builds/forest.rds")

# For 333 samples:
# 36 hours runtime for 1000 trees on 20 cores
# ~ 80 GB RAM required
# ~ 80 GB disk space

## Set options:
# Set parralell cluster size
# cluster.size <- 2
cluster.size <- 6
# cluster.size <- 20
# How many trees do you want to run this for? 2-1000?
# n.trees <- 2
# n.trees <- 6
# n.trees <- 6*10
n.trees <- 1000

# Set priors
prior <- list(G = list(G1 = list(V = 1, nu = 0.02)), 
              R = list(V = 1, nu = 0.02))
thin <- 100
burnin <- 3000


# Chain test --------------------------------------------------------------

# Set samples and iterations
# Run 333 for good chains for testing convergence
mcmc.samples <- 333
nitt <- burnin + mcmc.samples * thin

# For being able to rerun on the same data and get the same result
set.seed(42)

tree <- forest[[1]]
inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
chain.1 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = density.dataset, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE,
                    verbose = FALSE)
chain.2 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = density.dataset, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE,
                    verbose = FALSE)
chain.3 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = density.dataset, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE,
                    verbose = FALSE)
write_rds(chain.1, paste0("builds/mcmcglmms/tree", i, "chain1", alt, ".rds"))
write_rds(chain.2, paste0("builds/mcmcglmms/tree", i, "chain2", alt, ".rds"))
write_rds(chain.3, paste0("builds/mcmcglmms/tree", i, "chain3", alt, ".rds"))









# Set samples and iterations
# Run 1 sample per tree for 1000 trees
mcmc.samples <- 1
nitt <- burnin + mcmc.samples * thin


i = 1
mcmc.regression <- function(i) {
  set.seed(42)
  tree <- forest[[i]]
  inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
  chain.1 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = density.dataset, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.2 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = density.dataset, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.3 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = density.dataset, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)

  pred1 <- MCMC.predict(chain.1, df)
  pred2 <- MCMC.predict(chain.2, df)
  pred3 <- MCMC.predict(chain.3, df)

  post.pred <- rbind(pred1, pred2, pred3)
  
  solution <- rbind(chain.1$Sol[, 1:2],
                    chain.2$Sol[, 1:2],
                    chain.3$Sol[, 1:2])
  solution <- as.data.frame(solution)
  random.effects <- 3:ncol(chain.1$Sol)
  solution["random.effect"] <- c(rowMeans(chain.1$Sol[, random.effects]),
                                 rowMeans(chain.2$Sol[, random.effects]),
                                 rowMeans(chain.3$Sol[, random.effects]))
  solution["tree"] <- i
  solution["chain"] <- as.numeric(gl(3, mcmc.samples))
  
  return(list(solution, post.pred))
}

MCMC.predict <- function(object, newdata) {
  pred.object <- MCMCglmm(fixed=object$Fixed$formula, 
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
  
  W <- cbind(pred.object$X, pred.object$Z)
  
  post.pred <- t(apply(object$Sol, 1, function(x){(W %*% x)@x}))[, 1:n.mam]

  colnames(post.pred) <- mam$Binomial.1.2

  return(post.pred)
}

comb <- function(...) {
  args <- list(...)
  lapply(seq_along(args[[1]]), function(i)
    do.call('rbind', lapply(args, function(a) a[[i]])))
}

cl <- makeCluster(cluster.size)
registerDoSNOW(cl)
set.seed(42)
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
                   .multicombine = TRUE) %dopar% mcmc.regression(i)
toc()
stopCluster(cl)
gc()

write_csv(as_tibble(imputed[[1]]), paste0("builds/", mcmc.samples ,"_densities_fit.solution.csv"))
write_csv(as_tibble(imputed[[2]]), paste0("builds/", mcmc.samples ,"_densities_post.pred.csv"))
# write_csv(as_tibble(imputed[[2]]) %>% sample_n(9000), paste0("builds/", mcmc.samples ,"_densities_post.pred.9k.sample.csv"))