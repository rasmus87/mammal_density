library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)

# Run after data is loaded with "0.2 load.data.R"

# For 333 samples:
# 36 hours runtime for 1000 trees on 20 cores
# ~ 80 GB RAM required
# ~ 80 GB disk space

## Set options:
# Set parralell cluster size
# cluster.size <- 1
cluster.size <- 6
# cluster.size <- 20
# How many trees do you want to run this for? 2-1000?
# n.trees <- 1
# n.trees <- 6
# n.trees <- 6*10
n.trees <- 1000
# Number of mcmc samples per (1000 trees)
# Run 333 for good chains for testing convergence
# Run 3 samples for actual data is enough
mcmc.samples <- 3

prior <- list(G = list(G1 = list(V = 1, nu = 0.02)), 
              R = list(V = 1, nu = 0.02))
thin <- 75
burnin <- thin * 10
nitt <- mcmc.samples * thin + burnin
i = 1
mcmc.regression <- function(i) {
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
  if(i == 1 & mcmc.samples == 333) {
    saveRDS(chain.1, paste0("builds/mcmcglmms/tree", i, ".chain1.alt.rds"), compress = FALSE)
    saveRDS(chain.2, paste0("builds/mcmcglmms/tree", i, ".chain2.alt.rds"), compress = FALSE)
    saveRDS(chain.3, paste0("builds/mcmcglmms/tree", i, ".chain3.alt.rds"), compress = FALSE)
  }
  
  gc()
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

write_csv(as_tibble(imputed[[1]]), paste0("builds/", mcmc.samples ,"_densities_fit.solution.alt.csv"))
write_csv(as_tibble(imputed[[2]]), paste0("builds/", mcmc.samples ,"_densities_post.pred.alt.csv"))
