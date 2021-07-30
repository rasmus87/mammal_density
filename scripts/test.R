mcmc.samples <- 333

prior <- list(G = list(G1 = list(V = 1, nu = 0.02)), 
              R = list(V = 1, nu = 0.02))
thin <- 75
burnin <- 1000
nitt <- mcmc.samples * thin + burnin
i = 1
tree <- forest[[i]]
inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
chain.1 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = pantheria, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE,
                    verbose = TRUE)

pred1 <- MCMC.predict(chain.1, mam)

