# Run GLMM model
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)

# Normal PanTHERIA dataset
# Turn pantheria into data.frame for MCMCglmm
density.dataset <- read.csv("builds/imputation_dataset_PanTHERIA.csv")

# Load clean dataset of all species
mam <- read.csv("builds/imputation_dataset_mam.csv")
# Number of species
n.mam <- nrow(mam)

# Cross validation: 5-fold
# Set seed and make a random ordered data vector
set.seed(42)
n <- nrow(density.dataset)
density.dataset <- density.dataset[sample(n), ]

# Create 10 equally size folds
folds <- cut(1:n, breaks = 5, labels = FALSE)

# Load forest
forest <- read_rds("builds/forest.rds")

# For 333 samples:
# 36 hours runtime for 1000 trees on 20 cores
# ~ 80 GB RAM required
# ~ 80 GB disk space

## Set options:
# Set parralell cluster size
# cluster.size <- 2
cluster.size <- 30
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
burnin <- 5000 * 2


# Run actual chains -------------------------------------------------------

# Set samples and iterations
# Run 1 sample per chain per tree for 1000 trees
mcmc.samples <- 1
nitt <- burnin + (mcmc.samples - 1) * thin + 1

mcmc.regression <- function(i, train.data, test.data) {
  tree <- forest[[i]]
  inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
  chain.1 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = train.data, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.2 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = train.data, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.3 <- MCMCglmm(log10density ~ log10BM, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = train.data, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  
  pred1 <- MCMC.predict(chain.1, test.data, train.data)
  pred2 <- MCMC.predict(chain.2, test.data, train.data)
  pred3 <- MCMC.predict(chain.3, test.data, train.data)
  
  post.pred <- rbind(pred1, pred2, pred3)
  
  return(post.pred)
}

MCMC.predict <- function(object, test.data, train.data) {
  newdata <- rbind(test.data, train.data)
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
  
  post.pred <- t(apply(object$Sol, 1, function(x){(W %*% x)@x}))[, 1:nrow(test.data)]
  
  names(post.pred) <- test.data$Binomial.1.2
  
  return(post.pred)
}

cl <- makeCluster(cluster.size)
registerDoSNOW(cl)
set.seed(42)
timestamp()
tic()

# Perform 5 fold cross validation
res <- tibble()
for(fold in 1:5) {
  print(paste0("Cross-validation", fold, "/5"))
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds == fold)

  train.data <- density.dataset[-testIndexes, ]
  
  # Prediction dataset
  test.data <- mam[which(mam$Binomial.1.2 %in% density.dataset$Binomial.1.2[testIndexes]), ]

  # Impute
  pb <- txtProgressBar(max = n.trees, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  imputed <- foreach(i = 1:n.trees,
                     .packages = c('MCMCglmm', 'tibble'),
                     .inorder = FALSE,
                     .options.snow = opts,
                     .combine = rbind,
                     .multicombine = TRUE) %dopar% mcmc.regression(i, train.data, test.data)
  imputed.mean <- colMeans(imputed)
  imputed.long <- tibble(Binomial.1.2 = names(imputed.mean),
                         log10.density.mean = imputed.mean,
                         fold = fold)
  res <- rbind(res, imputed.long)
}
toc()
stopCluster(cl)
gc()

total.res <- res %>% 
  left_join(density.dataset %>% transmute(Binomial.1.2, log10.density.pantheria = log10density))


write_csv(total.res, "builds/densities_5xcross_val.csv")
total.res <- read_csv("builds/densities_5xcross_val.csv")

# Persons R-squared and RMSE
total.res %>% 
  group_by(fold) %>% 
  summarise(persons.r2 = cor(log10.density.mean, log10.density.pantheria)^2,
            rmse = sqrt(mean((log10.density.pantheria - log10.density.mean)^2))) %>% 
  summarise_at(c("persons.r2", "rmse"), mean)

total.res %>% 
  group_by(fold) %>% 
  summarise(persons.r2 = cor(10^log10.density.mean, 10^log10.density.pantheria)^2,
            rmse = sqrt(mean((10^log10.density.pantheria - 10^log10.density.mean)^2))) %>% 
  summarise_at(c("persons.r2", "rmse"), mean)

# total.res <- total.res %>% 
#   mutate(error = abs(log10.density.pantheria - log10.density.mean)^2,
#          error10 = abs(10^log10.density.pantheria - 10^log10.density.mean)^2)
# ggplot(total.res, aes(log10.density.mean, log10.density.pantheria, col = error10)) +
#   geom_point() +
#   geom_smooth() +
#   geom_abline(slope = 1) +
#   scale_color_viridis_c()
