library(tidyverse)
library(gridExtra)
library(coda)
library(MCMCglmm)
library(ggpmisc)

i = 1
chain.1 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain1.rds"))
chain.2 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain2.rds"))
chain.3 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain3.rds"))
# 
# chain.1 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain1.alt.rds"))
# chain.2 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain2.alt.rds"))
# chain.3 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain3.alt.rds"))

### Checking 3 chains for tree 1
sol <- bind_rows(as.data.frame(chain.1$Sol[, 1:2]), 
                 as.data.frame(chain.2$Sol[, 1:2]), 
                 as.data.frame(chain.3$Sol[, 1:2]))
sol["chain"] <- gl(3, 333)
sol["sample"] <- rep(1:333, 3)
sol <- gather(sol, key = "variable", value = "value", -chain, -sample)

left <- ggplot(sol, aes(x = sample, y = value, col = chain)) +
  geom_line() + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, lty = "dotted", col = "black") +
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + 
  ylab("") + 
  stat_poly_eq(aes(label = ..adj.rr.label..), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3, vstep = 0, hstep = 0.20)
right <- ggplot(sol, aes(x = value, col = chain)) +
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + 
  labs(x = "", y = "")
p.main <- grid.arrange(left, right, nrow = 1)
ggsave("output/appendix1_fig5.png", p.main, width = 25.6, height = 14.4, units = "cm")

VCV <- bind_rows(as.data.frame(chain.1$VCV), 
                 as.data.frame(chain.2$VCV), 
                 as.data.frame(chain.3$VCV))
VCV["chain"] <- gl(3, 333)
VCV["sample"] <- rep(1:333, 3)
VCV <- gather(VCV, key = "variable", value = "value", -chain, -sample)

left <- ggplot(VCV, aes(x = sample, y = value, col = chain)) +
  geom_line() + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, lty = "dotted", col = "black") +
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + 
  ylab("") + 
  stat_poly_eq(aes(label = ..adj.rr.label..), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3, vstep = 0, hstep = 0.20)
right <- ggplot(VCV, aes(x = value, col = chain)) +
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + 
  labs(x = "", y = "")
p.random <- grid.arrange(left, right, nrow = 1)
ggsave("output/appendix1_fig6.png", p.random, width = 25.6, height = 14.4, units = "cm")

# Checking convergence for our fixed factors
gelman.diag(mcmc.list(chain.1$Sol[, 1:2], chain.2$Sol[, 1:2], chain.3$Sol[, 1:2]), autoburnin = FALSE)

# Checking convergence for our random terms
gelman.diag(mcmc.list(chain.1$VCV, chain.2$VCV, chain.3$VCV), autoburnin = FALSE)

# Gelman plot:
gelman.plot(mcmc.list(chain.1$VCV, chain.2$VCV, chain.3$VCV), autoburnin = FALSE)
gelman.plot(mcmc.list(chain.1$Sol[, 1:2], chain.2$Sol[, 1:2], chain.3$Sol[, 1:2]), autoburnin = FALSE)


### Checking effective sample size
chain.1.2.3.Sol <- list(chain.1$Sol[, 1:2], chain.2$Sol[, 1:2], chain.3$Sol[, 1:2])
chain.1.2.3.VCV <- list(chain.1$VCV, chain.2$VCV, chain.3$VCV)
effectiveSize(chain.1.2.3.Sol)/3
# G/R structure
effectiveSize(chain.1.2.3.VCV)/3


# acf plot for the fixed estimates
acf(chain.1$Sol[, 1], lag.max = 20)
acf(chain.1$Sol[, 2], lag.max = 20)
# acf plot for the first random term in our model (the phyl term)
acf(chain.1$VCV[, 1], lag.max = 20)
acf(chain.1$VCV[, 2], lag.max = 20)

lambda <- chain.1$VCV[,'Binomial.1.2'] / (chain.1$VCV[,'Binomial.1.2'] + chain.1$VCV[,'units'])
mean(lambda)
median(lambda)

posterior.mode(lambda)
HPDinterval(lambda)

# unused:
# summary(chain.1)
# 
# HPDinterval(chain.1$Sol[,1:2])
# test0 <- predict.MCMCglmm(chain.1, df, interval = "confidence")[1:n.mam, ]
# test1 <- predict.MCMCglmm(chain.1, df, interval = "confidence", marginal = NULL)[1:n.mam, ]
# 
# se.0 <- (test1[,3] - test1[,2])/2/qnorm(0.975)
# 
# post.pred <- MCMC.predict(chain.1, df)
# post.pred.t <- t(post.pred)
# post.pred.t <- post.pred.t[1:n.mam, ]
# se.1 <- apply(post.pred.t, 1, sd)
# hist(se.1)
# 
# se.test <- tibble(se = c(se.0, se.1), group = gl(2, length(se.0), labels = c("se.0", "se.1")))
# ggplot(se.test, aes(se, col = group)) +
#   geom_density() +
#   geom_rug()
# 
# a <- as_tibble(t(post.pred.t))
# colnames(a) <- mam$Binomial.1.2
# a <- a %>% gather(key = "Species", value = "log10Density")
# 
# ggplot(a %>% filter(Species %in% sample(a$Species, 100)), 
#        aes(x = log10Density, col = Species)) +
#   geom_density() +
#   scale_color_discrete(guide=FALSE)
# 
# coda::HPDinterval(post.pred.t[1,])
# manual code: https://rdrr.io/cran/coda/src/R/HPDinterval.R
# 
# test1[1,]
# a <- ecdf(post.pred.t[1,])
# quantile(a, c(.5, .025, .975))
# 
# prob = .95
# obj <- as.matrix(post.pred.t[1,])
# vals <- apply(obj, 2, sort)
# nsamp <- nrow(vals)
# gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
# init <- 1:(nsamp - gap)
# inds <-  which.min(vals[init + gap, ] - vals[init, ])
# ans <- cbind(vals[inds],
#              vals[inds + gap])
# dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
# attr(ans, "Probability") <- gap/nsamp
# ans
# 
# plot(pantheria$log10density ~ pantheria$log10BM, col = "red")
# points(test0[,1] ~ mam$log10BM, col = "blue")
# points(test1[,1] ~ mam$log10BM, col = "purple")
