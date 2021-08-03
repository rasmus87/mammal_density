# GLMM model diagnostics across all data
# 03/08-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(gridExtra)
library(coda) # For HPDinterval()
library(MCMCglmm)
library(ggpmisc) # For stat_poly_eq()



imputed <- read_csv("builds/densities_post.pred.csv")

imputed <- tibble(log10density.mean = rowMeans(imputed),
                  tree = gl(1000, 3), 
                  chain = as_factor(rep(1:3, 1000)))

left <- ggplot(imputed, aes(x = as.numeric(tree), y = log10density.mean, col = chain)) +
  geom_line() + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, lty = "dotted", col = "black") +
  theme_bw() +
  theme(legend.position="none") + 
  ylab("") + 
  stat_poly_eq(aes(label = ..adj.rr.label..), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3, vstep = 0, hstep = 0.20)
right <- ggplot(imputed, aes(x = log10density.mean, col = chain)) +
  geom_density() +
  geom_rug() +
  theme_bw() +
  theme(legend.position="none") + 
  labs(x = "", y = "")
p.main <- grid.arrange(left, right, nrow = 1)
ggsave("output/appendix1_figXXX.png", p.main, width = 25.6, height = 14.4, units = "cm")

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
