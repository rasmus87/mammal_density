# GLMM model diagnostics
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(gridExtra)
library(coda) # For HPDinterval()
library(MCMCglmm)
library(ggpmisc) # For stat_poly_eq()


chain.1 <- read_rds(paste0("builds/mcmcglmms/chain1.rds"))
chain.2 <- read_rds(paste0("builds/mcmcglmms/chain2.rds"))
chain.3 <- read_rds(paste0("builds/mcmcglmms/chain3.rds"))

# chain.1 <- read_rds(paste0("builds/mcmcglmms/chain1.alt.rds"))
# chain.2 <- read_rds(paste0("builds/mcmcglmms/chain2.alt.rds"))
# chain.3 <- read_rds(paste0("builds/mcmcglmms/chain3.alt.rds"))

### Checking 3 chains
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



# Across all data ---------------------------------------------------------


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
ggsave("output/appendix1_figAllTrees1Sample.png", p.main, width = 25.6, height = 14.4, units = "cm")

# Checking convergence for our fixed factors
x1 <- as.matrix(imputed %>% filter(chain == 1) %>% pull(log10density.mean))
attr(x1, "class") <- "mcmc"
x2 <- as.matrix(imputed %>% filter(chain == 2) %>% pull(log10density.mean))
attr(x2, "class") <- "mcmc"
x3 <- as.matrix(imputed %>% filter(chain == 3) %>% pull(log10density.mean))
attr(x3, "class") <- "mcmc"

gelman.diag(mcmc.list(x1, x2, x3), autoburnin = FALSE)


### Checking effective sample size
chain.1.2.3 <- list(x1, x2, x3)
effectiveSize(chain.1.2.3.Sol)/3

