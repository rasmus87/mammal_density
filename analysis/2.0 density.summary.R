library(tidyverse)

dataset <- read_csv("builds/imputation_dataset.csv")

# imputed <- read_csv("builds/333_densities_post.pred.10k.sample.csv")
# imputed <- read_csv("builds/33_densities_post.pred.10k.sample.csv")
imputed <- read_csv("builds/3_densities_post.pred.csv")

dens <- imputed %>% gather("Binomial.1.2", "log10density")

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")

mam <- mam %>% filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
                      !Family.1.2 %in% c(whale.families, seal.families),
                      !Binomial.1.2 %in% marine.carnivores)

imputed <- as.matrix(imputed)
attr(imputed, "class") <- "mcmc"
imputed.ci <- coda::HPDinterval(imputed)
imputed.ci <- imputed.ci %>% as_tibble(rownames = "Binomial.1.2")
imputed.ci <- imputed.ci %>% transmute(Binomial.1.2,
                                       log10.lower.95hpd = lower,
                                       log10.upper.95hpd = upper)
imputed <- as.matrix(imputed)
imputed.q95 <- apply(imputed, 2, quantile, probs = c(0.025, 0.975))
imputed.q95 <- t(imputed.q95)
imputed.q95 <- as_tibble(imputed.q95, rownames = "Binomial.1.2")
colnames(imputed.q95)[2:3] <- c("q.025", "q.975")

test <- left_join(imputed.ci, imputed.q95)
test <- test %>% 
  mutate(low.diff = log10.lower.95hpd - q.025, 
         high.diff = log10.upper.95hpd - q.975)
ggplot(test, aes()) +
  geom_density(aes(log10.lower.95hpd, col = "95 % HPD")) +
  geom_density(aes(log10.upper.95hpd, col = "95 % HPD")) +
  geom_density(aes(q.025, col = "95 % Quantile")) +
  geom_density(aes(q.975, col = "95 % Quantile"))
  
  

dens.summary <- dens %>%
  group_by(Binomial.1.2) %>% 
  summarise(log10.density.median = median(log10density),
            log10.density.mean = mean(log10density),
            sd = sd(log10density)) %>% 
  left_join(imputed.ci, by = "Binomial.1.2")

mam.dens <- mam %>% 
  transmute(Binomial.1.2, Order.1.2, Family.1.2, log10BM = log10(Mass.g)) %>% 
  left_join(dens.summary, by = "Binomial.1.2")
mam.dens <- mam.dens %>% mutate(density.median= 10^log10.density.median,
                                density.mean = 10^log10.density.mean,
                                lower.95hpd = 10^log10.lower.95hpd,
                                upper.95hpd = 10^log10.upper.95hpd)

# write_csv(mam.dens, "builds/imputed.density_333.csv")
# write_csv(mam.dens, "builds/imputed.density_33.csv")
write_csv(mam.dens, "builds/imputed.density.csv")
mam.dens <- read_csv("builds/imputed.density.csv")








ggplot(mam.dens, aes(x = log10BM, col = Order.1.2)) +
  geom_linerange(aes(ymin = log10.lower.95hpd, ymax = log10.upper.95hpd)) +
  geom_point(aes(y = log10.density.mean)) +
  theme_bw() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~(km^2)))

col9 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
          "#ff7f00", "#cab2d6")
col27 <- rep(col9, each = 3)
pch3 <- c(1, 2, 3)
pch27 <- rep(pch3, times = 9)

ggplot(mam.dens, aes(x = log10BM, col = Order.1.2, shape = Order.1.2)) +
  geom_point(aes(y = log10.density.mean)) +
  theme_bw() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~(km^2))) +
  scale_color_manual(values = col27) +
  scale_shape_manual(values = pch27) +
  theme(legend.position = c(1, 1), legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(1,1)) +
  guides(col = guide_legend(nrow = 9)) +
  ylim(NA, 5.25) +
  xlim(NA, 7.5)
ggsave("figures/dens_plot.png", width = 25.6, height = 14.4, units = "cm")


ggplot(mam.dens, aes(x = log10BM, col = Binomial.1.2 %in% dataset$Binomial.1.2)) +
  geom_point(aes(y = log10.density.mean), pch = 19) +
  theme_bw() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~(km^-2))) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "none")
ggsave("output/appendix1_fig1.png")

full.data <- dataset %>%
  transmute(Binomial.1.2, log10.density.pantheria = log10density) %>% 
  right_join(mam.dens)

ggplot(full.data, aes(x = log10.density.pantheria, y = log10.density.median, col = log10BM)) +
  geom_point(pch = 19) +
  geom_abline(slope = 1, lty = 2, lwd = .7) +
  geom_smooth(method = lm, se = F, lty = 1, col = "black", lwd = .9) +
  theme_bw() +
  xlab(expression(log[10]~PanTHERIA~empirical~density~(km^-2))) +
  ylab(expression(log[10]~Imputed~median~density~(km^-2))) +
  scale_color_continuous(expression(log[10]~Body~mass~(g))) +
  coord_equal(xlim = c(-2, 5), ylim = c(-2, 5))
ggsave("output/appendix1_fig2.png")

ggplot(full.data, aes(x = log10BM, y = log10.density.median-log10.density.pantheria, 
                      col = Order.1.2, shape = Order.1.2)) +
  geom_point() +
  geom_abline(slope = 0, lty = 2, lwd = .7) +
  geom_smooth(method = lm, se = F, lty = 1, col = "black", lwd = .9, aes(shape = NULL)) +
  theme_bw() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~difference~(km^-2))) +
  scale_color_manual(values = col27) +
  scale_shape_manual(values = pch27)
ggsave("output/appendix1_fig3.png")

#### TESTS ####
test <- as.matrix(imputed)[,1]
exp((log(10) * sd(test))^2/2) * 10^mean(test)
mean(10^test)

test.mcmc <- test
attr(test.mcmc, "class") <- "mcmc"
10^coda::HPDinterval(test.mcmc)
coda::HPDinterval(10^test.mcmc)

10^quantile(test, c(.025, 0.975))
quantile(10^test, c(.025, 0.975))

ggplot(data_frame(), aes(x = 10^test)) +
  geom_density() +
  geom_vline(xintercept = 10^quantile(test, c(.025, 0.975)), lty = 2, col = "blue") +
  geom_vline(xintercept = coda::HPDinterval(10^test.mcmc)[1:2], lty = 2, col = "red") +
  geom_vline(xintercept = mean(10^test), lty = 1, col = "red", lwd = 1) +
  geom_vline(xintercept = median(10^test), lty = 1, col = "blue", lwd = 1)

  
ggplot(data_frame(), aes(x = test)) +
  geom_density() +
  geom_vline(xintercept = quantile(test, c(.025, 0.975)), lty = 2, col = "blue") +
  geom_vline(xintercept = coda::HPDinterval(test.mcmc)[1:2], lty = 2, col = "red") +
  geom_vline(xintercept = mean(test), lty = 1, col = "red", lwd = 1) +
  geom_vline(xintercept = median(test), lty = 1, col = "blue", lwd = 1)
