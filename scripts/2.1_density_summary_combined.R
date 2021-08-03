# Create density summary
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(ggpmisc)
library(ggtext)

# Load dataset only based on PanTHERIA
dataset <- read_csv("builds/imputation_dataset_PanTHERIA.csv")
imputed <- read_csv("builds/densities_post.pred.csv")

# Load dataset based on PanTHERIA mixed with TetraDENSITY
# dataset <- read_csv("builds/imputation_dataset_PanTHERIA_TetraDENSITY.csv")
# imputed <- read_csv("builds/densities_post.pred.alt.csv")

# Make dataset long
dens <- imputed %>% 
  pivot_longer(cols = everything(), 
               names_to = "Binomial.1.2",
               values_to = "log10density")

# Load PHYLACINE 1.2.1
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Filter to our imputed species list
mam <- mam %>% filter(Binomial.1.2 %in% dens$Binomial.1.2)

# Turn imputed data.frame into a mcmc object for calculating confidence intervals
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
# Almost the same thing! But we keep the HPD interval

# Estimate
dens.summary <- dens %>%
  group_by(Binomial.1.2) %>% 
  summarise(log10.density.median = median(log10density),
            log10.density.mean = mean(log10density),
            sd = sd(log10density)) %>% 
  left_join(imputed.ci, by = "Binomial.1.2")

mam.dens <- mam %>% 
  transmute(Binomial.1.2, Order.1.2, Family.1.2, log10BM = log10(Mass.g)) %>% 
  left_join(dens.summary, by = "Binomial.1.2")
mam.dens <- mam.dens %>% mutate(density.median = 10^log10.density.median,
                                density.mean = 10^log10.density.mean,
                                lower.95hpd = 10^log10.lower.95hpd,
                                upper.95hpd = 10^log10.upper.95hpd)

write_csv(mam.dens, "output/Table S4 Imputed density.csv")


# Build supplementary figures demonstrating test of the imputed result
mam.dens <- read_csv("output/Table S4 Imputed density.csv")
mam.dens0 <- read_csv("output/Table S4 Imputed density - Copy.csv")

# Test difference to old estimates
mam.dens.test <- left_join(mam.dens, mam.dens0, by = c("Binomial.1.2", "Order.1.2", "Family.1.2"))# %>%
  #filter(Order.1.2 == "Proboscidea")
test <- mam.dens.test[4:9] - mam.dens.test[14:19]
test <- test %>% pivot_longer(everything(), names_to = "parameter")
ggplot(test, aes(value, col = parameter)) +
  geom_density()
test %>% group_by(parameter) %>% summarise(sum(value), mean(value), median(value), sd(value))

theme_R <- function() {
  theme_bw() %+replace% 
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
}

# ggplot(mam.dens, aes(x = log10BM, col = Order.1.2)) +
#   geom_linerange(aes(ymin = log10.lower.95hpd, ymax = log10.upper.95hpd)) +
#   geom_point(aes(y = log10.density.mean)) +
#   theme_R() +
#   xlab(expression(log[10]~Body~mass~(g))) +
#   ylab(expression(log[10]~Density~(km^2)))

col9 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
          "#ff7f00", "#cab2d6")
col27 <- rep(col9, each = 3)
pch3 <- c(1, 2, 3)
pch27 <- rep(pch3, times = 9)
orders <- sort(unique(mam.dens$Order.1.2))
labels <- paste0("<span style='color:", col27, "'>", orders, "</span>")

ggplot(mam.dens, aes(x = log10BM, col = Order.1.2, shape = Order.1.2)) +
  geom_point(aes(y = log10.density.mean)) +
  theme_R() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~(km^-2))) +
  scale_color_manual(values = col27, breaks = orders, labels = labels) +
    scale_shape_manual(values = pch27, breaks = orders, labels = labels) +
  theme(legend.text = element_markdown()) +
  theme(legend.position = c(1, 1), 
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(1,1)) +
  guides(col = guide_legend(nrow = 9)) +
  ylim(NA, 4.30)
ggsave("output/appendix1_fig1.png", width = 25.6, height = 16.5, units = "cm")

ggplot(mam.dens, aes(x = log10BM, col = Binomial.1.2 %in% dataset$Binomial.1.2)) +
  geom_point(aes(y = log10.density.mean), pch = 19, cex = 1) +
  theme_R() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~(km^-2))) +
  scale_color_brewer(palette = "Paired", name = "Imputed density for species",
                     labels = c("without known density",
                                "with empirical known density")) +
  theme(legend.position = c(1, 1), 
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(1.5, 1.5)) +
  ylim(NA, 4.30)
ggsave("output/appendix1_fig2.png", width = 25.6, height = 16.5, units = "cm")

full.data <- dataset %>%
  transmute(Binomial.1.2, log10.density.pantheria = log10density) %>% 
  right_join(mam.dens)

ggplot(full.data %>% filter(!is.na(log10.density.pantheria)), 
       aes(x = log10.density.pantheria, y = log10.density.median, col = log10BM)) +
  geom_point(pch = 19) +
  geom_abline(slope = 1, lty = 2, lwd = .7) +
  geom_smooth(formula = y ~ x, method = lm, se = F, lty = 1, col = "black", lwd = .9) +
  theme_R() +
  xlab(expression(log[10]~PanTHERIA~empirical~density~(km^-2))) +
  ylab(expression(log[10]~Imputed~median~density~(km^-2))) +
  scale_color_continuous(expression(log[10]~Body~mass~(g))) +
  coord_equal(xlim = c(-2, 5), ylim = c(-2, 5)) + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*\' ,  \'*")), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3)
ggsave("output/appendix1_fig3.png", width = 25.6, height = 14.4, units = "cm")

ggplot(full.data %>% filter(!is.na(log10.density.pantheria)), 
       aes(x = log10BM, y = log10.density.median - log10.density.pantheria, 
           col = Order.1.2, shape = Order.1.2)) +
  geom_point() +
  geom_abline(slope = 0, lty = 2, lwd = .7) +
  # geom_smooth(formula = y ~ x, method = lm, se = F, lty = 1, col = "black", lwd = .9, aes(shape = NULL)) +
  theme_R() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~difference~(km^-2))) +
  scale_color_manual(values = col27, breaks = orders) +
  scale_shape_manual(values = pch27, breaks = orders) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*\' ,  \'*"),
                   col = NULL, shape = NULL), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3) +
  guides(col = guide_legend(ncol = 1))
ggsave("output/appendix1_fig4.png", width = 25.6, height = 14.4, units = "cm")

