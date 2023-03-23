# Check central tendency measures of our results

# Load libraries
library(tidyverse)

# Load dataset and imputed results
dataset <- read_csv("builds/imputation_dataset_PanTHERIA.csv")
imputed <- read_csv("builds/densities_post.pred.csv")

# Make dataset long
dens <- imputed %>%
  pivot_longer(cols = everything(), 
               names_to = "Binomial.1.2",
               values_to = "log10density")

# Summarise central tendency across 3000 samples per species
dens.means <- dens %>% 
  group_by(Binomial.1.2) %>% 
  summarise(ari.mean = mean(10^log10density),
            geo.mean = 10^mean(log10density),
            median = 10^median(log10density))

# Calculate deviation fraction from empirical value
stats <- dataset %>% 
  mutate(empirical.dens = 10^log10density) %>% 
  left_join(dens.means) %>% 
  mutate(ari.frac = 1 - ari.mean/empirical.dens,
         geo.frac = 1 - geo.mean/empirical.dens,
         med.frac = 1 - median/empirical.dens)

# Calculate pseudo R-squared
TSS <- sum((stats$log10density - mean(stats$log10density))^2)
RSS <- sum((stats$log10density - log10(stats$geo.mean))^2)
1 - RSS/TSS

# Make format long for plotting and further calculations
stats.long <- stats %>% 
  pivot_longer(cols = ends_with("frac"), names_to = "estimate")

# Plot showing distribution of deviations from empircal measures
ggplot(stats.long, aes(x = value, col = estimate)) + 
  geom_density() +
  geom_rug()
# Medians and geometric mean shows the same thing
# All of the outliers are overestiamtes

# Show the worst cases
stats.long %>% 
  arrange(-abs(value)) %>% 
  filter(!duplicated(Binomial.1.2)) %>% 
  filter(abs(value) > 10) %>% 
  print(n=100)

# Show stats
stats.long %>% 
  group_by(estimate) %>% 
  summarise(mean = mean(value),
            median = median(value),
            sd = sd(value),
            fraction.over = mean(value < 0))


stats %>% 
  pivot_longer(cols = c(empirical.dens, ari.mean, geo.mean), names_to = "estimate") %>% 
  ggplot(aes(x = value, col = estimate)) + 
  geom_density() +
  geom_rug() +
  scale_x_log10()



stats.long %>% 
  filter(Family.1.2 == "Talpidae",
         estimate == "geo.frac") %>% 
  summarise(median(value))


stats.long %>% 
  filter(Family.1.2 == "Canidae",
         estimate == "geo.frac")
