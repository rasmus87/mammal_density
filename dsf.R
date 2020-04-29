library(tidyverse)
imp.dens <- read_csv("builds/imputed.density.csv")
imp.dens <- imp.dens %>% mutate(imp.dens = 10^log10density.median)
mod.dens <- read_csv("output/animal.density.km2.csv")
mod.dens <- mod.dens %>% mutate(mod.dens = dens.est)

df <- full_join(imp.dens, mod.dens)

ggplot(df, aes(mod.dens, imp.dens, col = Order.1.2)) +
  geom_point() + 
  geom_smooth(method = "lm", aes(col = NULL)) +
  geom_abline(slope = 1) +
  coord_equal()

ggplot(df %>% filter(Order.1.2 == "Proboscidea"), aes(mod.dens, imp.dens, col = Binomial.1.2)) +
  geom_point() + 
  geom_smooth(method = "lm", aes(col = NULL)) +
  geom_abline(slope = 1) +
  coord_equal() +
  xlim(0, NA) +
  ylim(0, NA)

df %>% filter(Binomial.1.2 %in% c("Bos_primigenius", "Equus_ferus", "Cervus_elaphus")) %>% select(Binomial.1.2, mod.dens, imp.dens)

df <- df %>% mutate(diff = mod.dens-imp.dens)

ggplot(df, aes(mod.dens, diff, col = Order.1.2)) +
  geom_point()


ggplot(df, aes(log10BM, diff, col = Order.1.2)) +
  geom_point()

ggplot(df, aes(log10BM, log10density.median, col = Order.1.2 == "Carnivora")) +
  geom_point() + 
  geom_smooth(method = "lm", aes(col = NULL))

df %>% filter(Family.1.2 %in% c("Ursidae")) %>% select(Binomial.1.2, mod.dens, imp.dens)
