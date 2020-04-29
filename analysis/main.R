# Load libraries

library(tidyverse)

family.groups <- read_csv("data/ancestor.tree.DNA.WR05.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())
mam <- mam %>%
  # Order Cetacea (Whales s.l.)
  filter(!Family.1.2 %in% c("Balaenidae", "Balaenopteridae", "Delphinidae",
                            "Eschrichtiidae", "Iniidae", "Monodontidae",
                            "Neobalaenidae", "Phocoenidae", "Physeteridae",
                            "Ziphiidae")) %>%
  # Families in the clade Pinnipedia (Seal s.l.):
  # Odobenidae (walruses)
  # Otariidae (fur seals and sea lions)
  # Phocidae (true seals)
  filter(!Family.1.2 %in% c("Odobenidae", "Otariidae", "Phocidae")) %>%
  # Order Sirenia (Sea cows s.l.):
  filter(Order.1.2 != "Sirenia") %>%
  # Order Chiroptera (Bats)
  filter(Order.1.2 != "Chiroptera")

median.family.distance <- readRDS("builds/median.family.distance.rds")

does.not.exist <- which(!row.names(median.family.distance) %in% family.groups$families)

median.family.distance <- median.family.distance[-does.not.exist, ]
closest <- apply(median.family.distance, 2, which.min)
dist <- apply(median.family.distance, 2, min)
fam <- tibble(family = colnames(median.family.distance),
              closest.family = rownames(median.family.distance)[closest],
              dist = dist)

mam <- left_join(mam, fam, by = c("Family.1.2" = "family"))

taxonomy <- read_csv("data/ancestor.tree.DNA.WR05.csv", col_types = cols())

mam <- left_join(mam, taxonomy, by = c("closest.family" = "families"))

load("data/dens.DNA.tax.model.RData")

mam$log10_Mass <- log10(mam$Mass.g)

se.fit <- predict.lm(dens.DNA.tax.model, mam, se.fit = TRUE)$se.fit
res.ci <- predict.lm(dens.DNA.tax.model, mam, interval = "confidence")
res.pi <- predict.lm(dens.DNA.tax.model, mam, interval = "prediction")
res <- cbind(res.ci, res.pi[, 2:3], 10^res.ci, se.fit)
res <- as.tibble(res)
names(res) <- c("log10dens.est", 
                "log10dens.est_CIlwr", "log10dens.est_CIupr", 
                "log10dens.est_PIlwr", "log10dens.est_PIupr",
                "dens.est", 
                "dens.est_CIlwr", "dens.est_CIupr",
                "se.fit")

mam <- bind_cols(mam, res)

out <- mam %>% select(Binomial.1.2, 
                      log10dens.est,
                      se.fit,
                      log10dens.est_CIlwr, log10dens.est_CIupr,
                      log10dens.est_PIlwr, log10dens.est_PIupr,
                      dens.est, 
                      dens.est_CIlwr, dens.est_CIupr)

write_csv(out, "output/animal.density.km2.csv")




pantheria <- read_tsv("data/PanTHERIA_1-0_WR05_Aug2008.txt", col_types = cols())
names(pantheria) <- make.names(names(pantheria))
pantheria <- pantheria %>%
  filter(!is.na(X21.1_PopulationDensity_n.km2) & !is.na(X5.1_AdultBodyMass_g))

# Exclude non terrestrial mammals:
pantheria <- pantheria %>%
  # Order Cetacea (Whales s.l.)
  filter(MSW05_Order != "Cetacea") %>%
  # Families in the clade Pinnipedia (Seal s.l.):
  # Odobenidae (walruses)
  # Otariidae (fur seals and sea lions)
  # Phocidae (true seals)
  filter(!MSW05_Family %in% c("Odobenidae", "Otariidae", "Phocidae")) %>%
  # Order Sirenia (Sea cows s.l.):
  filter(MSW05_Order != "Sirenia") %>%
  # Other marine mammals:
  # Marine otter (Lontra felina)
  filter(MSW05_Binomial != "Lontra felina") %>%
  # Sea otter (Enhydra lutris)
  filter(MSW05_Binomial != "Enhydra lutris") %>%
  # Polar bear (Ursus maritimus)
  filter(MSW05_Binomial != "Ursus maritimus") %>%
  # Order Chiroptera (Bats)
  filter(MSW05_Order != "Chiroptera")

taxonomy <- read_csv("data/ancestor.tree.DNA.WR05.csv", col_types = cols())
# Find missing families
unique(pantheria$MSW05_Family[which(!pantheria$MSW05_Family %in% taxonomy$families)])
# [1] "Indriidae" "Aotidae"
# Add closely related families for replacements:
pantheria$MSW05_Family[pantheria$MSW05_Family == "Indriidae"] <- "Indridae"
pantheria$MSW05_Family[pantheria$MSW05_Family == "Aotidae"] <- "Cebidae"

# Add ancestral divergence points to all species
pantheria <- left_join(pantheria, taxonomy, by = c("MSW05_Family" = "families"))
pantheria["log10_Mass"] <- log10(pantheria$X5.1_AdultBodyMass_g)
pantheria["log10_Density"] <- log10(pantheria$X21.1_PopulationDensity_n.km2)


ggplot(mam, aes(log10_Mass, log10(fit))) +
  geom_point() +
  geom_point(data = pantheria, aes(log10_Mass, log10_Density), col = "red") 

ggplot(mam %>% filter(Genus.1.2 %in% c("Dendrolagus", "Macropus")), aes(log10_Mass, log10(fit), col = Genus.1.2)) +
  geom_line() +
  geom_point(data = pantheria %>% filter(MSW05_Genus %in% c("Dendrolagus", "Macropus")), aes(log10_Mass, log10_Density, col = MSW05_Genus)) 


ggplot(mam, aes(log10_Mass, (fit))) +
  geom_point() +
  geom_point(data = pantheria, aes(log10_Mass, X21.1_PopulationDensity_n.km2), col = "red") +
  scale_y_continuous(limits = c(0, 5000))

pantheria %>% filter(X21.1_PopulationDensity_n.km2 > 5000) %>% select(MSW05_Order, MSW05_Family, MSW05_Binomial, X5.1_AdultBodyMass_g, X21.1_PopulationDensity_n.km2) %>% 
  print(n = 25)