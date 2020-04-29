pantheria <- read_tsv("../metabolic_rate/data/PanTHERIA_1-0_WR05_Aug2008.txt", col_types = cols())
names(pantheria) <- make.names(names(pantheria))
pantheria <- pantheria %>% 
  filter(!is.na(X21.1_PopulationDensity_n.km2)) %>% 
  transmute(Binomial.Pantheria = paste(MSW05_Genus, MSW05_Species, sep = "_"),
            density = X21.1_PopulationDensity_n.km2,
            source = "PanTHERIA.2008")
# Load taxonomy solver
syn <- read_csv("../PHYLACINE_1.1/Data/Taxonomy/Synonymy_table_with_unaccepted_species.csv", col_types = cols(), guess_max = 5000)
syn <- syn %>% 
  transmute(Binomial.1.2,
            Binomial.EltonTraits.1.0 = paste(EltonTraits.1.0.Genus, EltonTraits.1.0.Species, sep = "_"),
            Binomial.1.0 = paste(Genus.1.0, Species.1.0, sep = "_"),
            Binomial.1.1 = paste(Genus.1.1, Species.1.0, sep = "_")) %>% 
  filter(!str_detect(Binomial.1.2, "000"))

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

syn <- syn %>% 
  gather(System, Binomial, -Binomial.1.2) %>% 
  filter(Binomial.1.2 != Binomial) %>% 
  filter(!duplicated(Binomial)) %>% 
  bind_rows(mam %>% transmute(Binomial.1.2,
                              System = "Binomial.1.2",
                              Binomial = Binomial.1.2))
anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial"))

# Cercopithecus_pogonias
# Possibly: Cercopithecus denti
# http://www.iucnredlist.org/details/136885/0

# Cebus_olivaceus is split into two species:
# Cebus brunneus and Cebus kaapori - pick a random

extra <- data_frame(Binomial.1.2 = c("Cercopithecus_denti", "Cebus_brunneus"), 
                    System = "Personal.judgement", 
                    Binomial = c("Cercopithecus_pogonias", "Cebus_olivaceus"))
syn <- syn %>% bind_rows(extra)
stopifnot(!anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial")) %>% nrow)

pantheria <- left_join(pantheria, syn %>% select(-System), by = c("Binomial.Pantheria" = "Binomial"))

pantheria <- pantheria %>% left_join(mam)

dup <- pantheria %>% filter(duplicated(Binomial.1.2)) %>% pull(Binomial.1.2)
pantheria %>% filter(Binomial.1.2 %in% dup) %>% select(Binomial.1.2, Binomial.Pantheria)

# Remove duplicated species and keep nominal
pantheria <- pantheria %>% 
  filter(!Binomial.Pantheria %in% c("Damaliscus_korrigum",
                                    "Felis_catus",
                                    "Alcelaphus_caama",
                                    "Alcelaphus_lichtensteinii"))
stopifnot(!anyDuplicated(pantheria$Binomial.1.2))

pantheria <- pantheria %>% 
  mutate(log10density = log10(density),
         log10BM = log10(Mass.g)) %>%
  select(Binomial.1.2, Order.1.2, Family.1.2, log10density, log10BM, source)

bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")

terrestrial <- mam %>% filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
                              !Family.1.2 %in% c(whale.families, seal.families),
                              !Binomial.1.2 %in% marine.carnivores) %>% pull(Binomial.1.2)

pantheria <- pantheria %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "PanTHERIA")

tetra.dens <- read_csv("data/TetraDENSITY.csv")
tetra.dens <- tetra.dens %>% filter(Class == "Mammalia")
tetra.dens <- tetra.dens %>% mutate(Binomial = paste(Genus, Species, sep = "_"))

mismatch <- which(!tetra.dens$Binomial %in% mam$Binomial.1.2)
tetra.dens[mismatch,]$Binomial %>% unique

# To do:
# - Fix mismatches
# - Remove duplicates
# - Redo imputation if data looks sound?...

tetra.dens <- tetra.dens %>% left_join(mam, by = c("Binomial" = "Binomial.1.2"))
tetra.dens <- tetra.dens %>% filter(Binomial %in% terrestrial)
tetra.dens <- tetra.dens %>%
  mutate(Binomial.1.2 = Binomial,
         log10density = log10(Density),
         log10BM = log10(Mass.g),
         dataset = "TetraDENSITY") %>%
  select(Binomial.1.2, Order.1.2, Family.1.2, log10density, log10BM, dataset)

df <- bind_rows(pantheria, tetra.dens)

ggplot(df, aes(log10BM, log10density, col = dataset)) +
  geom_point() +
  geom_smooth(method = "lm")
