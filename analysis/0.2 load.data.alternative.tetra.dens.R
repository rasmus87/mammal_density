library(tidyverse)
library(ape)

# Load Pantheria
pantheria <- read_tsv("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>% na_if(-999)
names(pantheria) <- make.names(names(pantheria))
pantheria <- pantheria %>% 
  filter(!is.na(X21.1_PopulationDensity_n.km2)) %>% 
  transmute(Binomial.Pantheria = paste(MSW05_Genus, MSW05_Species, sep = "_"),
            density = X21.1_PopulationDensity_n.km2,
            source = "PanTHERIA.2008")

# Load tetra density
tetra <- read_csv("data/TetraDENSITY.csv")
tetra <- tetra %>% 
  filter(Class == "Mammalia") %>% 
  transmute(Binomial.Tetra = paste(Genus, Species, sep = "_"),
            density = Density,
            source = "TetraDENSITY.2018")


# Load taxonomy solver
syn <- read_csv("../PHYLACINE_1.2/Data/Taxonomy/Synonymy_table_with_unaccepted_species.csv", col_types = cols(), guess_max = 5000)
syn <- syn %>% 
  transmute(Binomial.1.2,
            Binomial.EltonTraits.1.0 = paste(EltonTraits.1.0.Genus, EltonTraits.1.0.Species, sep = "_"),
            Binomial.1.0 = paste(Genus.1.0, Species.1.0, sep = "_"),
            Binomial.1.1 = paste(Genus.1.1, Species.1.0, sep = "_")) %>% 
  filter(!str_detect(Binomial.1.2, "000"))

mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

syn <- syn %>% 
  gather(System, Binomial, -Binomial.1.2) %>% 
  filter(Binomial.1.2 != Binomial) %>% 
  filter(!duplicated(Binomial)) %>% 
  bind_rows(mam %>% transmute(Binomial.1.2,
                              System = "Binomial.1.2",
                              Binomial = Binomial.1.2))


# Align Pantheria taxonomy >>>
anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial"))

# Cercopithecus_pogonias
# Possibly: Cercopithecus denti
# http://www.iucnredlist.org/details/136885/0

# Cebus_olivaceus is split into two species:
# Cebus brunneus and Cebus kaapori - pick a random

extra <- tibble(Binomial.1.2 = c("Cercopithecus_denti", "Cebus_brunneus"), 
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
# Align Pantheria taxonomy |||


# Align TetraDENSITY taxonomy >>>
anti_join(tetra, syn, by = c("Binomial.Tetra" = "Binomial")) %>% filter(!duplicated(Binomial.Tetra))

# Cercopithecus_pogonias
# Possibly: Cercopithecus denti
# http://www.iucnredlist.org/details/136885/0

# Cebus_olivaceus is split into two species:
# Cebus brunneus and Cebus kaapori - pick a random

# Gazella_arabica
# This is one of the two species that Gazella_gazella was split into

# Capricornis_sumatrensis 
# Misspelling of Capricornis_sumatraensis 

# Ovis musimon
# is Ovis aries

extra <- tibble(Binomial.1.2 = c("Capricornis_sumatraensis",
                                 "Gazella_gazella",
                                 "Ovis_aries",
                                 "Macropus_rufogriseus",
                                 "Cebus_brunneus",
                                 "Cercopithecus_denti",
                                 "Cercopithecus_denti",
                                 "Hylobates_mulleri",
                                 "Prosciurillus_rosenbergii"), 
                System = "Personal.judgement", 
                Binomial = c("Capricornis_sumatrensis",
                             "Gazella_arabica",
                             "Ovis_musimon",
                             "Notamacropus_rufogriseus",
                             "Cebus_olivaceus",
                             "Cercopithecus_pogonias",
                             "Cercopithecus_wolfi",
                             "Hylobates_mulleri x agilis",
                             "Prosciurillus_rosenbergi"))
syn <- syn %>% bind_rows(extra)
stopifnot(!anti_join(tetra, syn, by = c("Binomial.Tetra" = "Binomial")) %>% nrow)

tetra <- left_join(tetra, syn %>% select(-System), by = c("Binomial.Tetra" = "Binomial"))

tetra <- tetra %>% left_join(mam)

dup <- tetra %>% filter(duplicated(Binomial.1.2)) %>% pull(Binomial.1.2)
tetra %>% filter(Binomial.1.2 %in% dup) %>% select(Binomial.1.2, Binomial.Tetra)

tetra <- tetra %>% 
  mutate(log10density = log10(density),
         log10BM = log10(Mass.g)) %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10density, log10BM, source)
# Align TetraDENSITY taxonomy |||


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

density.dataset <- pantheria %>% bind_rows(tetra) %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "density")
write_csv(density.dataset, "builds/imputation_dataset_PanTetra.csv")

density.dataset <- as.data.frame(density.dataset)

mam <- mam %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  mutate(log10BM = log10(Mass.g), log10density = NA)
mam <- mam %>% mutate(dataset = "mam")
mam <- mam %>% select(names(density.dataset))
mam <- as.data.frame(mam) # Very important for MCMCglmm

n.mam <- nrow(mam)

df <- rbind(mam, density.dataset)
df <- as.data.frame(df)

forest <- readRDS("builds/forest.rds")
species <- forest[[1]]$tip.label
drop.species <- species[!species %in% terrestrial]
forest <- lapply(forest, drop.tip, tip = drop.species)
