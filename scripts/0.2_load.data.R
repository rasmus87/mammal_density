# Load mammal density data
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(ape)


# Make taxonomy solver ----------------------------------------------------

# Load PHYLACINE 1.2.1 trait data
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Load PHYLACINE 1.2.1 large synonomy
syn <- read_csv("../PHYLACINE_1.2/Data/Taxonomy/Synonymy_table_with_unaccepted_species.csv")
# Simplify syn table
syn <- syn %>% 
  transmute(Binomial.1.2,
            Binomial.EltonTraits.1.0 = paste(EltonTraits.1.0.Genus, EltonTraits.1.0.Species, sep = "_"),
            Binomial.1.0 = paste(Genus.1.0, Species.1.0, sep = "_"),
            Binomial.1.1 = paste(Genus.1.1, Species.1.0, sep = "_")) %>% 
  filter(!str_detect(Binomial.1.2, "000"))

# Make full synonomy matching table
syn <- syn %>% 
  pivot_longer(names_to = "System", values_to = "Binomial", -Binomial.1.2) %>% 
  filter(Binomial.1.2 != Binomial) %>% 
  filter(!duplicated(Binomial)) %>% 
  bind_rows(mam %>% transmute(Binomial.1.2,
                              System = "Binomial.1.2",
                              Binomial = Binomial.1.2))



# Load and align PanTHERIA ------------------------------------------------

# Load PanTHERIA
pantheria <- read_tsv("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>% 
  na_if(-999) %>% 
  filter(!is.na(`21-1_PopulationDensity_n/km2`)) %>% 
  transmute(Binomial.Pantheria = paste(MSW05_Genus, MSW05_Species, sep = "_"),
            density = `21-1_PopulationDensity_n/km2`,
            source = "PanTHERIA.2008")

# Check for any mismatches in synonomy
anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial"))

# Cercopithecus_pogonias
# Possibly: Cercopithecus denti
# http://www.iucnredlist.org/details/136885/0

# Cebus_olivaceus is split into two species:
# Cebus brunneus and Cebus kaapori - pick a random

# Expand synonomy table to include these
extra <- tibble(Binomial.1.2 = c("Cercopithecus_denti", "Cebus_brunneus"), 
                System = "Author.judgement", 
                Binomial = c("Cercopithecus_pogonias", "Cebus_olivaceus"))
syn <- syn %>% bind_rows(extra)

# Check that we are fully aligned
stopifnot(!anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial")) %>% nrow)

# Join synonomy and PHYLACINE traits to PanTHERIA
pantheria <- pantheria %>% 
  left_join(syn %>% 
              select(-System), 
            by = c("Binomial.Pantheria" = "Binomial")) %>% 
  left_join(mam, by = "Binomial.1.2")

# Which species are duplicates?
dup <- pantheria %>% filter(duplicated(Binomial.1.2)) %>% pull(Binomial.1.2)
pantheria %>% filter(Binomial.1.2 %in% dup) %>% select(Binomial.1.2, Binomial.Pantheria)

# Remove duplicated species and keep nominal
pantheria <- pantheria %>% 
  filter(!Binomial.Pantheria %in% c("Damaliscus_korrigum",
                                    "Felis_catus",
                                    "Alcelaphus_caama",
                                    "Alcelaphus_lichtensteinii"))
# Make sure we don't have any further duplicates
stopifnot(!anyDuplicated(pantheria$Binomial.1.2))

# Reduce dataset columns
pantheria <- pantheria %>% 
  mutate(log10density = log10(density),
         log10BM = log10(Mass.g)) %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10density, log10BM, source)

# Non terrestrials and humans
bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")
humans <- "Homo"

terrestrial <- mam %>% 
  filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
         !Family.1.2 %in% c(whale.families, seal.families),
         !Binomial.1.2 %in% marine.carnivores,
         !Genus.1.2 %in% humans) %>% 
  pull(Binomial.1.2)

# Filter to terrestrial species only
pantheria <- pantheria %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "density")

# Save dataset for imputation
write_csv(pantheria, "builds/imputation_dataset.csv")

# Turn pantheria into data.frame for MCMCglmm
pantheria <- as.data.frame(pantheria)

# Prepare mam for imputation
mam <- mam %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  mutate(log10BM = log10(Mass.g), 
         log10density = NA,
         dataset = "mam") %>% 
  select(names(pantheria))
# Turn mam into data.frame for MCMCglmm
mam <- as.data.frame(mam)

# Number of species
n.mam <- nrow(mam)

# Combine datasets
df <- rbind(mam, pantheria)

# Reduce forest to the designated species list
forest <- readRDS("builds/forest.rds")
species <- forest[[1]]$tip.label
drop.species <- species[!species %in% terrestrial]
forest <- lapply(forest, drop.tip, tip = drop.species)
