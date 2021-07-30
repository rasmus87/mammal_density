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


# Make terrestrial filter -------------------------------------------------

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
# (There shouldn't be duplicates in PanTHERIA, and duplicates are because of species misalignment)
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

# Filter to terrestrial species only
pantheria <- pantheria %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "density")

# Save dataset for imputation
write_csv(pantheria, "builds/imputation_dataset_PanTHERIA.csv")


# Load and align TetraDENSITY ---------------------------------------------

# Load tetra density
tetra <- read_csv("data/TetraDENSITY.csv") %>% 
  filter(Class == "Mammalia") %>% 
  transmute(Binomial.Tetra = paste(Genus, Species, sep = "_"),
            density = Density,
            source = "TetraDENSITY.2018")

# Check for any mismatches in synonomy
anti_join(tetra, syn, by = c("Binomial.Tetra" = "Binomial")) %>% 
  filter(!duplicated(Binomial.Tetra))

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

# Expand synonomy table to include these
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

# Check that we are fully aligned
stopifnot(!anti_join(tetra, syn, by = c("Binomial.Tetra" = "Binomial")) %>% nrow)

# Join synonomy and PHYLACINE traits to TetraDENSITY
tetra <- tetra %>% 
  left_join(syn %>%
              select(-System),
            by = c("Binomial.Tetra" = "Binomial")) %>% 
  left_join(mam)

# Which species are duplicates?
# Duplicates are acceptable here, since TeraDENSITY has multiple records of many species
dup <- tetra %>% filter(duplicated(Binomial.1.2)) %>% pull(Binomial.1.2)
tetra %>% filter(Binomial.1.2 %in% dup) %>% select(Binomial.1.2, Binomial.Tetra)

# Reduce dataset columns
tetra <- tetra %>% 
  mutate(log10density = log10(density),
         log10BM = log10(Mass.g)) %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10density, log10BM, source)

# Filter to terrestrial species only
tetra <- tetra %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "density")

# Combine datasets
density.dataset.alt <- pantheria %>% 
  bind_rows(tetra) %>% 
  mutate(dataset = "density.alt")
write_csv(density.dataset.alt, "builds/imputation_dataset_PanTHERIA_TetraDENSITY.csv")



# Reduce and prepare the traits dataset mam for imputation ----------------

# Prepare mam for imputation
mam <- mam %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  mutate(log10BM = log10(Mass.g), 
         log10density = NA,
         dataset = "mam") %>% 
  select(names(pantheria))
write_csv(mam, "builds/imputation_dataset_mam.csv")


# Prepare things for imputation -------------------------------------------
# 
# # Turn pantheria into data.frame for MCMCglmm
# pantheria <- as.data.frame(pantheria)
# 
# # Turn pantheria into data.frame for MCMCglmm
# density.dataset.alt <- as.data.frame(density.dataset.alt)
# 
# # Turn mam into data.frame for MCMCglmm
# mam <- as.data.frame(mam)
# 
# # Number of species
# n.mam <- nrow(mam)
# 
# # Combine datasets
# df <- rbind(mam, pantheria)
# df.alt <- rbind(mam, density.dataset.alt)