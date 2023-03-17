### Show family uncertainty ###


# Setup -------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ape)
library(ggtree)
# BiocManager::install("ggtreeExtra")
library(ggtreeExtra)

# Load forest
forest <- read_rds("builds/forest.rds")

# Subset to tree one
tree <- forest[[1]]

# Load PHYLACINE
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Load density
density <- read_csv("output/Table S4 Imputed density.csv")

# Read plotting info from metablic rate uncertainty plot
plot.info <- read_rds("../metabolic_rate/builds/uncertainty.plot.info.rds")

# Create family tree ------------------------------------------------------
# Replace all species names in tree with family names
index <- match(tree$tip.label, mam$Binomial.1.2)
tree$tip.label[] <- mam$Family.1.2[index]

# Subset to one tip per family
family.tree <- drop.tip(tree, which(duplicated(tree$tip.label)))

# Rearrange tree for plotting as metabolic rate plot
family.tree <- ape::rotateConstr(family.tree, plot.info$family.order)

# Plot tree ---------------------------------------------------------------
p.1 <- ggtree(family.tree, layout = "fan", ladderize = FALSE) +
  geom_tiplab(offset = 5)


# Highlight orders --------------------------------------------------------
node <- sapply(unique(density$Order.1.2), function(Order) getMRCA(family.tree, unique(density$Family.1.2[which(density$Order.1.2 == Order)]))) %>% unlist
node <- as_tibble(node, rownames = "Order.1.2")

order.of.orders <- unique(mam$Order.1.2[match(get_taxa_name(p.1), mam$Family.1.2)])
node$Order.1.2 <- as_factor(node$Order.1.2)
node$Order.1.2 <- fct_reorder(node$Order.1.2, match(node$Order.1.2, order.of.orders))

highlights <- scales::hue_pal()(length(plot.info$fill))
highlights <- highlights[sort(match(node$Order.1.2, plot.info$fill))]


p.2 <- p.1 +
  geom_hilight(data=node, aes(node=value, fill=Order.1.2), alpha = .5) +
  scale_fill_manual(values = highlights, guide = guide_legend(order = 2))


# Add point with color for mean uncertainty of FMR in family -------------
fam.density.sd <- density %>% group_by(Family.1.2) %>% summarise(fam.sd = mean(sd))

p.3 <- p.2  %<+% fam.density.sd +
  geom_tippoint(aes(color=fam.sd), size=3, position = position_nudge(x = 3)) +
  scale_color_continuous("Mean logDensity standard deviation", low = "black", high = "red")
p.3


# Save plot --------------------------------------------------------------
ggsave("output/appendix1_density_uncertainty.png", p.3, width = 35, height = 30, units = "cm")
