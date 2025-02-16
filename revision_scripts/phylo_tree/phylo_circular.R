# Load necessary libraries
library(ggtree)
library(viridis)
library(ape)
library(phytools)
library(tidytree)

# Read in the complete tree
full_tree <- read.tree("revision_scripts/phylo_tree/PhytoPhylo")

# Clean up species names
full_tree$tip.label <- gsub("_", " ", full_tree$tip.label)

# Create basic circular tree
circular_tree <- ggtree(full_tree, layout="circular") 

# Add all tips as small grey points
circular_tree <- circular_tree + 
  geom_tippoint(size=0.5, color="grey80")

# Add colored points and labels just for your subset
circular_tree <- circular_tree %<+% distance_data + 
  geom_tippoint(data=subset(circular_tree$data, label %in% distance_data$species),
                aes(color=label), size=3) +
  geom_tiplab(data=subset(circular_tree$data, label %in% distance_data$species),
              aes(color=label), size=5) +
  scale_color_manual(values=setNames(distance_data$color, distance_data$species)) +
  theme_tree2() +
  theme(legend.position="none")

# Display the plot
print(circular_tree)
