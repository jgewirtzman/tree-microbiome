library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library("RColorBrewer")
library(viridis)
library("phytools")
library(patchwork)
library(cowplot)

data<-GC_data
data <- data %>% filter(nchar(Species.ID) == 4) %>%
  group_by(Species.ID) %>%
  mutate(n_obs = n()) %>%
  ungroup()


#example<-read.csv("example.splist",header=T, sep="\t")       # read in the example species list.    

source("S.PhyloMaker-master/S.PhyloMaker.R")

my_splist<-read.csv("S.PhyloMaker-master/updated_species_data.csv")       # read in the example species list.    
phylo<-read.tree("S.PhyloMaker-master/PhytoPhylo")      # read in the megaphylogeny.    
nodes<-read.csv("S.PhyloMaker-master/nodes",header=T, sep="\t")     # read in the nodes information of the megaphylogeny.    
#result<-S.PhyloMaker(spList=example, tree=phylo, nodes=nodes)      # run the function S.PhyloMaker.    
result<-S.PhyloMaker(spList=my_splist, tree=phylo, nodes=nodes)      # run the function S.PhyloMaker.    

# Assume each tree is loaded from the result object
tree_scenario1 <- result$Scenario.1

# Calculate cophenetic distances (phylogenetic distances) for each scenario
dist_matrix1 <- cophenetic(tree_scenario1)

# Prepare the distance matrix as a matrix (if it is not already)
dist_matrix1 <- as.matrix(dist_matrix1)


# Calculate average phylogenetic distance for each species
average_distances <- rowMeans(dist_matrix1)  # Mean distance to other species for each species

average_distances <- data.frame(
  Species = names(average_distances),
  Distance = as.numeric(average_distances)
)

average_distances <- average_distances %>%
  mutate(
    Species.ID = toupper(paste0(substr(Species, 1, 2), substr(gsub(".*_", "", Species), 1, 2))),
    Species = gsub("_", " ", Species) # Remove underscores from species names
  )

# Assuming `data` is already loaded, merge `average_distances` with `data` on Species.ID
merged_data <- data %>%
  left_join(average_distances, by = "Species.ID")

# Sort average_distances by Distance and Species to ensure a unique order
ordered_species <- average_distances %>%
  arrange(Distance, Species) %>%
  pull(Species)

# Set Species factor levels in merged_data to follow this order
merged_data <- merged_data %>%
  mutate(Species = factor(Species, levels = ordered_species))

# Check the structure to ensure no NA values
str(merged_data)

# Generate a color palette based on the number of species
species_colors <- viridis::viridis(n = nrow(merged_data))

# Map colors to each species based on distance
names(species_colors) <- merged_data$Species

# Add a Color column to average_distances by matching Species with species_colors
average_distances$Color <- species_colors[average_distances$Species]

# Create a named vector for the species palette
# Generate a color palette based on the unique number of species
species_palette <- viridis::viridis(n = length(unique(merged_data$Species)))
names(species_palette) <- levels(merged_data$Species)  # Assign species names to colors
species_palette <- species_palette[levels(merged_data$Species)]

# Plot each individual chart with Species legend color-coded by distance, and remove Distance color bar
pa <- ggplot(merged_data, aes(x = O2_concentration / 1000000 * 100, y = CH4_concentration / 1000000 * 100)) +
  geom_point(aes(color = Species), alpha = 0.7, show.legend = TRUE) +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab(expression("O"[2] ~ "(%)"))+
  ylab(expression("CH"[4] ~ "(%)"))+ 
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  #ylim(0,7)+
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  scale_color_manual(values = species_palette, guide = guide_legend(title = NULL, nrow = 4, order = 1))

pb <- ggplot(merged_data, aes(x = O2_concentration / 1000000 * 100, y = N2O_concentration)) +
  geom_point(aes(color = Species), alpha = 0.7, show.legend = TRUE) +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab(expression("O"[2] ~ "(%)"))+
  ylab(expression("N"[2]*"O (%)"))+ 
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) + 
  ylim(0,3)+
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  scale_color_manual(values = species_palette, guide = guide_legend(title = NULL, nrow = 4, order = 1))

pc <- ggplot(merged_data, aes(x = O2_concentration / 1000000 * 100, y = CO2_concentration / 1000000 * 100)) +
  geom_point(aes(color = Species), alpha = 0.7, show.legend = TRUE) +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  xlab(expression("O"[2] ~ "(%)"))+
  ylab(expression("CO"[2] ~ "(%)"))+ 
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  scale_color_manual(values = species_palette, guide = guide_legend(title = NULL, nrow = 4, order = 1))

# Combine the plots with patchwork and keep only the Species legend
combined_plot <- (pc | pa | pb) + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom", 
    legend.box = "horizontal", 
    legend.box.just = "left",
    legend.key.size = unit(0.2, "cm"),               # Smaller legend key size for compactness
    legend.key.width = unit(0.15, "cm"),             # Further compress horizontally
    legend.spacing.x = unit(0.02, "cm"),             # Minimal horizontal spacing between items
    legend.spacing.y = unit(0.05, "cm"),             # Minimal vertical spacing for multi-row legends
    legend.text = element_text(size = 10),           # Slightly larger text for readability
    legend.title = element_text(size = 0),          # Adjust title size to balance with text
    legend.margin = margin(t = 2, r = 5, b = 2, l = 0, unit = "pt"),  # Minimal padding around legend box
    plot.margin = margin(1,1,1,1, unit = "pt")       # Minimal outer plot margin
  )

combined_plot

# Add the label "j" in the top left
ic_plot <- plot_grid(combined_plot, labels = "j", label_size = 14, label_x = 0.00, label_y = 0.95, hjust = 0, vjust = 1)

# Display the plot
ic_plot

