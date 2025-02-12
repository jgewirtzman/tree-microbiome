source("revision_scripts/phylo/S.phylomaker.R")


S.PhyloMaker(tree, spList, nodes, output.spList = T, scenarios = c("S1", "S2", "S3")) 
  


library("phytools")                       # load the "phytools" package.    
#example<-read.csv("example.splist",header=T, sep="\t")       # read in the example species list.    
my_splist<-read.csv("revision_scripts/phylo/updated_species_data.csv")       # read in the example species list.    
phylo<-read.tree("revision_scripts/phylo/PhytoPhylo")      # read in the megaphylogeny.    
nodes<-read.csv("revision_scripts/phylo/nodes",header=T, sep="\t")     # read in the nodes information of the megaphylogeny.    
#result<-S.PhyloMaker(spList=example, tree=phylo, nodes=nodes)      # run the function S.PhyloMaker.    
result<-S.PhyloMaker(spList=my_splist, tree=phylo, nodes=nodes)      # run the funxction S.PhyloMaker.    
str(result)       # the structure of the ouput of S.PhyloMaker.    
par(mfrow=c(3,1),mar=c(0,0,1,0))       # show the phylogenies of the three scenarios.    
plot(result$Scenario.1,cex=1.1,main="Scenarion One")    
plot(result$Scenario.2,cex=1.1,main="Scenarion Two")    
plot(result$Scenario.3,cex=1.1,main="Scenarion Three")  
dev.off()

# Load necessary libraries
library(ape)
library(phytools)
library(vegan)

# Assume each tree is loaded from the result object
tree_scenario1 <- result$Scenario.1

# Calculate cophenetic distances (phylogenetic distances) for each scenario
dist_matrix1 <- cophenetic(tree_scenario1)


# Load necessary libraries
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("umap", quietly = TRUE)) install.packages("umap")
if (!requireNamespace("Rtsne", quietly = TRUE)) install.packages("Rtsne")

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("umap", quietly = TRUE)) install.packages("umap")
if (!requireNamespace("Rtsne", quietly = TRUE)) install.packages("Rtsne")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

library(ggplot2)
library(vegan)
library(umap)
library(Rtsne)
library(ggrepel)
library(cowplot)  # For combining plots

# Prepare the distance matrix as a matrix (if it is not already)
dist_matrix1 <- as.matrix(dist_matrix1)

# 1. PCA
pca_result <- prcomp(dist_matrix1, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x[, 1:2])  # Take the first two principal components
pca_data$species <- rownames(dist_matrix1)  # Add species names for labeling

# 2. NMDS
nmds_result <- metaMDS(dist_matrix1, k = 2)
nmds_data <- as.data.frame(nmds_result$points)
nmds_data$species <- rownames(dist_matrix1)

# 3. UMAP
umap_result <- umap(dist_matrix1, n_neighbors = 15, min_dist = 0.1, n_components = 2)
umap_data <- as.data.frame(umap_result$layout)
umap_data$species <- rownames(dist_matrix1)

# 4. t-SNE
set.seed(123)  # For reproducibility
tsne_result <- Rtsne(dist_matrix1, dims = 2, perplexity = 5)
tsne_data <- as.data.frame(tsne_result$Y)
tsne_data$species <- rownames(dist_matrix1)


# Calculate average phylogenetic distance for each species
average_distances <- rowMeans(dist_matrix1)  # Mean distance to other species for each species

# Add this distance as a color variable to each data frame
pca_data$distance <- average_distances
nmds_data$distance <- average_distances
umap_data$distance <- average_distances
tsne_data$distance <- average_distances


library(ggrepel)
library(ggplot2)
library(cowplot)

# 1. PCA Plot with ggrepel and color gradient
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = species, color = log(distance))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = species), size = 3, max.overlaps = Inf) +
  labs(title = "PCA", x = "PC1", y = "PC2") +
  scale_color_gradient(low = "blue", high = "red", name = "Phylogenetic\nDistance") +
  theme_minimal()

# 2. NMDS Plot with ggrepel and color gradient
nmds_plot <- ggplot(nmds_data, aes(x = MDS1, y = MDS2, label = species, color = log(distance))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = species), size = 3, max.overlaps = Inf) +
  labs(title = "NMDS", x = "NMDS1", y = "NMDS2") +
  scale_color_gradient(low = "blue", high = "red", name = "Phylogenetic\nDistance") +
  theme_minimal()

# 3. UMAP Plot with ggrepel and color gradient
umap_plot <- ggplot(umap_data, aes(x = V1, y = V2, label = species, color = log(distance))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = species), size = 3, max.overlaps = Inf) +
  labs(title = "UMAP", x = "UMAP1", y = "UMAP2") +
  scale_color_gradient(low = "blue", high = "red", name = "Phylogenetic\nDistance") +
  theme_minimal()

# 4. t-SNE Plot with ggrepel and color gradient
tsne_plot <- ggplot(tsne_data, aes(x = V1, y = V2, label = species, color = log(distance))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = species), size = 3, max.overlaps = Inf) +
  labs(title = "t-SNE", x = "t-SNE1", y = "t-SNE2") +
  scale_color_gradient(low = "blue", high = "red", name = "Phylogenetic\nDistance") +
  theme_minimal()

# Combine plots into a 2x2 grid
combined_plot <- plot_grid(pca_plot, nmds_plot, umap_plot, tsne_plot, 
                           labels = c("A", "B", "C", "D"), 
                           ncol = 2, nrow = 2)

# Display the combined plot
print(combined_plot)

# Assuming tree_scenario1 is your tree and average_distances contains average distances
library(ggtree)
library(viridis)

# Convert the average_distances to a data frame and ensure species names match
distance_data <- data.frame(species = names(average_distances), distance = average_distances)

# Check that the species names in `distance_data` match the tree's tip labels
tree_scenario1$tip.label <- gsub(" ", "_", tree_scenario1$tip.label)  # Ensure matching names


library(ggtree)

# Use ggtree to plot, add color based on distance, and add colored points at each tip
tree_plot <- ggtree(tree_scenario1) %<+% distance_data +
  geom_tiplab(aes(color = log10(distance)), fontface = "bold") +  # Make labels bold
  geom_tippoint(aes(color = log10(distance)), size = 2) +  # Add colored points at each tip
  scale_color_viridis_c(name = "Log10 Phylogenetic Distance") +
  theme_tree2() +
  #labs(title = "Phylogenetic Tree Colored by Average Distance") +
  theme(legend.position = "top") +
  ggplot2::xlim(0, 500)

# Display the plot
print(tree_plot)


# Load necessary libraries
library(ggtree)
library(viridis)
library(ape)
library(phytools)

# Function to create a gradient color palette with viridis
create_viridis_palette <- function(data) {
  # Remove underscores from species names
  names(data) <- gsub("_", " ", names(data))
  
  # Convert data to a data frame
  df <- data.frame(Species = names(data), Value = as.numeric(data))
  
  # Generate a viridis color palette
  color_palette <- setNames(viridis(length(data)), df$Species[order(df$Value)])
  
  # Return the named color palette
  return(color_palette)
}

# Define species data with phylogenetic distances
data <- distance_data[,1:2]

# Generate the viridis color palette for the species
color_palette_phylo <- create_viridis_palette(data)

# Read in the phylogenetic tree
tree_scenario1 <- read.tree("revision_scripts/phylo/PhytoPhylo")  

# Convert average distances to a data frame
distance_data <- data.frame(species = names(data), distance = data)

# Ensure species names match in both datasets (convert underscores to spaces)
distance_data$species <- gsub("_", " ", distance_data$species)
tree_scenario1$tip.label <- gsub("_", " ", tree_scenario1$tip.label)

# Assign viridis colors based on species names
distance_data$color <- color_palette_phylo[distance_data$species]

# Check that colors are correctly assigned
print(distance_data)

# Plot the phylogenetic tree using ggtree with exact viridis colors
tree_plot <- ggtree(tree_scenario1) %<+% distance_data +
  geom_tiplab(aes(color = species), fontface = "bold") +  # Label species with color
  geom_tippoint(aes(color = species), size = 2) +  # Add colored points
  scale_color_manual(values = color_palette_phylo, name = "Phylogenetic Distance") +
  theme_tree2() +
  theme(legend.position = "top") +
  ggplot2::xlim(0, 500)  # Adjust the x-axis if needed

# Display the plot
print(tree_plot)
