# Load required libraries
# Load or install necessary packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

# Load required packages
p_load(
  phyloseq, ape, dplyr, purrr, tidyr, ggplot2, multcompView, 
  RColorBrewer, pheatmap, emmeans, lme4, multcomp, lme4
)


#### 16S IMPORT #### 

#### Load OTU table ####
otu_tab <- read.delim("16S/OTU_table.txt", header=TRUE, row.names = 1, stringsAsFactors = FALSE)

# Extract taxonomy from last columns of OTU table
bastard_tax <- otu_tab[, 590:596]  # Assumes taxonomy is in last columns
bastard_tax[bastard_tax == ""] <- NA  # Replace empty strings with NA

# Convert taxonomy table into Phyloseq-compatible format
tax_tab_pre <- tax_table(as.matrix(bastard_tax))
taxa_names(tax_tab_pre) <- sub("sp", "seq", taxa_names(tax_tab_pre))  # Fix name formatting

# Remove taxonomy columns from OTU table to retain abundance matrix only
otu_tab_corr <- otu_tab[, 1:589]  # Keep only abundance data
otu_table_pre <- otu_table(as.matrix(otu_tab_corr), taxa_are_rows = TRUE)

#### Load Phylogenetic Tree ####
phylo_tree <- read_tree("16S/rooted_tree.nwk")

#### Load Sample Metadata ####
samp_data <- read.delim("16S/tree_16s_mapping_dada2_corrected.txt.no_gz", row.names = 1, stringsAsFactors = FALSE)
samp_data$RowName <- row.names(samp_data)

#### Merge Additional Metadata ####
ddpcr <- read.csv(file = "Other-Metadata/ddPCR_meta_all_data.csv", stringsAsFactors = FALSE)
metadata <- read.csv("Other-Metadata/annotated_metadata/16S_tree_sample_table_with_meta.csv", stringsAsFactors = FALSE)

#### Metadata Preparation and Filtering ####
metadata_filtered <- data.frame(
  SampleID   = metadata$X.1,
  Species    = substr(metadata$Inner.Core.Sample.ID.x, 1, 4),
  SampleType = metadata$core_type
)
# Filter for Inner/Outer samples
metadata_filtered <- subset(metadata_filtered, SampleType %in% c("Inner", "Outer"))
rownames(metadata_filtered) <- metadata_filtered$SampleID

# Merge metadata
samp_data_merged <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)
samp_data_merged <- merge(samp_data_merged, metadata_filtered, by.x = "RowName", by.y = "SampleID", all.x = TRUE)

# Remove duplicates
dups <- which(duplicated(samp_data_merged$RowName))
samp_data_merged <- samp_data_merged[-c(dups), ]
row.names(samp_data_merged) <- samp_data_merged$RowName

#### Construct Phyloseq Object ####
raw_ps <- phyloseq(otu_table_pre, tax_tab_pre, phylo_tree, sample_data(samp_data_merged))

#### Remove Plastids ####
pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

mitochondria <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 5] == "Mitochondria")]
chloroplast <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 4] == "Chloroplast")]
badTaxa <- c(mitochondria, chloroplast)

no_mito <- pop_taxa(raw_ps, badTaxa)
taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))

#### Begin Filtering Pipeline ####
# Initial state
print(paste("Initial counts - Samples:", nsamples(no_mito), "ASVs:", ntaxa(no_mito)))

# 1. Sample-level filtering
# Second Inner/Outer filter at phyloseq level
ps.filt <- prune_samples(sample_data(no_mito)$core_type %in% c("Inner", "Outer"), no_mito)
print(paste("After Inner/Outer filtering - Samples:", nsamples(ps.filt), "ASVs:", ntaxa(ps.filt)))

# Remove samples with insufficient reads
min_reads <- 1000
ps.filt <- prune_samples(sample_sums(ps.filt) >= min_reads, ps.filt)
print(paste("After read depth filtering - Samples:", nsamples(ps.filt), "ASVs:", ntaxa(ps.filt)))

# 2. Taxa-level filtering
# Remove zero-abundance taxa
ps.filt <- prune_taxa(taxa_sums(ps.filt) > 0, ps.filt)
print(paste("After removing zero-abundance taxa - Samples:", nsamples(ps.filt), "ASVs:", ntaxa(ps.filt)))

# Remove low abundance taxa (at least 3 reads total)
ps.filt <- prune_taxa(taxa_sums(ps.filt) >= 3, ps.filt)
print(paste("After abundance filtering - Samples:", nsamples(ps.filt), "ASVs:", ntaxa(ps.filt)))

# Remove low prevalence taxa (present in at least 2 samples)
ps.filt <- prune_taxa(rowSums(otu_table(ps.filt) > 0) >= 2, ps.filt)
print(paste("After prevalence filtering - Samples:", nsamples(ps.filt), "ASVs:", ntaxa(ps.filt)))

# 3. Perform rarefaction
set.seed(46814)
ps.rare <- rarefy_even_depth(ps.filt, sample.size = 3500)
print(paste("After rarefaction - Samples:", nsamples(ps.rare), "ASVs:", ntaxa(ps.rare)))

# 4. Final tree pruning
remaining_taxa <- taxa_names(ps.rare)
ps.rare <- prune_taxa(remaining_taxa, ps.rare)
pruned_tree <- prune_taxa(remaining_taxa, phy_tree(ps.rare))
phy_tree(ps.rare) <- pruned_tree

# Store final object
phyloseq_clean <- ps.rare # The final cleaned phyloseq object

phyloseq_clean

###

#### Calculate Distance Matrices ####
wunifrac_dist_clean_16s <- phyloseq::distance(phyloseq_clean, method = "wunifrac")
range(wunifrac_dist_clean_16s); mean(wunifrac_dist_clean_16s)

uunifrac_dist_clean_16s <- phyloseq::distance(phyloseq_clean, method = "unifrac")
range(uunifrac_dist_clean_16s); mean(uunifrac_dist_clean_16s)

bray_dist_clean_16s <- phyloseq::distance(phyloseq_clean, method = "bray")
range(bray_dist_clean_16s); mean(bray_dist_clean_16s)

###

metadata_clean <- data.frame(sample_data(phyloseq_clean)) %>%
  filter(!is.na(Species) & !is.na(SampleType))

metadata_clean$Group <- paste(metadata_clean$Species, metadata_clean$SampleType, sep = "_")
group_names <- unique(metadata_clean$Group)

unifrac_matrix <- as.matrix(wunifrac_dist_clean_16s)

# List of distance matrices
distance_matrices <- list(
  wunifrac = as.matrix(wunifrac_dist_clean_16s),
  uunifrac = as.matrix(uunifrac_dist_clean_16s),
  braycurtis = as.matrix(bray_dist_clean_16s)
)

# Function to process each distance matrix and generate heatmap & barplot
# List of distance matrices
distance_matrices <- list(
  wunifrac = as.matrix(wunifrac_dist_clean_16s),
  uunifrac = as.matrix(uunifrac_dist_clean_16s),
  braycurtis = as.matrix(bray_dist_clean_16s)
)

# Function to process each distance matrix and generate heatmap & barplot
analyze_distance_matrix <- function(dist_matrix, metric_name) {
  # Define comparison order
  comparison_order <- c(
    "Sapwood_Heartwood_Within",
    "Sapwood_Sapwood_Within",
    "Heartwood_Heartwood_Within",
    "Sapwood_Heartwood_Between",
    "Sapwood_Sapwood_Between",
    "Heartwood_Heartwood_Between"
  )
  
  # Convert distance matrix to long format
  dist_df <- as.data.frame(as.matrix(dist_matrix))
  dist_df$Sample1 <- rownames(dist_df)
  
  dist_long <- dist_df %>%
    pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance") %>%
    filter(Sample1 != Sample2) %>%
    left_join(metadata_clean, by = c("Sample1" = "RowName")) %>%
    rename(Species1 = Species, SampleType1 = SampleType) %>%
    left_join(metadata_clean, by = c("Sample2" = "RowName")) %>%
    rename(Species2 = Species, SampleType2 = SampleType) %>%
    mutate(
      Sample1 = factor(Sample1),
      Sample2 = factor(Sample2),
      Comparison = factor(case_when(
        Species1 == Species2 & SampleType1 == "Outer" & SampleType2 == "Inner" ~ "Sapwood_Heartwood_Within",
        Species1 == Species2 & SampleType1 == "Outer" & SampleType2 == "Outer" ~ "Sapwood_Sapwood_Within",
        Species1 == Species2 & SampleType1 == "Inner" & SampleType2 == "Inner" ~ "Heartwood_Heartwood_Within",
        Species1 != Species2 & SampleType1 == "Outer" & SampleType2 == "Inner" ~ "Sapwood_Heartwood_Between",
        Species1 != Species2 & SampleType1 == "Outer" & SampleType2 == "Outer" ~ "Sapwood_Sapwood_Between",
        Species1 != Species2 & SampleType1 == "Inner" & SampleType2 == "Inner" ~ "Heartwood_Heartwood_Between",
        TRUE ~ NA_character_
      ), levels = comparison_order)
    ) %>%
    filter(!is.na(Comparison))
  
  # Fit mixed model
  mixed_model <- lmer(Distance ~ Comparison + (1|Sample1) + (1|Sample2),
                      data = dist_long)
  print(summary(mixed_model))
  
  # Calculate EMMs
  emm <- emmeans(mixed_model, specs = "Comparison")
  
  # Generate compact letter display
  cld <- multcomp::cld(emm, 
                       alpha = 0.05,
                       adjust = "tukey",
                       Letters = letters[1:6])
  
  # Prepare plot data
  plot_data <- as.data.frame(cld) %>%
    dplyr::select(Comparison, emmean, SE, .group) %>%
    rename(Letters = .group)
  
  # Define comparison labels
  comparison_labels <- c(
    "Sapwood_Heartwood_Within" = "Sapwood-\nHeartwood",
    "Sapwood_Sapwood_Within" = "Sapwood-\nSapwood",
    "Heartwood_Heartwood_Within" = "Heartwood-\nHeartwood",
    "Sapwood_Heartwood_Between" = "Sapwood-\nHeartwood",
    "Sapwood_Sapwood_Between" = "Sapwood-\nSapwood",
    "Heartwood_Heartwood_Between" = "Heartwood-\nHeartwood"
  )
  
  # Calculate y-axis limits
  y_max <- max(plot_data$emmean + plot_data$SE, na.rm = TRUE)
  bracket_height <- y_max * 0.1
  text_height <- bracket_height * 1.2
  
  # Generate barplot
  barplot_plot <- ggplot(plot_data, aes(x = Comparison, y = emmean, fill = Comparison)) +
    geom_bar(stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
    geom_text(aes(y = emmean + SE + y_max * 0.05, label = Letters), size = 5, hjust = 0.5) +
    
    # Bracket for "Within Species"
    geom_segment(aes(x = 1, xend = 3, 
                     y = y_max + bracket_height, 
                     yend = y_max + bracket_height), size = 0.5) +
    annotate("text", x = 2, 
             y = y_max + text_height, 
             label = "Within Species", size = 4) +
    
    # Bracket for "Between Species"
    geom_segment(aes(x = 4, xend = 6, 
                     y = y_max + bracket_height, 
                     yend = y_max + bracket_height), size = 0.5) +
    annotate("text", x = 5, 
             y = y_max + text_height, 
             label = "Between Species", size = 4) +
    
    scale_x_discrete(labels = comparison_labels) +
    scale_fill_manual(values = c(
      "Sapwood_Heartwood_Within" = "grey50",
      "Sapwood_Sapwood_Within" = "#377EB8",
      "Heartwood_Heartwood_Within" = "#8B4513",
      "Sapwood_Heartwood_Between" = "grey50",
      "Sapwood_Sapwood_Between" = "#377EB8",
      "Heartwood_Heartwood_Between" = "#8B4513"
    )) +
    labs(x = "Comparison", 
         y = switch(metric_name,
                    "wunifrac" = "WUnifrac Distance",
                    "uunifrac" = "Unweighted Unifrac Distance",
                    "braycurtis" = "Bray-Curtis Distance")) +
    theme_classic() +
    theme(legend.position = "none")
  
  # Store the plot in the global environment
  assign(paste0("barplot_", metric_name), barplot_plot, envir = .GlobalEnv)
  
  # Print the plot
  print(barplot_plot)
}

# Run analysis for each distance metric
for (metric in names(distance_matrices)) {
  analyze_distance_matrix(distance_matrices[[metric]], metric)
}

# If you want to combine plots
library(patchwork)
combined_plot <- barplot_wunifrac + barplot_uunifrac
print(combined_plot)




generate_heatmap <- function(dist_matrix, metric_name) {
  
  metadata_clean <- data.frame(sample_data(phyloseq_clean)) %>%
    filter(!is.na(Species) & !is.na(SampleType))
  
  metadata_clean$Group <- paste(metadata_clean$Species, metadata_clean$SampleType, sep = "_")
  group_names <- unique(metadata_clean$Group)
  
  # Initialize distance matrix
  group_distance_matrix <- matrix(
    0,
    nrow = length(group_names),
    ncol = length(group_names),
    dimnames = list(group_names, group_names)
  )
  
  # Calculate mean distances for groups
  for (i in seq_along(group_names)) {
    for (j in seq_along(group_names)) {
      group1_samples <- rownames(subset(metadata_clean, Group == group_names[i]))
      group2_samples <- rownames(subset(metadata_clean, Group == group_names[j]))
      group_distances <- dist_matrix[group1_samples, group2_samples]
      group_distance_matrix[i, j] <- mean(group_distances, na.rm = TRUE)
    }
  }
  
  # Count observations per group
  group_counts <- metadata_clean %>%
    group_by(Group) %>%
    summarise(n = n(), .groups = "drop")
  
  # Filter groups with more than 5 samples
  filtered_groups <- group_counts %>%
    filter(n > 5) %>%
    pull(Group)
  
  # Subset distance matrix to filtered groups
  filtered_distance_matrix <- group_distance_matrix[filtered_groups, filtered_groups]
  
  # Rename row and column names for clarity (Heartwood & Sapwood)
  rownames(filtered_distance_matrix) <- gsub("Inner", "Heartwood", 
                                             gsub("Outer", "Sapwood", 
                                                  rownames(filtered_distance_matrix)))
  colnames(filtered_distance_matrix) <- gsub("Inner", "Heartwood", 
                                             gsub("Outer", "Sapwood", 
                                                  colnames(filtered_distance_matrix)))
  
  # Create annotation for SampleType
  filtered_annotation <- data.frame(
    SampleType = gsub(".*_(Heartwood|Sapwood)", "\\1", rownames(filtered_distance_matrix))
  )
  rownames(filtered_annotation) <- rownames(filtered_distance_matrix)
  
  # Define colors for annotation
  annotation_colors <- list(
    SampleType = c("Heartwood" = "#8B4513", "Sapwood" = "#377EB8") # Distinct colors
  )
  
  # Create the heatmap using pheatmap
  heatmap_plot <- pheatmap(
    filtered_distance_matrix,
    color = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(50)),  # Diverging palette
    clustering_distance_rows = "euclidean",
    fontsize = 10,
    annotation_row = filtered_annotation,
    annotation_colors = annotation_colors
  )
  
  assign(paste0("heatmap_", metric_name), heatmap_plot, envir = .GlobalEnv)
}

# Run heatmap generation for each distance metric
for (metric in names(distance_matrices)) {
  generate_heatmap(distance_matrices[[metric]], metric)
}

# MINIMAL ADDITION 1: After running 16S analysis, save the plots with specific names
barplot_16S_wunifrac <- barplot_wunifrac
heatmap_16S_wunifrac <- heatmap_wunifrac