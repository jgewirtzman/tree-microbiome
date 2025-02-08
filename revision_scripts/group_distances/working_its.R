# Load required libraries
# Load or install necessary packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

# Load required packages
p_load(
  phyloseq, ape, dplyr, purrr, tidyr, ggplot2, multcompView, 
  RColorBrewer, pheatmap, emmeans, lme4
)


#### ITS IMPORT #### 

#### Load OTU table ####
otu_tab <- read.delim("ITS/OTU_table.txt", header=TRUE, row.names = 1, stringsAsFactors = FALSE)

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
phylo_tree <- read_tree("ITS/rooted_tree.nwk")

#### Load Sample Metadata ####
samp_data <- read.delim("ITS/tree_its_mapping_dada2_corrected.txt.no_gz", row.names = 1, stringsAsFactors = FALSE)
samp_data$RowName <- row.names(samp_data)

#### Merge Additional Metadata ####
ddpcr <- read.csv(file = "Other-Metadata/ddPCR_meta_all_data.csv", stringsAsFactors = FALSE)
metadata <- read.csv("Other-Metadata/annotated_metadata/ITS_tree_sample_table_with_meta.csv", stringsAsFactors = FALSE)

metadata_filtered <- data.frame(
  SampleID   = metadata$X.1,
  Species    = substr(metadata$Inner.Core.Sample.ID.x, 1, 4),
  SampleType = metadata$core_type
)
metadata_filtered <- subset(metadata_filtered, SampleType %in% c("Inner", "Outer"))
rownames(metadata_filtered) <- metadata_filtered$SampleID


# 1. Rename TreatmentGroup to core_type in samp_data
colnames(samp_data)[colnames(samp_data) == "TreatmentGroup"] <- "core_type"

# 3. Remove duplicate seq_id entries in both datasets
#ddpcr <- ddpcr[!duplicated(ddpcr$seq_id), ]
#samp_data <- samp_data[!duplicated(samp_data$seq_id), ]

# 4. Remove rows with NA values in seq_id
samp_data <- samp_data[!is.na(samp_data$seq_id), ]

# 5. Ensure matching data types for merging
ddpcr$seq_id <- as.character(ddpcr$seq_id)
samp_data$seq_id <- as.character(samp_data$seq_id)

ddpcr$core_type <- as.character(ddpcr$core_type)
samp_data$core_type <- as.character(samp_data$core_type)

# 6. Merge datasets by seq_id and core_type
samp_data_merged <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)

# View merged data
head(samp_data_merged)

samp_data_merged <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)
samp_data_merged <- merge(samp_data_merged, metadata_filtered, by.x = "RowName", by.y = "SampleID", all.x = TRUE)

# Remove duplicates
dups <- which(duplicated(samp_data_merged$RowName))
samp_data_merged <- samp_data_merged[-c(dups), ]

# Set row names
row.names(samp_data_merged) <- samp_data_merged$RowName

#### Construct Phyloseq Object ####
raw_ps <- phyloseq(otu_table_pre, tax_tab_pre, phylo_tree, sample_data(samp_data_merged))

set.seed(46814)
ps.rare <- rarefy_even_depth(raw_ps, sample.size = 4000)

# Filter only Inner & Outer samples
ps.rare <- prune_samples(sample_data(ps.rare)$core_type %in% c("Inner", "Outer"), ps.rare)
ps.rare <- prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)
ps.rare <- prune_taxa(taxa_sums(ps.rare) >= 10, ps.rare)
ps.rare <- prune_taxa(rowSums(otu_table(ps.rare) > 0) >= 2, ps.rare)

phyloseq_clean <- ps.rare  # The final cleaned phyloseq object

# Ensure tree tips match the remaining OTUs
remaining_taxa <- taxa_names(ps.rare)
ps.rare <- prune_taxa(remaining_taxa, ps.rare)
pruned_tree <- prune_taxa(remaining_taxa, phy_tree(ps.rare))
phy_tree(ps.rare) <- pruned_tree
ps.rare <- prune_samples(sample_sums(ps.rare) > 0, ps.rare)

# Define a minimum read count threshold (e.g., 1,000 reads per sample)
min_reads <- 1000
ps.rare <- prune_samples(sample_sums(ps.rare) >= min_reads, ps.rare)

phyloseq_clean <- ps.rare

phyloseq_clean

#### Calculate Distance Matrices ####
wunifrac_dist_clean_ITS <- phyloseq::distance(phyloseq_clean, method = "wunifrac")
range(wunifrac_dist_clean_ITS); mean(wunifrac_dist_clean_ITS)

uunifrac_dist_clean_ITS <- phyloseq::distance(phyloseq_clean, method = "unifrac")
range(uunifrac_dist_clean_ITS); mean(uunifrac_dist_clean_ITS)

bray_dist_clean_ITS <- phyloseq::distance(phyloseq_clean, method = "bray")
range(bray_dist_clean_ITS); mean(bray_dist_clean_ITS)



metadata_clean <- data.frame(sample_data(phyloseq_clean)) %>%
  filter(!is.na(Species) & !is.na(SampleType))

metadata_clean$Group <- paste(metadata_clean$Species, metadata_clean$SampleType, sep = "_")
group_names <- unique(metadata_clean$Group)

unifrac_matrix <- as.matrix(wunifrac_dist_clean_ITS)







# List of distance matrices
distance_matrices <- list(
  wunifrac = as.matrix(wunifrac_dist_clean_ITS),
  uunifrac = as.matrix(uunifrac_dist_clean_ITS),
  braycurtis = as.matrix(bray_dist_clean_ITS)
)

# Function to process each distance matrix and generate heatmap & barplot
analyze_distance_matrix <- function(dist_matrix, metric_name) {
  
  metadata_clean <- data.frame(sample_data(phyloseq_clean)) %>%
    filter(!is.na(Species) & !is.na(SampleType))
  
  metadata_clean$Group <- paste(metadata_clean$Species, metadata_clean$SampleType, sep = "_")
  group_names <- unique(metadata_clean$Group)
  
  # Convert distance matrix to long format
  dist_df <- as.data.frame(as.matrix(dist_matrix))
  dist_df$Sample1 <- rownames(dist_df)
  
  dist_long <- dist_df %>%
    pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance") %>%
    filter(Sample1 != Sample2) %>%
    left_join(metadata_clean, by = c("Sample1" = "RowName")) %>%
    rename(Species1 = Species, SampleType1 = SampleType) %>%
    left_join(metadata_clean, by = c("Sample2" = "RowName")) %>%
    rename(Species2 = Species, SampleType2 = SampleType)
  
  # Define comparison groups
  dist_long <- dist_long %>%
    mutate(
      Comparison = case_when(
        Species1 == Species2 & SampleType1 == "Outer" & SampleType2 == "Inner" ~ "Sapwood_Heartwood_Within",
        Species1 == Species2 & SampleType1 == "Outer" & SampleType2 == "Outer" ~ "Sapwood_Sapwood_Within",
        Species1 == Species2 & SampleType1 == "Inner" & SampleType2 == "Inner" ~ "Heartwood_Heartwood_Within",
        Species1 != Species2 & SampleType1 == "Outer" & SampleType2 == "Inner" ~ "Sapwood_Heartwood_Between",
        Species1 != Species2 & SampleType1 == "Outer" & SampleType2 == "Outer" ~ "Sapwood_Sapwood_Between",
        Species1 != Species2 & SampleType1 == "Inner" & SampleType2 == "Inner" ~ "Heartwood_Heartwood_Between",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Comparison))
  
  # Aggregate distances
  group_means <- dist_long %>%
    group_by(Comparison) %>%
    summarise(
      MeanDistance = mean(Distance, na.rm = TRUE),
      SD = sd(Distance, na.rm = TRUE),
      UniquePairs = n_distinct(pmap_chr(list(Sample1, Sample2), ~paste(sort(c(...)), collapse = "_"))),
      .groups = "drop"
    ) %>%
    mutate(SE = SD / sqrt(UniquePairs))
  
  # Perform ANOVA and Tukey's HSD only if there are enough comparisons
  if (nrow(group_means) > 1) {
    anova_model <- aov(Distance ~ Comparison, data = dist_long)
    tukey_result <- TukeyHSD(anova_model)
    tukey_letters <- multcompLetters4(anova_model, tukey_result, reversed=TRUE)
    
    if (!is.null(tukey_letters$Comparison)) {
      group_letters <- data.frame(
        Comparison = names(tukey_letters$Comparison$Letters),
        Letters = tukey_letters$Comparison$Letters
      )
      group_means <- group_means %>%
        left_join(group_letters, by = "Comparison")
    } else {
      group_means$Letters <- ""  # If no letters were generated, assign an empty string
    }
  } else {
    group_means$Letters <- ""  # If ANOVA cannot be performed, assign empty letters
  }
  
  # Ensure the correct order for comparisons
  comparison_order <- c(
    "Sapwood_Heartwood_Within",
    "Sapwood_Sapwood_Within",
    "Heartwood_Heartwood_Within",
    "Sapwood_Heartwood_Between",
    "Sapwood_Sapwood_Between",
    "Heartwood_Heartwood_Between"
  )
  group_means$Comparison <- factor(group_means$Comparison, levels = comparison_order)
  
  # Create a mapping for new x-axis labels (preserve 6 groups, just relabel)
  comparison_labels <- c(
    "Sapwood_Heartwood_Within" = "Sapwood-Heartwood",
    "Sapwood_Sapwood_Within" = "Sapwood-Sapwood",
    "Heartwood_Heartwood_Within" = "Heartwood-Heartwood",
    "Sapwood_Heartwood_Between" = "Sapwood-Heartwood",
    "Sapwood_Sapwood_Between" = "Sapwood-Sapwood",
    "Heartwood_Heartwood_Between" = "Heartwood-Heartwood"
  )
  
  # Define the correct factor levels to preserve plotting order
  comparison_order <- c(
    "Sapwood_Heartwood_Within",
    "Sapwood_Sapwood_Within",
    "Heartwood_Heartwood_Within",
    "Sapwood_Heartwood_Between",
    "Sapwood_Sapwood_Between",
    "Heartwood_Heartwood_Between"
  )
  
  group_means$Comparison <- factor(group_means$Comparison, levels = comparison_order)
  
  # Define bracket positions
  y_max <- max(group_means$MeanDistance + group_means$SE, na.rm = TRUE)  # Highest value for placement
  bracket_height <- 0.05  # Adjust if needed
  
  # Define y-axis labels for each metric
  y_axis_labels <- c(
    "wunifrac" = "WUnifrac Distance",
    "uunifrac" = "Unweighted Unifrac Distance",
    "braycurtis" = "Bray-Curtis Distance"
  )
  
  # Generate barplot
  barplot_plot <- ggplot(group_means, aes(x = Comparison, y = MeanDistance, fill = Comparison)) +
    geom_bar(stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = MeanDistance - SE, ymax = MeanDistance + SE), width = 0.2) +
    geom_text(aes(y = MeanDistance + SE + 0.02, label = Letters), size = 5, hjust = 0.5) +
    
    # Bracket for "Within Species"
    geom_segment(aes(x = 1, xend = 3, y = y_max + bracket_height, yend = y_max + bracket_height), size = 0.7) +
    annotate("text", x = 2, y = y_max + bracket_height + 0.02, label = "Within Species", size = 5) +
    
    # Bracket for "Between Species"
    geom_segment(aes(x = 4, xend = 6, y = y_max + bracket_height, yend = y_max + bracket_height), size = 0.7) +
    annotate("text", x = 5, y = y_max + bracket_height + 0.02, label = "Between Species", size = 5) +
    
    scale_x_discrete(labels = comparison_labels) +  # Apply new labels while keeping 6 bars
    scale_fill_manual(values = c(
      "Sapwood_Heartwood_Within" = "grey50",
      "Sapwood_Sapwood_Within" = "#377EB8",
      "Heartwood_Heartwood_Within" = "#8B4513",
      "Sapwood_Heartwood_Between" = "grey50",
      "Sapwood_Sapwood_Between" = "#377EB8",
      "Heartwood_Heartwood_Between" = "#8B4513"
    )) +
    labs(x = "Comparison", y = y_axis_labels[[metric_name]]) +  # Dynamically set y-axis label
    theme_classic() +
    theme(legend.position = "none")
  
  assign(paste0("barplot_", metric_name), barplot_plot, envir = .GlobalEnv)
  
  # Explicitly print the barplot
  print(barplot_plot)
  
  
}

# Run for each distance metric
for (metric in names(distance_matrices)) {
  analyze_distance_matrix(distance_matrices[[metric]], metric)
}





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
