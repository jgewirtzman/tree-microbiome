# Load required libraries
library(tidyverse)

# ------------------------------
# Step 1: Data Preprocessing
# ------------------------------

# Function: Read and preprocess OTU table at the genus level
prepare_genus_level_data <- function(otu_table_path) {
  read.table(otu_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    pivot_longer(cols = where(is.numeric), names_to = "Sample", values_to = "Abundance") %>%
    mutate(Media = case_when(
      str_detect(Sample, "HEART") ~ "Heartwood",
      str_detect(Sample, "SAP") ~ "Sapwood",
      str_detect(Sample, "LITTER") ~ "Leaf Litter",
      str_detect(Sample, "BARK") ~ "Bark",
      str_detect(Sample, "FOLIAGE") ~ "Foliage",
      str_detect(Sample, "MINERAL") ~ "Mineral Soil",
      str_detect(Sample, "ORGANIC") ~ "Organic Soil",
      str_detect(Sample, "BRANCH") ~ "Branch",
      str_detect(Sample, "COARSE") ~ "Coarse Root",
      str_detect(Sample, "FINE") ~ "Fine Root",
      str_detect(Sample, "ROT") ~ "Rot",
      TRUE ~ "Other"
    )) %>%
    filter(Media != "Other") %>% # Remove "Other" samples
    filter(Abundance > 0) %>% # Remove zero-abundance rows
    filter(!Genus %in% c("Mitochondria", "Chloroplast")) %>% # Filter mitochondria and chloroplasts
    select(Genus, Media, Abundance) %>% # Keep relevant columns
    group_by(Genus, Media) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") # Summarize by genus and media
}

# Step 2: Prepare relative abundance data
prepare_relative_abundance <- function(genus_level_data) {
  genus_level_data %>%
    group_by(Media) %>% # Group by media to calculate relative abundance
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>% # Normalize within each medium
    ungroup()
}

# Step 3: Retain only top 5 genera per medium
filter_top_genera <- function(relative_abundance_data) {
  relative_abundance_data %>%
    group_by(Media) %>% # Group by medium
    arrange(desc(RelativeAbundance)) %>% # Sort genera by relative abundance within each medium
    mutate(Rank = row_number()) %>% # Assign rank to each genus within the medium
    filter(Rank <= 5) %>% # Keep only the top 5 genera per medium
    ungroup() %>%
    select(-Rank) # Remove the rank column
}

# Step 4: Rank genera by mean relative abundance
rank_genera <- function(relative_abundance_data) {
  relative_abundance_data %>%
    group_by(Genus) %>%
    summarize(MeanRelativeAbundance = mean(RelativeAbundance, na.rm = TRUE)) %>% # Calculate mean abundance
    arrange(desc(MeanRelativeAbundance)) %>% # Sort descending by mean abundance
    pull(Genus) # Extract the ranked genus names
}

# ------------------------------
# Step 5: Load Data and Preprocess
# ------------------------------

# Read and preprocess data
genus_level_data <- prepare_genus_level_data("Black Oak/16S/OTU_table.txt")
genus_level_data <- prepare_relative_abundance(genus_level_data)
genus_level_data <- filter_top_genera(genus_level_data)

# Rank genera by mean relative abundance
ranked_genera <- rank_genera(genus_level_data)

# Determine max relative abundance for scaling
max_relative_abundance <- max(genus_level_data$RelativeAbundance)

# Pivot data to wide format
heatmap_data <- genus_level_data %>%
  pivot_wider(names_from = Media, values_from = RelativeAbundance, values_fill = 0)

# Perform hierarchical clustering on media (columns)
media_dist <- dist(t(as.matrix(heatmap_data[-1])))  # Distance matrix on transposed data
media_cluster <- hclust(media_dist)                # Hierarchical clustering

# Reorder media based on clustering
ordered_media <- media_cluster$labels[media_cluster$order]

# ------------------------------
# Step 6: Final Preprocessing
# ------------------------------

# Reformat data for ggplot
heatmap_data_long <- genus_level_data %>%
  mutate(Media = factor(Media, levels = ordered_media), # Order media
         Genus = factor(Genus, levels = rev(ranked_genera))) %>% # Reverse order of genera
  complete(Genus, Media, fill = list(RelativeAbundance = 0)) %>% # Ensure all combinations exist
  filter(!is.na(Genus) & Genus != "Unknown Genus" & Genus != "") %>% # Remove unknown/blank genera
  filter(Media != "Abundance") # Remove the "Abundance" column

# ------------------------------
# Step 7: Create the Heatmap
# ------------------------------

ggplot(heatmap_data_long, aes(x = Media, y = Genus, fill = RelativeAbundance)) +
  geom_tile(color = "white") + # Add borders between tiles for better separation
  scale_fill_gradientn(
    colors = c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"),
    values = scales::rescale(c(0.01, 0.20, max_relative_abundance)), # Rescale the range for 1%-20%
    limits = c(0.01, 0.10), # Limit scale to 1% to 20%
    oob = scales::squish, # Squish values above 20% to the maximum shade
    na.value = "white",
    name = "Relative Abundance (%)" # Add legend title
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Rotate and enlarge x-axis labels
    axis.text.y = element_text(size = 10), # Enlarge y-axis labels
    axis.title = element_text(size = 14), # Enlarge axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Center and enlarge title
    legend.title = element_text(size = 12), # Enlarge legend title
    legend.text = element_text(size = 10) # Enlarge legend text
  )
