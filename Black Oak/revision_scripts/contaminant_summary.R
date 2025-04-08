# Load required libraries
library(tidyverse)
library(vegan)  # For rarefaction

# ------------------------------
# Step 1: Data Preprocessing and Contaminant Analysis
# ------------------------------

# Function: Read and analyze OTU table for contaminants before filtering
analyze_contaminants <- function(otu_table_path) {
  # Read the data
  raw_data <- read.table(otu_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get numeric columns (samples)
  sample_cols <- select_if(raw_data, is.numeric) %>% colnames()
  abundance_matrix <- raw_data[, sample_cols]
  
  # Identify chloroplast and mitochondrial sequences
  chloroplast_rows <- which(str_detect(raw_data$Order, "Chloroplast") | 
                              str_detect(raw_data$Family, "Chloroplast") | 
                              str_detect(raw_data$Genus, "Chloroplast"))
  
  mitochondria_rows <- which(str_detect(raw_data$Order, "Mitochondria") | 
                               str_detect(raw_data$Family, "Mitochondria") | 
                               str_detect(raw_data$Genus, "Mitochondria"))
  
  # Extract abundance data for chloroplasts and mitochondria
  chloroplast_abundance <- matrix(0, nrow = length(sample_cols), ncol = 1)
  if(length(chloroplast_rows) > 0) {
    chloroplast_abundance <- colSums(abundance_matrix[chloroplast_rows, , drop = FALSE])
  }
  
  mitochondria_abundance <- matrix(0, nrow = length(sample_cols), ncol = 1)
  if(length(mitochondria_rows) > 0) {
    mitochondria_abundance <- colSums(abundance_matrix[mitochondria_rows, , drop = FALSE])
  }
  
  # Total abundance per sample
  total_abundance <- colSums(abundance_matrix)
  
  # Calculate percentages
  contaminant_data <- data.frame(
    Sample = sample_cols,
    TotalReads = total_abundance,
    ChloroplastReads = chloroplast_abundance,
    MitochondriaReads = mitochondria_abundance,
    ChloroplastPercent = (chloroplast_abundance / total_abundance) * 100,
    MitochondriaPercent = (mitochondria_abundance / total_abundance) * 100
  )
  
  # Add media information
  contaminant_data <- contaminant_data %>%
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
    filter(Media != "Other")
  
  return(contaminant_data)
}

# ------------------------------
# Step 2: Calculate Summary Statistics by Media Type
# ------------------------------

summarize_contaminants_by_media <- function(contaminant_data) {
  contaminant_summary <- contaminant_data %>%
    group_by(Media) %>%
    summarize(
      MeanChloroplastPercent = mean(ChloroplastPercent, na.rm = TRUE),
      SDChloroplastPercent = sd(ChloroplastPercent, na.rm = TRUE),
      MedianChloroplastPercent = median(ChloroplastPercent, na.rm = TRUE),
      MaxChloroplastPercent = max(ChloroplastPercent, na.rm = TRUE),
      
      MeanMitochondriaPercent = mean(MitochondriaPercent, na.rm = TRUE),
      SDMitochondriaPercent = sd(MitochondriaPercent, na.rm = TRUE),
      MedianMitochondriaPercent = median(MitochondriaPercent, na.rm = TRUE),
      MaxMitochondriaPercent = max(MitochondriaPercent, na.rm = TRUE),
      
      SampleCount = n()
    )
  
  return(contaminant_summary)
}

# ------------------------------
# Step A: Analyze contaminants before filtering
# ------------------------------

# Process data for contaminant analysis
contaminant_data <- analyze_contaminants("Black Oak/16S/OTU_table.txt")

# Calculate summary by media type
contaminant_summary <- summarize_contaminants_by_media(contaminant_data)

# ------------------------------
# Step B: Visualization
# ------------------------------

# Create boxplot of chloroplast percentages by media type
chloro_plot <- ggplot(contaminant_data, aes(x = reorder(Media, ChloroplastPercent, FUN = median), y = ChloroplastPercent)) +
  geom_boxplot(fill = "lightgreen", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkgreen") +
  theme_minimal() +
  labs(
    title = "Chloroplast DNA Percentage by Component",
    x = NULL,
    y = "Percentage of Chloroplast DNA (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Create boxplot of mitochondria percentages by media type
mito_plot <- ggplot(contaminant_data, aes(x = reorder(Media, MitochondriaPercent, FUN = median), y = MitochondriaPercent)) +
  geom_boxplot(fill = "lightcoral", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkred") +
  theme_minimal() +
  labs(
    title = "Mitochondrial DNA Percentage by Component",
    x = NULL,
    y = "Percentage of Mitochondrial DNA (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Combined scatter plot
combined_plot <- ggplot(contaminant_data, aes(x = ChloroplastPercent, y = MitochondriaPercent, color = Media)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Chloroplast vs Mitochondrial DNA Percentages",
    x = "Percentage of Chloroplast DNA (%)",
    y = "Percentage of Mitochondrial DNA (%)"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Print summary results
print("Summary of contaminant percentages by media type:")
print(contaminant_summary)

# Save plots
ggsave("chloroplast_boxplot.pdf", plot = chloro_plot, width = 10, height = 6)
ggsave("mitochondria_boxplot.pdf", plot = mito_plot, width = 10, height = 6)
ggsave("contaminant_scatter.pdf", plot = combined_plot, width = 10, height = 8)

# ------------------------------
# Step C: Original analysis pipeline continues below
# ------------------------------

# Function: Read and preprocess OTU table at the Class level
prepare_Genus_level_data <- function(otu_table_path) {
  # To store samples removed
  removed_samples <- character(0)
  # Read the data
  raw_data <- read.table(otu_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Filter out mitochondria and chloroplasts at Order, Family, and Genus levels
  filtered_data <- raw_data %>%
    filter(!str_detect(Order, "Chloroplast|Mitochondria") &
             !str_detect(Family, "Mitochondria|Chloroplast") &
             !str_detect(Genus, "Mitochondria|Chloroplast"))
  
  # Get numeric columns (samples)
  sample_cols <- select_if(filtered_data, is.numeric) %>% colnames()
  abundance_matrix <- filtered_data[, sample_cols]
  
  # Remove samples with less than 3500 reads
  sample_sums <- colSums(abundance_matrix)
  low_depth_samples <- names(sample_sums[sample_sums < 3500])
  
  if(length(low_depth_samples) > 0) {
    message(paste("Removing samples with < 3500 reads:", 
                  paste(low_depth_samples, collapse = ", ")))
    abundance_matrix <- abundance_matrix[, !colnames(abundance_matrix) %in% low_depth_samples]
    removed_samples <- low_depth_samples
  }
  
  # Perform rarefaction
  set.seed(123)  # For reproducibility
  rarefied_matrix <- t(rrarefy(t(abundance_matrix), sample = 3500))
  
  # Aggregate by Genus
  Genus_aggregated <- rarefied_matrix %>%
    data.frame() %>%
    mutate(Genus = filtered_data$Genus) %>%
    group_by(Genus) %>%
    summarize(across(everything(), sum)) %>%
    pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance") %>%
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
    filter(Media != "Other")
  
  return(list(
    data = Genus_aggregated,
    removed_samples = removed_samples
  ))
}