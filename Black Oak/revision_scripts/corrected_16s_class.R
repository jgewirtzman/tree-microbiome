# Load required libraries
library(tidyverse)
library(vegan)  # For rarefaction

# ------------------------------
# Step 1: Data Preprocessing
# ------------------------------

# Function: Read and preprocess OTU table at the Class level
prepare_class_level_data <- function(otu_table_path) {
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
  
  # Aggregate by Class
  class_aggregated <- rarefied_matrix %>%
    data.frame() %>%
    mutate(Class = filtered_data$Class) %>%
    group_by(Class) %>%
    summarize(across(everything(), sum)) %>%
    pivot_longer(-Class, names_to = "Sample", values_to = "Abundance") %>%
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
    data = class_aggregated,
    removed_samples = removed_samples
  ))
}

# ------------------------------
# Step 2: Calculate Relative Abundance
# ------------------------------

prepare_relative_abundance <- function(class_data) {
  class_data %>%
    group_by(Media, Sample) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance) * 100) %>%
    group_by(Class, Media) %>%
    summarize(
      RelativeAbundance = mean(RelativeAbundance),
      .groups = "drop"
    )
}

# ------------------------------
# Step 3: Process Data
# ------------------------------

# Read and preprocess data
processed_data <- prepare_class_level_data("Black Oak/16S/OTU_table.txt")
class_level_data <- prepare_relative_abundance(processed_data$data)
low_depth_samples <- processed_data$removed_samples

# Get top 25 classes by total relative abundance across all components
top_25_classes <- class_level_data %>%
  group_by(Class) %>%
  summarize(
    TotalAbundance = sum(RelativeAbundance)
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25) %>%
  pull(Class)

# Format data for plotting
plot_data <- class_level_data %>%
  filter(Class %in% top_25_classes) %>%
  mutate(
    Class = factor(Class, levels = rev(top_25_classes)),
    Media = factor(Media, levels = c("Heartwood", "Sapwood", "Bark", "Branch", 
                                     "Foliage", "Leaf Litter", "Coarse Root", 
                                     "Fine Root", "Rot", "Organic Soil", 
                                     "Mineral Soil"))
  ) %>%
  complete(Class, Media, fill = list(RelativeAbundance = 0))

# ------------------------------
# Step 4: Create Plot
# ------------------------------

p <- ggplot(plot_data, aes(x = Media, y = Class, fill = RelativeAbundance)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"),
    name = "Relative\nAbundance (%)",
    trans = "log1p",
    breaks = c(0, 0.5, 1, 2, 5, 10, 20, 50),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 12, margin = margin(t = 20)),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm")
  ) +
  labs(x = NULL) +
  scale_y_discrete(position = "left")

p

# ------------------------------
# Step 5: Save Plot
# ------------------------------

# Save with specified dimensions
ggsave("heatmap.pdf", 
       plot = p,
       width = 12,
       height = 10,
       limitsize = FALSE,
       device = cairo_pdf)

# Print diagnostic information
print("Samples removed due to low depth:")
print(low_depth_samples)

print("\nTop 25 classes by total abundance:")
print(top_25_classes)

# Print summary of relative abundances
print("\nSummary of relative abundances by media type:")
plot_data %>%
  group_by(Media) %>%
  summarise(
    Total = sum(RelativeAbundance),
    Mean = mean(RelativeAbundance),
    Max = max(RelativeAbundance)
  ) %>%
  print(n = Inf)