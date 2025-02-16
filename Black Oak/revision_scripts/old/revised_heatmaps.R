# Load required libraries
library(tidyverse)

# ------------------------------
# Step 1: Data Preprocessing
# ------------------------------

prepare_genus_level_data <- function(otu_table_path) {
  read.table(
    file = otu_table_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  ) %>%
    # Pivot ONLY the sample columns.
    # Adjust the starts_with("QUVE") pattern to match your actual sample column prefixes.
    pivot_longer(
      cols = starts_with("QUVE"),    # Modify this pattern based on your actual data
      names_to = "Sample",
      values_to = "Abundance"
    ) %>%
    # Normalize Sample names by removing suffixes like '.S512'
    mutate(
      BaseSample = str_remove(Sample, "\\.S\\d+$"),  # Removes '.S' followed by digits at the end
      Media = case_when(
        str_detect(BaseSample, regex("HEART", ignore_case = TRUE))   ~ "Heartwood",
        str_detect(BaseSample, regex("SAP", ignore_case = TRUE))     ~ "Sapwood",
        str_detect(BaseSample, regex("LITTER", ignore_case = TRUE))  ~ "Leaf Litter",
        str_detect(BaseSample, regex("BARK", ignore_case = TRUE))    ~ "Bark",
        str_detect(BaseSample, regex("FOLIAGE", ignore_case = TRUE)) ~ "Foliage",
        str_detect(BaseSample, regex("MINERAL", ignore_case = TRUE)) ~ "Mineral Soil",
        str_detect(BaseSample, regex("ORGANIC", ignore_case = TRUE)) ~ "Organic Soil",
        str_detect(BaseSample, regex("BRANCH", ignore_case = TRUE))  ~ "Branch",
        str_detect(BaseSample, regex("COARSE", ignore_case = TRUE))  ~ "Coarse Root",
        str_detect(BaseSample, regex("FINE", ignore_case = TRUE))    ~ "Fine Root",
        str_detect(BaseSample, regex("ROT", ignore_case = TRUE))     ~ "Rot",
        TRUE ~ NA_character_
      )
    ) %>%
    # Filter out rows with NA Media
    filter(!is.na(Media)) %>%
    # Replace Abundance < 1 with 0 instead of filtering them out
    mutate(Abundance = ifelse(Abundance > 0, Abundance, 0)) %>%
    # Exclude unwanted genera
    filter(!Class %in% c("Mitochondria", "Chloroplast")) %>%
    # Remove rows with unknown or blank genera
    filter(!is.na(Genus) & Genus != "") %>%
    # Aggregate Abundance for the same Genus and BaseSample
    group_by(Genus, Media, BaseSample) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    # Rename BaseSample back to Sample for consistency
    rename(Sample = BaseSample) %>%
    select(Genus, Media, Sample, Abundance)
}


# Step 2: Prepare relative abundance data
prepare_relative_abundance <- function(genus_level_data) {
  genus_level_data %>%
    group_by(Sample) %>% # Group by individual samples to calculate relative abundance
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>% # Normalize within each sample
    ungroup()
}

# ------------------------------
# Step 3: Filter Top Genera
# ------------------------------

filter_top_genera <- function(genus_level_data) {
  genus_level_data %>%
    group_by(Media, Genus) %>%
    summarize(TotalAbundance = sum(RelativeAbundance), .groups = "drop") %>%
    group_by(Media) %>%
    arrange(desc(TotalAbundance)) %>%
    mutate(Rank = row_number()) %>%
    filter(Rank <= 3) %>%
    pull(Genus)
}

# ------------------------------
# Step 4: Load Data and Preprocess
# ------------------------------

# Read and preprocess data
genus_level_data <- prepare_genus_level_data("Black Oak/16S/OTU_table.txt")
genus_level_data <- prepare_relative_abundance(genus_level_data)

# Get the top genera across media (ensure uniqueness)
top_genera <- unique(filter_top_genera(genus_level_data))

# Keep only the top genera for plotting
heatmap_data <- genus_level_data %>%
  filter(Genus %in% top_genera) %>%
  mutate(
    Sample = factor(Sample, levels = unique(Sample)),      # Retain individual sample order
    Genus  = factor(Genus, levels = rev(top_genera)),     # Reverse order of unique top genera
    Media  = factor(Media, levels = unique(Media))
  ) %>%
  # Optional: If you still want to include all combinations, ensure Media is included in complete()
  # Otherwise, remove complete() to prevent introducing duplicated or unwanted samples
  # complete(Genus, Sample, fill = list(RelativeAbundance = 0)) 
  distinct(Genus, Sample, .keep_all = TRUE)  # Ensure no duplicate Genus-Sample pairs

# ------------------------------
# Step 5: Create the Heatmap
# ------------------------------

# Modify the Media labels to include line breaks
heatmap_data <- heatmap_data %>%
  mutate(Media = str_replace_all(Media, " ", "\n"))

# Generate the heatmap
ggplot(heatmap_data, aes(x = Sample, y = Genus, fill = RelativeAbundance)) +
  geom_tile(color = "white") + # Add borders between tiles
  scale_fill_gradientn(
    colors = c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"),
    # Dynamic scaling based on relative abundance
    values = scales::rescale(c(0.0, 0.1, max(heatmap_data$RelativeAbundance, na.rm = TRUE))),
    limits = c(0.01, 0.25), # Adjust color scale range as needed
    oob = scales::squish,    # Handle out-of-bound values
    na.value = "#4575B4",
    name = "Relative Abundance (%)"
  ) +
  facet_grid(. ~ Media, scales = "free_x", space = "free_x") + # Group samples by Media
  theme_minimal() +
  theme(
    axis.text.x    = element_text(angle = 90, hjust = 1, size = 8), # Rotate x-axis labels for samples
    axis.text.y    = element_text(size = 10),
    axis.title     = element_text(size = 14),
    plot.title     = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title   = element_text(size = 12),
    legend.text    = element_text(size = 10),
    strip.text.x   = element_text(size = 10, face = "bold") # Style facet labels
  ) +
  labs(
    x = "Sample",
    y = "Genus"
  )








###







# Load required libraries
library(tidyverse)

# ------------------------------
# Step 1: Data Preprocessing
# ------------------------------

prepare_Class_level_data <- function(otu_table_path) {
  read.table(
    file = otu_table_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  ) %>%
    # Pivot ONLY the sample columns.
    # Adjust the starts_with("QUVE") pattern to match your actual sample column prefixes.
    pivot_longer(
      cols = starts_with("QUVE"),    # Modify this pattern based on your actual data
      names_to = "Sample",
      values_to = "Abundance"
    ) %>%
    # Normalize Sample names by removing suffixes like '.S512'
    mutate(
      BaseSample = str_remove(Sample, "\\.S\\d+$"),  # Removes '.S' followed by digits at the end
      Media = case_when(
        str_detect(BaseSample, regex("HEART", ignore_case = TRUE))   ~ "Heartwood",
        str_detect(BaseSample, regex("SAP", ignore_case = TRUE))     ~ "Sapwood",
        str_detect(BaseSample, regex("LITTER", ignore_case = TRUE))  ~ "Leaf Litter",
        str_detect(BaseSample, regex("BARK", ignore_case = TRUE))    ~ "Bark",
        str_detect(BaseSample, regex("FOLIAGE", ignore_case = TRUE)) ~ "Foliage",
        str_detect(BaseSample, regex("MINERAL", ignore_case = TRUE)) ~ "Mineral Soil",
        str_detect(BaseSample, regex("ORGANIC", ignore_case = TRUE)) ~ "Organic Soil",
        str_detect(BaseSample, regex("BRANCH", ignore_case = TRUE))  ~ "Branch",
        str_detect(BaseSample, regex("COARSE", ignore_case = TRUE))  ~ "Coarse Root",
        str_detect(BaseSample, regex("FINE", ignore_case = TRUE))    ~ "Fine Root",
        str_detect(BaseSample, regex("ROT", ignore_case = TRUE))     ~ "Rot",
        TRUE ~ NA_character_
      )
    ) %>%
    # Filter out rows with NA Media
    filter(!is.na(Media)) %>%
    # Replace Abundance < 1 with 0 instead of filtering them out
    mutate(Abundance = ifelse(Abundance > 0, Abundance, 0)) %>%
    # Exclude unwanted genera
    filter(!Class %in% c("Mitochondria", "Chloroplast")) %>%
    # Remove rows with unknown or blank genera
    filter(!is.na(Class) & Class != "") %>%
    # Aggregate Abundance for the same Class and BaseSample
    group_by(Class, Media, BaseSample) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    # Rename BaseSample back to Sample for consistency
    rename(Sample = BaseSample) %>%
    select(Class, Media, Sample, Abundance)
}


# Step 2: Prepare relative abundance data
prepare_relative_abundance <- function(Class_level_data) {
  Class_level_data %>%
    group_by(Sample) %>% # Group by individual samples to calculate relative abundance
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>% # Normalize within each sample
    ungroup()
}

# ------------------------------
# Step 3: Filter Top Genera
# ------------------------------

filter_top_genera <- function(Class_level_data) {
  Class_level_data %>%
    group_by(Media, Class) %>%
    summarize(TotalAbundance = sum(RelativeAbundance), .groups = "drop") %>%
    group_by(Media) %>%
    arrange(desc(TotalAbundance)) %>%
    mutate(Rank = row_number()) %>%
    filter(Rank <= 5) %>%
    pull(Class)
}

# ------------------------------
# Step 4: Load Data and Preprocess
# ------------------------------

# Read and preprocess data
Class_level_data <- prepare_Class_level_data("Black Oak/16S/OTU_table.txt")
Class_level_data <- prepare_relative_abundance(Class_level_data)

# Get the top genera across media (ensure uniqueness)
top_genera <- unique(filter_top_genera(Class_level_data))

# Keep only the top genera for plotting
heatmap_data <- Class_level_data %>%
  filter(Class %in% top_genera) %>%
  mutate(
    Sample = factor(Sample, levels = unique(Sample)),      # Retain individual sample order
    Class  = factor(Class, levels = rev(top_genera)),     # Reverse order of unique top genera
    Media  = factor(Media, levels = unique(Media))
  ) %>%
  # Optional: If you still want to include all combinations, ensure Media is included in complete()
  # Otherwise, remove complete() to prevent introducing duplicated or unwanted samples
  # complete(Class, Sample, fill = list(RelativeAbundance = 0)) 
  distinct(Class, Sample, .keep_all = TRUE)  # Ensure no duplicate Class-Sample pairs

# ------------------------------
# Step 5: Create the Heatmap
# ------------------------------

# Modify the Media labels to include line breaks
heatmap_data <- heatmap_data %>%
  mutate(Media = str_replace_all(Media, " ", "\n"))

# Generate the heatmap
ggplot(heatmap_data, aes(x = Sample, y = Class, fill = RelativeAbundance)) +
  geom_tile(color = "white") + # Add borders between tiles
  scale_fill_gradientn(
    colors = c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"),
    # Dynamic scaling based on relative abundance
    values = scales::rescale(c(0.0, 0.15, max(heatmap_data$RelativeAbundance, na.rm = TRUE))),
    limits = c(0.01, 0.5), # Adjust color scale range as needed
    oob = scales::squish,    # Handle out-of-bound values
    na.value = "#4575B4",
    name = "Relative Abundance (%)"
  ) +
  facet_grid(. ~ Media, scales = "free_x", space = "free_x") + # Group samples by Media
  theme_minimal() +
  theme(
    axis.text.x    = element_text(angle = 90, hjust = 1, size = 8), # Rotate x-axis labels for samples
    axis.text.y    = element_text(size = 10),
    axis.title     = element_text(size = 14),
    plot.title     = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title   = element_text(size = 12),
    legend.text    = element_text(size = 10),
    strip.text.x   = element_text(size = 10, face = "bold") # Style facet labels
  ) +
  labs(
    x = "Sample",
    y = "Class"
  )

