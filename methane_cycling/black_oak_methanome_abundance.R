# Load libraries
library(tidyverse)
library(viridis)
library(cluster) # For clustering columns (media)

# Step 1: Load and Prepare Data
otu_table <- read.table("/Users/jongewirtzman/Downloads/OTU_table (2).txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define methanogen families
methanogen_families <- c(
  "Bathyarchaeia Class", 
  "Methanobacteriaceae", 
  "Methanocellaceae", 
  "Methanocorpusculaceae", 
  "Methanomassiliicoccaceae", 
  "Methanomethylophilaceae", 
  "Methanomicrobiaceae", 
  "Methanoregulaceae", 
  "Methanosaetaceae", 
  "Methanosarcinaceae", 
  "Rice Cluster II", 
  "Methanopyraceae", 
  "Methanothermaceae", 
  "Methanocaldococcaceae", 
  "Methanoplanaceae", 
  "Methanospirillaceae"
)

# Define methanotroph families
methanotroph_families <- c(
  "Methylococcaceae", 
  "Methylacidiphilaceae", 
  "Beijerinckiaceae", 
  "Methylocystaceae", 
  "Methylomonadaceae", 
  "Methylophilaceae", 
  "Methylopilaceae", 
  "Methyloligellaceae", 
  "Methylomirabilaceae"
)

# Add explicit methanotroph genera
methanotroph_genera <- c(
  "Methylobacterium-Methylorubrum", 
  "Methylocella", 
  "Methylorosula", 
  "Methylocapsa", 
  "1174-901-12", # Uncertain
  "Roseiarcus"   # Uncertain
)

# Filter for relevant families and genera
otu_table <- otu_table %>%
  filter(!is.na(Family)) %>% 
  mutate(Group = case_when(
    Family %in% methanogen_families ~ "Methanogen",
    (Family %in% methanotroph_families & !(Family %in% c("Methylopilaceae", "Beijerinckiaceae") & Genus %in% c("Bosea", "Microvirga", "Psychroglaciecola"))) |
      Genus %in% methanotroph_genera ~ "Methanotroph",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

# Pivot and aggregate data
data_cols <- otu_table %>%
  select(where(is.numeric)) %>%
  colnames()

otu_table_long <- otu_table %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Sample", values_to = "Abundance") %>%
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

aggregated_data <- otu_table_long %>%
  group_by(Group, Family, Genus, Media) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Media, values_from = Total_Abundance, values_fill = 0)

# Log transform the data
log_transformed_data <- aggregated_data %>%
  mutate(across(where(is.numeric), ~ log10(. + 1)))

# Prepare data for ggplot (without Group in row labels)
heatmap_data <- log_transformed_data %>%
  pivot_longer(cols = -c(Group, Family, Genus), names_to = "Media", values_to = "Abundance") %>%
  mutate(Row_Label = paste(Family, Genus, sep = " | ")) # Exclude Group

# Cluster columns (media types)
media_dist <- dist(t(as.matrix(log_transformed_data[, -(1:3)])))
media_clust <- hclust(media_dist)
media_order <- media_clust$labels[media_clust$order]

heatmap_data <- heatmap_data %>%
  mutate(Media = factor(Media, levels = media_order))

# Plot the heatmap
ggplot(heatmap_data, aes(x = Media, y = Row_Label, fill = Abundance)) +
  geom_tile() +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") + # Split by Methanogen/Methanotroph
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    strip.text.y = element_text(size = 10, angle = 0)
  ) +
  labs(
    x = "",
    y = "Taxa (Family | Genus)",
    title = "Methanogens and Methanotrophs Abundance (Log Scale)"
  ) +
  scale_fill_viridis_c(option = "C") # Unified scale if necessary





ggplot(heatmap_data, aes(x = Media, y = Row_Label, fill = Abundance)) +
  geom_tile(color = "white") + # Add white borders between cells
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") + # Split by Methanogen/Methanotroph
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    strip.text.y = element_text(size = 10, angle = 0),
    axis.title.x = element_blank() # Remove x-axis label
  ) +
  labs(
    x = "",
    y = "Taxa (Family | Genus)"
  ) +
  scale_fill_viridis_c(option = "C") # Unified color scale

