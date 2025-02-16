# Calculate relative abundance
otu_table_long_relative <- otu_table_long %>%
  group_by(Sample) %>% 
  mutate(Relative_Abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>% # Normalize within each sample
  ungroup()

# Aggregate data for visualization
aggregated_relative_data <- otu_table_long_relative %>%
  group_by(Group, Family, Genus, Media) %>%
  summarize(Total_Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Media, values_from = Total_Relative_Abundance, values_fill = 0)

# Log transform the relative abundance data (optional)
log_transformed_relative_data <- aggregated_relative_data %>%
  mutate(across(where(is.numeric), ~ log10(. + 1)))

# Prepare data for ggplot (without Group in row labels)
heatmap_relative_data <- log_transformed_relative_data %>%
  pivot_longer(cols = -c(Group, Family, Genus), names_to = "Media", values_to = "Abundance") %>%
  mutate(Row_Label = paste(Family, Genus, sep = " | ")) # Exclude Group

# Cluster columns (media types)
media_dist_relative <- dist(t(as.matrix(log_transformed_relative_data[, -(1:3)])))
media_clust_relative <- hclust(media_dist_relative)
media_order_relative <- media_clust_relative$labels[media_clust_relative$order]

heatmap_relative_data <- heatmap_relative_data %>%
  mutate(Media = factor(Media, levels = media_order_relative))

# Plot heatmap for relative abundance
ggplot(heatmap_relative_data, aes(x = Media, y = Row_Label, fill = Abundance)) +
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
