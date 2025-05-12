library(ggplot2)
library(patchwork)
library(dplyr)
library(forcats)
library(tidyr)

# Filter the data to include only relevant columns and rows
df_filtered <- GC_data %>%
  dplyr::select(Species, Tree.No, Tissue, CH4_concentration) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(CH4_concentration)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(CH4_concentration = as.numeric(CH4_concentration)) %>%  # Convert CH4_concentration to numeric
  filter(!is.na(CH4_concentration)) %>%  # Remove rows where CH4_concentration couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(CH4_concentration = mean(CH4_concentration), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = CH4_concentration) %>%
  filter(!is.na(heartwood) & !is.na(sapwood))  # Keep only rows where both heartwood and sapwood values are present

# Ensure proper ordering by heartwood values (descending)
df_wide <- df_wide %>%
  arrange(Species, desc(heartwood)) %>%  # Sort by species and descending heartwood values
  group_by(Species) %>%
  mutate(Tree.No = as.character(row_number())) %>%  # Assign ordered tree numbers
  ungroup() %>%
  mutate(Tree.No = fct_reorder(Tree.No, heartwood, .desc = TRUE))  # Reorder factor by descending heartwood

# Set Species as a factor with specific order
df_wide$Species <- factor(df_wide$Species, 
                          levels = c("Acer rubrum", "Acer saccharum", "Betula lenta", 
                                     "Pinus strobus", "Quercus rubra", "Tsuga canadensis"))

# Function to create base plot with an inset for log-transformed values
create_facet_with_inset <- function(species_data) {
  species_name <- unique(species_data$Species)
  
  base_plot <- ggplot(species_data, aes(y = Tree.No)) +
    geom_segment(aes(x = heartwood, xend = sapwood, yend = Tree.No), color = "grey30") +
    geom_point(aes(x = heartwood, color = "Heartwood"), size = 2.5, alpha = 0.75) +
    geom_point(aes(x = sapwood, color = "Sapwood"), size = 2.5, alpha = 0.75) +
    scale_color_manual(
      values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),
      name = "Tissue"
    ) +
    labs(x = expression(CH[4] * " (ppm)"), y = "Trees (Ranked)") +
    labs(title = species_name) +
    theme_classic() +
    theme(
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(size = 10, hjust = 0),
      axis.title.x = element_text(size = 8)
    ) + coord_flip()
  
  inset_plot <- ggplot(species_data, aes(y = Tree.No)) +
    geom_segment(aes(x = log10(pmax(heartwood, 1, na.rm = TRUE)), xend = log10(pmax(sapwood, 1, na.rm = TRUE)), yend = Tree.No), color = "grey30") +
    geom_point(aes(x = log10(pmax(sapwood, 1, na.rm = TRUE)), color = "Sapwood"), size = 1.5, alpha = 0.75) +
    geom_point(aes(x = log10(pmax(heartwood, 1, na.rm = TRUE)), color = "Heartwood"), size = 1.5, alpha = 0.75) +
    scale_color_manual(
      values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(0, ceiling(max(log10(pmax(species_data$heartwood, 1, na.rm = TRUE)), 
                                log10(pmax(species_data$sapwood, 1, na.rm = TRUE))))),
      breaks = seq(0, ceiling(max(log10(pmax(species_data$heartwood, 1, na.rm = TRUE)), 
                                  log10(pmax(species_data$sapwood, 1, na.rm = TRUE)))), 1)
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 8),
      plot.margin = margin(2, 2, 2, 2),
      plot.background = element_rect(fill = "white", color = "grey30", size = 1),
      panel.background = element_rect(fill = "white", color = NA),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) + coord_flip()
  
  base_plot + inset_element(inset_plot, left = 0.5, bottom = 0.25, right = 1, top = 1)
}

# Split data by species and apply function to create individual facet plots
species_plots <- lapply(split(df_wide, df_wide$Species), create_facet_with_inset)

# Arrange plots in alphabetical order using wrap_plots with shared legend
final_plot <- wrap_plots(species_plots, ncol = 2) + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -10, r = -10, b = -10, l = -10)
  )

# Display the final combined plot
print(final_plot)
ggsave("dumbell_plot_with_species_labels.png", final_plot, width = 5, height = 5)

ggsave("dumbell_plot_with_species_labels.svg", final_plot, width = 5, height = 5)


# Use plot_grid to add the label "i"
dumbbell_plot <- plot_grid(final_plot, labels = "i", label_size = 14, label_x = 0.00, label_y = 0.95, hjust = 0, vjust = 1)

# Display the plot
dumbbell_plot