# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)


library(ggplot2)
library(patchwork)
library(dplyr)
library(forcats)

# Read the original data from CSV
#data <- read.csv("Downloads/GC_241010 - Run 1 Compiled (1).csv", stringsAsFactors = FALSE)
data <- GC_data


###CH4###

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
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

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(heartwood)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = heartwood, xend = sapwood, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = heartwood, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = sapwood, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
  scale_color_manual(
    values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),  # Custom colors for tissues
    name = "Tissue"  # Legend title
  ) +
  facet_wrap(~ Species, scales = "free") +  # Free axes for each species
  labs(
    #title = "CH4 Concentrations by Tree and Tissue Type",
    x = "CH4 (ppm)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),  # Customize facet label size
    axis.text = element_text(size = 10),   # Customize axis text size
    axis.title = element_text(size = 12)   # Customize axis title size
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "top")




###CH4 Log###

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
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

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(heartwood)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

df_wide$heartwood<-log10(df_wide$heartwood)
df_wide$sapwood<-log10(df_wide$sapwood)

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = heartwood, xend = sapwood, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = heartwood, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = sapwood, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
  scale_color_manual(
    values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),  # Custom colors for tissues
    name = "Tissue"  # Legend title
  ) +
  facet_wrap(~ Species, scales = "free") +  # Free axes for each species
  labs(
    #title = "CH4 Concentrations by Tree and Tissue Type",
    x = "Log10 CH4 (ppm)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),  # Customize facet label size
    axis.text = element_text(size = 10),   # Customize axis text size
    axis.title = element_text(size = 12)   # Customize axis title size
  )+  theme(axis.text.x = element_text(angle = 45, hjust = 1))





###O2###

data$O2_concentration<-data$O2_concentration/1000000*100

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  dplyr::select(Species, Tree.No, Tissue, O2_concentration) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(O2_concentration)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(O2_concentration = as.numeric(O2_concentration)) %>%  # Convert O2_concentration to numeric
  filter(!is.na(O2_concentration)) %>%  # Remove rows where O2_concentration couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(O2_concentration = mean(O2_concentration), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = O2_concentration) %>%
  filter(!is.na(heartwood) & !is.na(sapwood))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, (heartwood)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = heartwood, xend = sapwood, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = heartwood, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = sapwood, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
  scale_color_manual(
    values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),  # Custom colors for tissues
    name = "Tissue"  # Legend title
  ) +
  facet_wrap(~ Species, scales = "free") +  # Free axes for each species
  labs(
    #title = "O2 Concentrations by Tree and Tissue Type",
    x = "O2 (%)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),  # Customize facet label size
    axis.text = element_text(size = 10),   # Customize axis text size
    axis.title = element_text(size = 12)   # Customize axis title size
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





###N2O###

#data$N2O_concentration<-log(data$N2O_concentration)

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  dplyr::select(Species, Tree.No, Tissue, N2O_concentration) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(N2O_concentration)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(N2O_concentration = as.numeric(N2O_concentration)) %>%  # Convert N2O_concentration to numeric
  filter(!is.na(N2O_concentration)) %>%  # Remove rows where N2O_concentration couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(N2O_concentration = mean(N2O_concentration), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = N2O_concentration) %>%
  filter(!is.na(heartwood) & !is.na(sapwood))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(heartwood)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = heartwood, xend = sapwood, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = heartwood, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = sapwood, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
  scale_color_manual(
    values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),  # Custom colors for tissues
    name = "Tissue"  # Legend title
  ) +
  facet_wrap(~ Species, scales = "free") +  # Free axes for each species
  labs(
    #title = "N2O Concentrations by Tree and Tissue Type",
    x = "N2O (ppm)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),  # Customize facet label size
    axis.text = element_text(size = 10),   # Customize axis text size
    axis.title = element_text(size = 12)   # Customize axis title size
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





###CO2###

#data$CO2_concentration<-log(data$CO2_concentration)

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  dplyr::select(Species, Tree.No, Tissue, CO2_concentration) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(CO2_concentration)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(CO2_concentration = as.numeric(CO2_concentration)) %>%  # Convert CO2_concentration to numeric
  filter(!is.na(CO2_concentration)) %>%  # Remove rows where CO2_concentration couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(CO2_concentration = mean(CO2_concentration), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = CO2_concentration) %>%
  filter(!is.na(heartwood) & !is.na(sapwood))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(heartwood)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = heartwood, xend = sapwood, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = heartwood, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = sapwood, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
  scale_color_manual(
    values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),  # Custom colors for tissues
    name = "Tissue"  # Legend title
  ) +
  facet_wrap(~ Species, scales = "free") +  # Free axes for each species
  labs(
    #title = "CO2 Concentrations by Tree and Tissue Type",
    x = "CO2 (ppm)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),  # Customize facet label size
    axis.text = element_text(size = 10),   # Customize axis text size
    axis.title = element_text(size = 12)   # Customize axis title size
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




df_wide <- df_wide %>%
  # Arrange data by species and descending Tree.No to set the new order
  arrange(Species, (Tree.No)) %>%
  # Group by Species for renumbering within each species group
  group_by(Species) %>%
  # Create new Tree.No by combining spcode with descending sequence numbers
  mutate(Tree.No = paste0(row_number())) %>%
  ungroup() %>%
  # Reorder Tree.No as a factor based on the new order in the dataframe
  mutate(Tree.No = fct_reorder(Tree.No, row_number(), .desc = TRUE))

# Set Species as a factor with levels in alphabetical order
df_wide$Species <- factor(df_wide$Species, 
                          levels = c("Acer rubrum", "Acer saccharum", "Betula lenta", 
                                     "Pinus strobus", "Quercus rubra", "Tsuga canadensis"))

# Function to create base plot with non-log values and an inset with log-transformed values, plus species annotation
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
    #annotate("text", x = Inf, y = Inf, label = species_name, hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold") +
    labs(title = species_name)+
    theme_classic() +
    theme(
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(size = 10, hjust = 0),
      axis.title.x = element_text(size = 8)  # Set x-axis label font size
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
