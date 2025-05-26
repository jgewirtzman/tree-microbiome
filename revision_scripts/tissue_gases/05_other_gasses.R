library(ggplot2)
library(patchwork)
library(dplyr)
library(forcats)
library(tidyr)
library(cowplot)

# Function to create dumbbell plot for any gas
create_gas_dumbbell_plot <- function(data, gas_column, x_label_expression, 
                                     convert_o2 = FALSE) {
  
  # Convert O2 from ppm to percentage if specified
  if (convert_o2) {
    data[[gas_column]] <- data[[gas_column]] / 1000000 * 100
  }
  
  # Filter the data to include only relevant columns and rows
  df_filtered <- data %>%
    dplyr::select(Species, Tree.No, Tissue, !!sym(gas_column)) %>%
    filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(!!sym(gas_column))) %>%
    filter(Tissue != "") %>%
    mutate(concentration = as.numeric(!!sym(gas_column))) %>%
    filter(!is.na(concentration)) %>%
    group_by(Species, Tree.No, Tissue) %>%
    summarize(concentration = mean(concentration), .groups = 'drop')
  
  # Reshape the data to have separate columns for heartwood and sapwood
  df_wide <- df_filtered %>%
    pivot_wider(names_from = Tissue, values_from = concentration) %>%
    filter(!is.na(heartwood) & !is.na(sapwood))
  
  # For O2, sort by ascending heartwood (lowest O2 first), for others descending
  if (convert_o2) {
    df_wide <- df_wide %>%
      arrange(Species, heartwood) %>%
      group_by(Species) %>%
      mutate(Tree.No = as.character(row_number())) %>%
      ungroup() %>%
      mutate(Tree.No = fct_reorder(Tree.No, heartwood, .desc = FALSE))
  } else {
    df_wide <- df_wide %>%
      arrange(Species, desc(heartwood)) %>%
      group_by(Species) %>%
      mutate(Tree.No = as.character(row_number())) %>%
      ungroup() %>%
      mutate(Tree.No = fct_reorder(Tree.No, heartwood, .desc = TRUE))
  }
  
  # Set Species as a factor with specific order
  df_wide$Species <- factor(df_wide$Species, 
                            levels = c("Acer rubrum", "Acer saccharum", "Betula lenta", 
                                       "Pinus strobus", "Quercus rubra", "Tsuga canadensis"))
  
  # Function to create facet plot
  create_facet_plot <- function(species_data) {
    species_name <- unique(species_data$Species)
    
    ggplot(species_data, aes(y = Tree.No)) +
      geom_segment(aes(x = heartwood, xend = sapwood, yend = Tree.No), color = "grey30") +
      geom_point(aes(x = heartwood, color = "Heartwood"), size = 2.5, alpha = 0.75) +
      geom_point(aes(x = sapwood, color = "Sapwood"), size = 2.5, alpha = 0.75) +
      scale_color_manual(
        values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4"),
        name = "Tissue"
      ) +
      labs(x = x_label_expression, y = "Trees (Ranked)") +
      labs(title = species_name) +
      theme_classic() +
      theme(
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(size = 10, hjust = 0),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)
      ) + coord_flip()
  }
  
  # Split data by species and apply function to create individual facet plots
  species_plots <- lapply(split(df_wide, df_wide$Species), create_facet_plot)
  
  # Arrange plots using wrap_plots with shared legend
  final_plot <- wrap_plots(species_plots, ncol = 2) + 
    plot_layout(guides = "collect") & 
    theme(
      legend.position = "bottom",
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = -10, r = -10, b = -10, l = -10)
    )
  
  return(final_plot)
}

# Create O2 dumbbell plot (convert to percentage)
o2_plot <- create_gas_dumbbell_plot(
  data = GC_data,
  gas_column = "O2_concentration",
  x_label_expression = expression(O[2] * " (%)"),
  convert_o2 = TRUE
)

# Create CO2 dumbbell plot
co2_plot <- create_gas_dumbbell_plot(
  data = GC_data,
  gas_column = "CO2_concentration",
  x_label_expression = expression(CO[2] * " (ppm)"),
  convert_o2 = FALSE
)

# Create N2O dumbbell plot
n2o_plot <- create_gas_dumbbell_plot(
  data = GC_data,
  gas_column = "N2O_concentration",
  x_label_expression = expression(N[2]*O * " (ppm)"),
  convert_o2 = FALSE
)

# Display the plots
print("O2 Dumbbell Plot:")
print(o2_plot)

print("CO2 Dumbbell Plot:")
print(co2_plot)

print("N2O Dumbbell Plot:")
print(n2o_plot)

# Save the plots (PNG only)
ggsave("o2_dumbbell_plot.png", o2_plot, width = 5, height = 5)
ggsave("co2_dumbbell_plot.png", co2_plot, width = 5, height = 5)
ggsave("n2o_dumbbell_plot.png", n2o_plot, width = 5, height = 5)

# Combine all three plots in one row
combined_plot <- plot_grid(o2_plot, co2_plot, n2o_plot, ncol = 2, nrow = 2)

# Save the combined plot
ggsave("all_gases_dumbbell_plots.png", combined_plot, width = 10, height = 10)
