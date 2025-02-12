# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the original data from CSV
#data <- read.csv("Downloads/GC_241010 - Run 1 Compiled (1).csv", stringsAsFactors = FALSE)
data <- GC_data

# Ensure all concentration columns are numeric and apply pmax to set minimum to 0
#correct sample concentration to tree concentration (samples were diluted 1:3)
data$CH4_Area <- pmax(0, data$CH4_concentration * 3)  # Multiplies non-negative values by 3, sets negatives to 0
data$O2_Area  <- pmax(0, data$O2_concentration * 3)
data$O2_Area  <- pmin(210000, data$O2_Area)
data$CO2_Area <- pmax(0, data$CO2_concentration * 3)
data$N2O_Area <- pmax(0, data$N2O_concentration * 3)


#data$CH4_Area<-log(data$CH4_Area)

# dplyr::recode the species based on the mapping provided
species_mapping <- c("BB" = "BELE", "H" = "TSCA", "RM" = "ACRU", "RO" = "QURU", "SM" = "ACSA", "WP" = "PIST")
data$Species <- dplyr::recode(data$Species, !!!species_mapping)


###CH4###

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  dplyr::select(Species, Tree.No, Tissue, CH4_Area) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(CH4_Area)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(CH4_Area = as.numeric(CH4_Area)) %>%  # Convert CH4_Area to numeric
  filter(!is.na(CH4_Area)) %>%  # Remove rows where CH4_Area couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(CH4_Area = mean(CH4_Area), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = CH4_Area) %>%
  filter(!is.na(H) & !is.na(S))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(H)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = H, xend = S, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = H, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = S, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
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
  dplyr::select(Species, Tree.No, Tissue, CH4_Area) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(CH4_Area)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(CH4_Area = as.numeric(CH4_Area)) %>%  # Convert CH4_Area to numeric
  filter(!is.na(CH4_Area)) %>%  # Remove rows where CH4_Area couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(CH4_Area = mean(CH4_Area), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = CH4_Area) %>%
  filter(!is.na(H) & !is.na(S))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(H)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

df_wide$H<-log10(df_wide$H)
df_wide$S<-log10(df_wide$S)

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = H, xend = S, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = H, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = S, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
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

#data$O2_Area<-log(data$O2_Area)

# dplyr::recode the species based on the mapping provided
species_mapping <- c("BB" = "BELE", "H" = "TSCA", "RM" = "ACRU", "RO" = "QURU", "SM" = "ACSA", "WP" = "PIST")
data$Species <- dplyr::recode(data$Species, !!!species_mapping)

data$O2_Area<-data$O2_Area/1000000*100

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  dplyr::select(Species, Tree.No, Tissue, O2_Area) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(O2_Area)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(O2_Area = as.numeric(O2_Area)) %>%  # Convert O2_Area to numeric
  filter(!is.na(O2_Area)) %>%  # Remove rows where O2_Area couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(O2_Area = mean(O2_Area), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = O2_Area) %>%
  filter(!is.na(H) & !is.na(S))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, (H)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = H, xend = S, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = H, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = S, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
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

#data$N2O_Area<-log(data$N2O_Area)

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  dplyr::select(Species, Tree.No, Tissue, N2O_Area) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(N2O_Area)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(N2O_Area = as.numeric(N2O_Area)) %>%  # Convert N2O_Area to numeric
  filter(!is.na(N2O_Area)) %>%  # Remove rows where N2O_Area couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(N2O_Area = mean(N2O_Area), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = N2O_Area) %>%
  filter(!is.na(H) & !is.na(S))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(H)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = H, xend = S, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = H, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = S, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
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

#data$CO2_Area<-log(data$CO2_Area)

# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  dplyr::select(Species, Tree.No, Tissue, CO2_Area) %>%  # Adjust column names if necessary
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(CO2_Area)) %>%  # Remove rows with missing values
  filter(Tissue != "") %>%  # Ensure no empty strings in Tissue column
  mutate(CO2_Area = as.numeric(CO2_Area)) %>%  # Convert CO2_Area to numeric
  filter(!is.na(CO2_Area)) %>%  # Remove rows where CO2_Area couldn't be converted to numeric
  group_by(Species, Tree.No, Tissue) %>%  # Group by the key columns
  summarize(CO2_Area = mean(CO2_Area), .groups = 'drop')  # Average duplicates if they exist

# Reshape the data to have separate columns for heartwood and sapwood
df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = CO2_Area) %>%
  filter(!is.na(H) & !is.na(S))  # Keep only rows where both heartwood and sapwood values are present

# Reorder Tree.No based on heartwood concentration (H), in descending order
df_wide <- df_wide %>%
  arrange(Species, desc(H)) %>%  # Sort by Species and then by descending H
  mutate(Tree.No = factor(Tree.No, levels = rev(unique(Tree.No))))  # Reorder Tree.No factor in descending order

# Create the dumbbell plot using ggplot2
ggplot(df_wide, aes(y = Tree.No)) +
  geom_segment(aes(x = H, xend = S, yend = Tree.No), color = "gray") +  # Draw lines between H and S
  geom_point(aes(x = H, color = "Heartwood"), size = 3, alpha = 0.7) +  # Points for heartwood
  geom_point(aes(x = S, color = "Sapwood"), size = 3, alpha = 0.7) +   # Points for sapwood
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

