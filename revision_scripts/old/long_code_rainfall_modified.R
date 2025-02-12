# Function to unload all non-base packages
unload_all_packages <- function() {
  # Get the list of currently loaded packages
  loaded_packages <- setdiff(loadedNamespaces(), c("base", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base"))
  
  # Detach each package
  for (pkg in loaded_packages) {
    try(detach(paste("package:", pkg, sep = ""), character.only = TRUE, unload = TRUE), silent = TRUE)
  }
}

# Call the function to unload all non-base packages
unload_all_packages()


# Load necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(broom)
library(cowplot)

# Load the data files
O2_standards <- read.csv("Downloads/GC_241010 - O2 Standards.csv")
GHG_standards <- read.csv("Downloads/GC_241010 - GHG Standards.csv")
GC_data <- read.csv("Downloads/GC_241010 - Run 1 Compiled (6).csv")

GC_data$CH4_Area[331]<-NA
GC_data$O2_Area[which(GC_data$Sample.ID=="N2")]<-NA

# Rename columns for easier processing (remove spaces)
colnames(GC_data) <- make.names(colnames(GC_data))

# Merge the O2 and GHG standards with the GC data based on Sample/Sample ID
# O2 Standards Merge
O2_data <- O2_standards %>%
  left_join(GC_data, by = c("Sample" = "Sample.ID"))

# GHG Standards Merge
GHG_data <- GHG_standards %>%
  left_join(GC_data, by = c("Sample" = "Sample.ID"))

# Convert relevant columns to numeric to avoid coercion issues
# Remove commas from concentration columns and convert to numeric
O2_data$O2_Area <- as.numeric(O2_data$O2_Area)
O2_data$O2_ppm <- as.numeric(gsub(",", "", O2_data$X.O2...ppm.))  # Convert O2 concentration to numeric

GHG_data$N2O_Area <- as.numeric(GHG_data$N2O_Area)
GHG_data$N2O_ppm <- as.numeric(gsub(",", "", GHG_data$X.N2O...ppm.))  # Convert N2O concentration to numeric

GHG_data$CH4_Area <- as.numeric(GHG_data$CH4_Area)
GHG_data$CH4_ppm <- as.numeric(gsub(",", "", GHG_data$X.CH4...ppm.))  # Convert CH4 concentration to numeric

GHG_data$CO2_Area <- as.numeric(GHG_data$CO2_Area)
GHG_data$CO2_ppm <- as.numeric(gsub(",", "", GHG_data$X.CO2...ppm.))  # Convert CO2 concentration to numeric

# Filter out any NAs or missing values
GHG_data <- GHG_data %>%
  filter(!is.na(N2O_Area) & !is.na(N2O_ppm)) %>%
  filter(!is.na(CH4_Area) & !is.na(CH4_ppm)) %>%
  filter(!is.na(CO2_Area) & !is.na(CO2_ppm))

# Define the low-range and high-range standards for O2
low_range_o2_standards <- c("O1", "O1.5", "O2", "O2.5")
high_range_o2_standards <- c( "O3", "O4", "O5", "O6")

# Filter O2 data to separate low-range and high-range standards
low_range_o2_data <- O2_data %>%
  filter(Sample %in% low_range_o2_standards)

high_range_o2_data <- O2_data %>%
  filter(Sample %in% high_range_o2_standards)

# Combine both low and high range standards for high-range curve
all_high_range_o2_data <- bind_rows(low_range_o2_data, high_range_o2_data)

# Create low-range and high-range standard curves for O2
O2_low_curve <- lm(O2_ppm ~ O2_Area, data = low_range_o2_data)
O2_high_curve <- lm(O2_ppm ~ O2_Area, data = all_high_range_o2_data)

# Get the maximum peak area for the low-range O2 curve
max_O2_low_area <- max(low_range_o2_data$O2_Area, na.rm = TRUE)

# Define the low-range and high-range standards for GHGs
low_range_standards <- c("SB1", "SB3a", "SB3b", "SB3", "SB4a")
high_range_standards <- c("SB4", "SB5a", "SB5", "S3a", "S3b", "S3c")

# Filter GHG data to separate low-range and high-range standards
low_range_data <- GHG_data %>%
  filter(Sample %in% low_range_standards)

high_range_data <- GHG_data %>%
  filter(Sample %in% high_range_standards)

# Create low-range and high-range standard curves for low-range and all data for high-range
N2O_low_curve <- lm(N2O_ppm ~ N2O_Area, data = low_range_data)
N2O_high_curve <- lm(N2O_ppm ~ N2O_Area, data = GHG_data)

CH4_low_curve <- lm(CH4_ppm ~ CH4_Area, data = low_range_data)
CH4_high_curve <- lm(CH4_ppm ~ CH4_Area, data = GHG_data)

CO2_low_curve <- lm(CO2_ppm ~ CO2_Area, data = low_range_data)
CO2_high_curve <- lm(CO2_ppm ~ CO2_Area, data = GHG_data)

# Get the maximum peak areas for the low-range curves
max_N2O_low_area <- max(low_range_data$N2O_Area, na.rm = TRUE)
max_CH4_low_area <- max(low_range_data$CH4_Area, na.rm = TRUE)
max_CO2_low_area <- max(low_range_data$CO2_Area, na.rm = TRUE)

# Function to calculate concentration based on peak area and curve selection
calculate_concentration <- function(peak_area, max_low_area, low_curve, high_curve) {
  if (peak_area <= max_low_area) {
    return(predict(low_curve, newdata = data.frame(peak_area = peak_area)))
  } else {
    return(predict(high_curve, newdata = data.frame(peak_area = peak_area)))
  }
}

# Apply the calculation without rowwise(), using if_else and vectorized operations
GC_data <- GC_data %>%
  mutate(
    O2_concentration = if_else(
      O2_Area <= max_O2_low_area,
      predict(O2_low_curve, newdata = data.frame(O2_Area = O2_Area)),
      predict(O2_high_curve, newdata = data.frame(O2_Area = O2_Area))
    ),
    
    N2O_concentration = if_else(
      N2O_Area <= max_N2O_low_area,
      predict(N2O_low_curve, newdata = data.frame(N2O_Area = N2O_Area)),
      predict(N2O_high_curve, newdata = data.frame(N2O_Area = N2O_Area))
    ),
    CH4_concentration = if_else(
      CH4_Area <= max_CH4_low_area,
      predict(CH4_low_curve, newdata = data.frame(CH4_Area = CH4_Area)),
      predict(CH4_high_curve, newdata = data.frame(CH4_Area = CH4_Area))
    ),
    CO2_concentration = if_else(
      CO2_Area <= max_CO2_low_area,
      predict(CO2_low_curve, newdata = data.frame(CO2_Area = CO2_Area)),
      predict(CO2_high_curve, newdata = data.frame(CO2_Area = CO2_Area))
    )
  )


# Function to plot standard curve with equation and R^2
plot_standard_curve <- function(low_data, high_data, low_curve, high_curve, analyte, x, y, range_type) {
  low_summary <- summary(low_curve)
  high_summary <- summary(high_curve)
  
  low_equation <- paste("Low:", round(coef(low_curve)[2], 3), "*x +", round(coef(low_curve)[1], 3))
  high_equation <- paste("High:", round(coef(high_curve)[2], 3), "*x +", round(coef(high_curve)[1], 3))
  
  low_r2 <- paste("R² (low) =", round(low_summary$r.squared, 3))
  high_r2 <- paste("R² (high) =", round(high_summary$r.squared, 3))
  
  ggplot() +
    geom_point(data = high_data, aes_string(x = x, y = y), color = "red") +
    geom_smooth(data = high_data, aes_string(x = x, y = y), method = "lm", se = FALSE, color = "red") +
    geom_point(data = low_data, aes_string(x = x, y = y), color = "blue") +
    geom_smooth(data = low_data, aes_string(x = x, y = y), method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, label = paste(low_equation, low_r2, sep = "\n"), 
             hjust = 1.2, vjust = 2, size = 3.5, color = "blue", parse = FALSE) +
    annotate("text", x = Inf, y = Inf, label = paste(high_equation, high_r2, sep = "\n"), 
             hjust = 1.2, vjust = 1, size = 3.5, color = "red", parse = FALSE) +
    labs(title = paste(analyte, range_type, "Standard Curve"),
         x = "Peak Area",
         y = "Concentration (ppm)") +
    theme_minimal()
}

# Create separate plots for low-range and high-range standards for each analyte
# O2 Low-Range and High-Range Curves
O2_low_plot <- plot_standard_curve(low_range_o2_data, O2_data, O2_low_curve, O2_high_curve, "O2", "O2_Area", "O2_ppm", "Low-Range")
O2_high_plot <- plot_standard_curve(all_high_range_o2_data, O2_data, O2_low_curve, O2_high_curve, "O2", "O2_Area", "O2_ppm", "High-Range")

# Low-Range Curves for GHGs
N2O_low_plot <- plot_standard_curve(low_range_data, GHG_data, N2O_low_curve, N2O_high_curve, "N2O", "N2O_Area", "N2O_ppm", "Low-Range")
CH4_low_plot <- plot_standard_curve(low_range_data, GHG_data, CH4_low_curve, CH4_high_curve, "CH4", "CH4_Area", "CH4_ppm", "Low-Range")
CO2_low_plot <- plot_standard_curve(low_range_data, GHG_data, CO2_low_curve, CO2_high_curve, "CO2", "CO2_Area", "CO2_ppm", "Low-Range")

# High-Range Curves for GHGs (including all standards)
N2O_high_plot <- plot_standard_curve(GHG_data, GHG_data, N2O_low_curve, N2O_high_curve, "N2O", "N2O_Area", "N2O_ppm", "High-Range")
CH4_high_plot <- plot_standard_curve(GHG_data, GHG_data, CH4_low_curve, CH4_high_curve, "CH4", "CH4_Area", "CH4_ppm", "High-Range")
CO2_high_plot <- plot_standard_curve(GHG_data, GHG_data, CO2_low_curve, CO2_high_curve, "CO2", "CO2_Area", "CO2_ppm", "High-Range")

# Display the plots in RStudio
print(O2_low_plot)
print(O2_high_plot)

print(N2O_low_plot)
print(N2O_high_plot)

print(CH4_low_plot)
print(CH4_high_plot)

print(CO2_low_plot)
print(CO2_high_plot)


# Function to plot standard curve with equation, R^2, and axis limits for low-range
plot_low_range_curve <- function(low_data, low_curve, analyte, x, y) {
  low_summary <- summary(low_curve)
  low_equation <- paste("Low:", round(coef(low_curve)[2], 3), "*x +", round(coef(low_curve)[1], 3))
  low_r2 <- paste("R² (low) =", round(low_summary$r.squared, 3))
  
  # Define the x-axis limits based on the low-range data
  x_min <- min(low_data[[x]], na.rm = TRUE)
  x_max <- max(low_data[[x]], na.rm = TRUE)
  
  ggplot() +
    geom_point(data = low_data, aes_string(x = x, y = y), color = "blue") +
    geom_smooth(data = low_data, aes_string(x = x, y = y), method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, label = paste(low_equation, low_r2, sep = "\n"), 
             hjust = 1.2, vjust = 2, size = 3.5, color = "blue", parse = FALSE) +
    scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits for low-range
    labs(title = paste(analyte, "Standard Curve (Low-Range)"),
         x = "Peak Area",
         y = "Concentration (ppm)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
}

# Low-Range Plots for O2, N2O, CH4, and CO2
O2_low_plot <- plot_low_range_curve(low_range_o2_data, O2_low_curve, "O2", "O2_Area", "O2_ppm")
N2O_low_plot <- plot_low_range_curve(low_range_data, N2O_low_curve, "N2O", "N2O_Area", "N2O_ppm")
CH4_low_plot <- plot_low_range_curve(low_range_data, CH4_low_curve, "CH4", "CH4_Area", "CH4_ppm")
CO2_low_plot <- plot_low_range_curve(low_range_data, CO2_low_curve, "CO2", "CO2_Area", "CO2_ppm")

O2_low_plot
N2O_low_plot
CH4_low_plot
CO2_low_plot




# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the original data from CSV
#data <- read.csv("Downloads/GC_241010 - Run 1 Compiled (1).csv", stringsAsFactors = FALSE)
data <- GC_data


# Recode the species based on the mapping provided
species_mapping <- c("BB" = "BELE", "H" = "TSCA", "RM" = "ACRU", "RO" = "QURU", "SM" = "ACSA", "WP" = "PIST")
data <- data %>%   mutate(Species = species_mapping[match(Species, names(species_mapping))])

tissue_mapping <- c("S" = "sapwood", "H" = "heartwood")
data <- data %>%   mutate(Tissue = tissue_mapping[match(Tissue, names(tissue_mapping))])


# Ensure all concentration columns are numeric and apply pmax to set minimum to 0
#correct sample concentration to tree concentration (samples were diluted 1:3)
data$CH4_Area <- pmax(0, data$CH4_concentration * 3)  # Multiplies non-negative values by 3, sets negatives to 0
data$O2_Area  <- pmax(0, data$O2_concentration * 3)
data$O2_Area  <- pmin(210000, data$O2_Area)
data$CO2_Area <- pmax(0, data$CO2_concentration * 3)
data$N2O_Area <- pmax(0, data$N2O_concentration * 3)

data$O2_Area<-data$O2_Area/1000000*100
data$CH4_Area<-data$CH4_Area/1000000*100


# Filter the data to include only relevant columns and rows
df_filtered <- data %>%
  filter(Tissue != "") 


density_plot <- ggplot(df_filtered, aes(x=O2_Area, fill=Tissue))+
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4"))  + # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+  ggtitle("pmoA")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(df_filtered, aes(x = O2_Area, y = Tissue, fill = Tissue)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) + 
  # Custom colors and labels
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab(expression(O[2] ~ "(%)"))


aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

PA<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))

df_filtered$log10ch4<-log10(df_filtered$CH4_Area*1000000/100)

density_plot <- ggplot(df_filtered, aes(x=log10ch4, fill=Tissue))+
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4"))  + # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+  ggtitle("pmoA")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(df_filtered, aes(x =log10ch4, y = Tissue, fill = Tissue)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) + 
  # Custom colors and labels
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab(expression(log[10] ~ "(" * CH[4] ~ "ppm)"))


aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

PB<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))

tissue_plot<-cowplot::plot_grid(PA, PB)










library(tidyverse)


###DDPCR PLOT###

ddpcr<-read.csv("/Users/jongewirtzman/Downloads/ddPCR_tree_transposed_data.csv")

ddpcr <- ddpcr %>%
  mutate(core_type = dplyr::recode(core_type, "Inner" = "heartwood", "Outer" = "sapwood"))

# Overlapping density plot
density_plot <- ggplot(ddpcr, aes(x = log10(pmoa_loose), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+  ggtitle("pmoA")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(ddpcr, aes(x = log10(pmoa_loose), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) + 
  # Custom colors and labels
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab(expression(log[10] ~ "(pmoA)"))

# Combine the plots vertically using cowplot
library(cowplot)

# Align the density plot on top and the boxplot with jittered points below
combined_plot1 <- cowplot::plot_grid(density_plot, box_jitter_plot, ncol = 1, align = "v", rel_heights = c(1, 1.5))
#combined_plot1


aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

combined_plot1<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))


###




# Overlapping density plot
density_plot <- ggplot(ddpcr, aes(x = log10(mcra_probe_loose), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+  ggtitle("mcrA")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(ddpcr, aes(x = log10(mcra_probe_loose), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  #  scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab(expression(log[10] ~ "(mcrA)"))

# Align the density plot on top and the boxplot with jittered points below
combined_plot2 <- cowplot::plot_grid(density_plot, box_jitter_plot, ncol = 1, align = "v", rel_heights = c(1, 1.5))
#combined_plot2

aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

combined_plot2<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))





###FAPROTAX 1###
faprotax<-read.csv("/Users/jongewirtzman/Downloads/metabolisms/16S_metabolisms_weighted.csv")

names(faprotax)[1]<-"SampleID"

# Extract the sample type (Inner, Outer, Organic, Mineral) from the sample name
faprotax <- faprotax %>%
  mutate(Sample_Type = case_when(
    grepl("Inner", SampleID) ~ "Inner",
    grepl("Outer", SampleID) ~ "Outer",
    grepl("Organic", SampleID) ~ "Organic",
    grepl("Mineral", SampleID) ~ "Mineral",
    TRUE ~ "Unknown"
  ))

faprotax <- faprotax %>%
  filter(!(Sample_Type %in% c("Organic", "Mineral")))


faprotax <- faprotax %>%
  mutate(core_type = dplyr::recode(Sample_Type, "Inner" = "heartwood", "Outer" = "sapwood"))

# Overlapping density plot
density_plot <- ggplot(faprotax, aes(x = (aerobic_chemoheterotrophy), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+ggtitle("Aerobic Chemoheterotrophy")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(faprotax, aes(x = aerobic_chemoheterotrophy, y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) + 
  # Custom colors and labels
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab("% aerobic chemoheterotrophy")

# Combine the plots vertically using cowplot
library(cowplot)

# Align the density plot on top and the boxplot with jittered points below
combined_plot3 <- cowplot::plot_grid(density_plot, box_jitter_plot, ncol = 1, align = "v", rel_heights = c(1, 1.5))
#combined_plot1


aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

combined_plot3<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))


###




# Overlapping density plot
density_plot <- ggplot(faprotax, aes(x = fermentation, fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity', bw = 5) +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+  ggtitle("Fermentation")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(faprotax, aes(x = fermentation, y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  #  scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab("% fermentation")

# Align the density plot on top and the boxplot with jittered points below
combined_plot4 <- cowplot::plot_grid(density_plot, box_jitter_plot, ncol = 1, align = "v", rel_heights = c(1, 1.5))
#combined_plot2

aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

combined_plot4<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))




###FAPROTAX 2###
faprotax<-read.csv("/Users/jongewirtzman/Downloads/metabolisms/16S_metabolisms_weighted.csv")

names(faprotax)[1]<-"SampleID"

# Extract the sample type (Inner, Outer, Organic, Mineral) from the sample name
faprotax <- faprotax %>%
  mutate(Sample_Type = case_when(
    grepl("Inner", SampleID) ~ "Inner",
    grepl("Outer", SampleID) ~ "Outer",
    grepl("Organic", SampleID) ~ "Organic",
    grepl("Mineral", SampleID) ~ "Mineral",
    TRUE ~ "Unknown"
  ))

faprotax <- faprotax %>%
  filter(!(Sample_Type %in% c("Organic", "Mineral")))


faprotax <- faprotax %>%
  mutate(core_type = dplyr::recode(Sample_Type, "Inner" = "heartwood", "Outer" = "sapwood"))

# Overlapping density plot
density_plot <- ggplot(faprotax, aes(x = log10(hydrogenotrophic_methanogenesis), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+  ggtitle("Hyderogenotrophic Methanogenesis")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(faprotax, aes(x = log10(hydrogenotrophic_methanogenesis), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) + 
  # Custom colors and labels
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab(expression(log[10] ~ "(% " * H[2] * " methanogenesis)"))

# Combine the plots vertically using cowplot
library(cowplot)

# Align the density plot on top and the boxplot with jittered points below
combined_plot5 <- cowplot::plot_grid(density_plot, box_jitter_plot, ncol = 1, align = "v", rel_heights = c(1, 1.5))
#combined_plot1


aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

combined_plot5<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))


###




# Overlapping density plot
density_plot <- ggplot(faprotax, aes(x = log10(methanotrophy), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +  # Overlapping density plots
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  labs(fill = NULL) + # Custom colors
  #scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(10, 5, 0, 5)) 
#+  ggtitle("Methanotrophy")

# Boxplot and jitter plot
box_jitter_plot <- ggplot(faprotax, aes(x = log10(methanotrophy), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width=0, alpha = 0.4, size = 1) +  # Jittered points
  geom_boxplot(width = 0.4, alpha = 0.7) +  # Boxplot
  scale_fill_manual(values = c("#a6611a", "#1f78b4")) +  # Custom colors and labels
  #  scale_x_continuous(limits = c(170, 230)) +  # Set x-axis limits
  #facet_grid(rows = vars(core_type), scales = "free_y") +  # Faceting similar to Vega-Lite
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )+ 
  theme(plot.margin = margin(0, 5, 10, 5)) +
  xlab(expression(log[10] ~ "(% methanotrophy)"))

# Align the density plot on top and the boxplot with jittered points below
combined_plot6 <- cowplot::plot_grid(density_plot, box_jitter_plot, ncol = 1, align = "v", rel_heights = c(1, 1.5))
#combined_plot2

aligned<-align_plots(density_plot, box_jitter_plot, align = "v", axis="tblr")

p1<-ggdraw(aligned[[1]])
p2<-ggdraw(aligned[[2]])

combined_plot6<-cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1,2))




###

#final_plot<-cowplot::plot_grid(combined_plot1, combined_plot2)
#final_plot

rainfall<-cowplot::plot_grid(combined_plot1, combined_plot2, 
          combined_plot6, combined_plot5, 
          combined_plot3, combined_plot4, 
          PA, PB,
          nrow=4)

rainfall

ggsave("rainfall_plot.png", rainfall, width = 5.33, height = 8)

#ggsave("final_plot.png", plot = final_plot, width = 10, height = 5, units = "in", dpi = 300)




