# NEW CODE FROM SCRATCH - 8 Panel Comprehensive Plot
# Using unique object names to avoid conflicts with existing code

# Load required libraries (if not already loaded)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

# Define colors (redefine to be safe)
CORE_COLORS_NEW <- c("Heartwood" = "#a6611a", "Sapwood" = "#dfc27d", 
                     "Mineral" = "#80cdc1", "Organic" = "#018571")

###########################################
# TISSUE TYPE PLOTS (from original code) #
###########################################

# These should already exist from your original code:
# tissue_16s_chao1 = p3 (16S Chao1 by tissue type)
# tissue_16s_shannon = p4 (16S Shannon by tissue type) 
# tissue_its_chao1 = p5 (ITS Chao1 by tissue type)
# tissue_its_shannon = p6 (ITS Shannon by tissue type)

# Assign to new names to be safe and add panel letters, bigger text
tissue_16s_chao1_plot <- p3 + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, label = "a", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

tissue_16s_shannon_plot <- p4 + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, label = "e", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

tissue_its_chao1_plot <- p5 + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, label = "b", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

tissue_its_shannon_plot <- p6 + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, label = "f", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

############################################
# SPECIES PLOTS (recreate individually)   #
############################################

# Create individual species plots with unique names (remove titles, add panel letters)
species_16s_chao1_plot <- ggplot(filter(alpha_div_long, measure == "Chao1"), 
                                 aes(x = species, y = value, color = core_type)) +
  geom_point(size = 2, alpha = 0.5, position = position_identity()) +
  stat_summary(aes(group = core_type),
               fun = mean, 
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               position = position_identity(),
               size = 1, alpha = 0.7) +
  scale_color_manual(values = CORE_COLORS_NEW) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  labs(x = NULL, y = "Chao1") +
  annotate("text", x = -Inf, y = Inf, label = "c", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

species_16s_shannon_plot <- ggplot(filter(alpha_div_long, measure == "Shannon"), 
                                   aes(x = species, y = value, color = core_type)) +
  geom_point(size = 2, alpha = 0.5, position = position_identity()) +
  stat_summary(aes(group = core_type),
               fun = mean, 
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               position = position_identity(),
               size = 1, alpha = 0.7) +
  scale_color_manual(values = CORE_COLORS_NEW) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  labs(x = NULL, y = "Shannon") +
  annotate("text", x = -Inf, y = Inf, label = "g", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

species_its_chao1_plot <- ggplot(filter(alpha_div_its_long, measure == "Chao1"), 
                                 aes(x = species, y = value, color = core_type)) +
  geom_point(size = 2, alpha = 0.5, position = position_identity()) +
  stat_summary(aes(group = core_type),
               fun = mean, 
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               position = position_identity(),
               size = 1, alpha = 0.7) +
  scale_color_manual(values = CORE_COLORS_NEW) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  labs(x = NULL, y = "Chao1") +
  annotate("text", x = -Inf, y = Inf, label = "d", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

# Last plot with legend and smaller legend points
species_its_shannon_plot <- ggplot(filter(alpha_div_its_long, measure == "Shannon"), 
                                   aes(x = species, y = value, color = core_type)) +
  geom_point(size = 2, alpha = 0.5, position = position_identity()) +
  stat_summary(aes(group = core_type),
               fun = mean, 
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               position = position_identity(),
               size = 1, alpha = 0.7) +
  scale_color_manual(values = CORE_COLORS_NEW) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = NULL, y = "Shannon", color = "Sample Type") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
  annotate("text", x = -Inf, y = Inf, label = "h", hjust = -0.5, vjust = 1.5, 
           size = 6, fontface = "bold")

########################################
# CREATE 8-PANEL COMPREHENSIVE PLOT   #
########################################

# Create the 8-panel layout: 4 rows x 2 columns
# Left column: All 16S | Right column: All ITS
# Top to bottom: Tissue Chao1, Species Chao1, Tissue Shannon, Species Shannon
eight_panel_plot <- (tissue_16s_chao1_plot | tissue_its_chao1_plot) / 
  (species_16s_chao1_plot | species_its_chao1_plot) / 
  (tissue_16s_shannon_plot | tissue_its_shannon_plot) / 
  (species_16s_shannon_plot | species_its_shannon_plot) +
  plot_layout(heights = c(1, 1, 1, 1.2))  # Extra height for bottom row with legend

# Display the plot
print(eight_panel_plot)

# Save as A4 size .tif (210 x 297 mm = 8.27 x 11.69 inches)
ggsave("comprehensive_alpha_diversity_8panels.tif", 
       plot = eight_panel_plot,
       width = 8.27, 
       height = 11.69, 
       units = "in",
       dpi = 300,
       compression = "lzw")

# Print confirmation
cat("\n8-panel comprehensive plot created and saved!\n")
cat("File: comprehensive_alpha_diversity_8panels.tif\n")
cat("Layout:\n")
cat("Row 1: Tissue Chao1 (16S | ITS)\n")
cat("Row 2: Species Chao1 (16S | ITS)\n")
cat("Row 3: Tissue Shannon (16S | ITS)\n")
cat("Row 4: Species Shannon (16S | ITS + legend)\n")