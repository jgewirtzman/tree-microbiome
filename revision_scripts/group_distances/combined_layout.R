# Load required libraries
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("ggplotify", quietly = TRUE)) install.packages("ggplotify")
library(patchwork)
library(ggplotify)

# MINIMAL ADDITION 3: Create the combined 2x2 figure
# Create text labels for barplots
text_16S <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "16S", size = 6, fontface = "bold") +
  theme_void()

text_ITS <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "ITS", size = 6, fontface = "bold") +
  theme_void()

# Create barplot columns with text above
barplot_col_16S <- text_16S / barplot_16S_wunifrac + 
  plot_layout(heights = c(0.15, 1))

barplot_col_ITS <- text_ITS / barplot_ITS_wunifrac + 
  plot_layout(heights = c(0.15, 1))

# Combine everything in a 2x2 grid
combined_figure <- wrap_plots(
  barplot_col_16S, as.ggplot(heatmap_16S_wunifrac),
  barplot_col_ITS, as.ggplot(heatmap_ITS_wunifrac),
  ncol = 2, nrow = 2,
  widths = c(1, 2)
)

# Display the combined figure
print(combined_figure)

# Save the figure
ggsave("Combined_16S_ITS_WUnifrac_Figure.tif", 
       combined_figure, 
       width = 12, 
       height = 12, 
       units = "in",
       dpi = 600)
