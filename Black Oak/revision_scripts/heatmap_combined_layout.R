
# Create combined plot
combined_plot <- p_16S + p_ITS + 
  plot_layout(ncol = 1) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

# Display combined plot
combined_plot

# Save combined plot
ggsave("combined_heatmaps.tif", 
       plot = combined_plot,
       height = 297,
       width = 210,
       units = "mm",
       dpi = 600)
