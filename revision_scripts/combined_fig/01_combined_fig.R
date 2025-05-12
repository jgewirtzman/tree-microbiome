library(patchwork)

# Arrange the layout with specified proportions
final_layout <- (rainfall_labeled | (ic_plot / dumbbell_plot + plot_layout(heights = c(1, 2)))) +
  plot_layout(widths = c(1, 1))  # Rainfall takes full left column, right column has 2:1 height ratio

# Display the final layout
print(final_layout)

# Save the layout to a file
ggsave("final_layout.png", final_layout, width = 12, height = 8, dpi = 300)

ggsave("final_layout.svg", final_layout, width = 12, height = 8, dpi = 300)
