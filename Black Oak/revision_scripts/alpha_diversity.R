# Load required libraries
library(tidyverse)
library(vegan)       # For diversity indices
library(data.table)  # For efficient reshaping
library(patchwork)   # For combining plots
library(viridis)     # For modern color palettes

# -----------------------------
# A. Determine Media Order and Create Global Color Mapping
# -----------------------------
# 1. Read the 16S OTU table
sixteen_s_otu_table <- read.table(
  "Black Oak/16S/OTU_table.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 2. Filter out Mitochondria/Chloroplast
sixteen_s_otu_table <- sixteen_s_otu_table %>%
  filter(!grepl("Mitochondria|Chloroplast", Class, ignore.case = TRUE))

# 3. Pivot to long format, label Media, remove zero abundances
sixteen_s_long <- sixteen_s_otu_table %>%
  pivot_longer(cols = where(is.numeric), names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Media = case_when(
      str_detect(Sample, "HEART")   ~ "Heartwood",
      str_detect(Sample, "SAP")     ~ "Sapwood",
      str_detect(Sample, "LITTER")  ~ "Leaf Litter",
      str_detect(Sample, "BARK")    ~ "Bark",
      str_detect(Sample, "FOLIAGE") ~ "Foliage",
      str_detect(Sample, "MINERAL") ~ "Mineral Soil",
      str_detect(Sample, "ORGANIC") ~ "Organic Soil",
      str_detect(Sample, "BRANCH")  ~ "Branch",
      str_detect(Sample, "COARSE")  ~ "Coarse Root",
      str_detect(Sample, "FINE")    ~ "Fine Root",
      str_detect(Sample, "ROT")     ~ "Rot",
      TRUE                          ~ "Other"
    )
  ) %>%
  filter(Media != "Other") %>%
  filter(Abundance > 0)

# 4. Compute Chao1 per SAMPLE+MEDIA, then average Chao1 by MEDIA
diversity_16S <- sixteen_s_long %>%
  group_by(Sample, Media) %>%
  summarize(
    Chao1 = estimateR(Abundance)[1],
    .groups = "drop"
  ) %>%
  group_by(Media) %>%
  summarize(
    Mean_Chao1 = mean(Chao1, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Sort from highest to lowest Chao1
diversity_16S <- diversity_16S %>%
  arrange(desc(Mean_Chao1))  # Highest on top

# 6. Create global color mapping based on this order
media_order <- diversity_16S$Media
color_mapping <- setNames(viridis::viridis(length(media_order), option = "plasma"), media_order)

cat("Media order (Highest → Lowest 16S Chao1):\n")
print(media_order)

# ------------------------------------------------
# B. Function to process data and create alpha-diversity plots
# ------------------------------------------------
process_data_and_create_plots <- function(otu_table_path, title_prefix, color_mapping) {
  # 1. Read and preprocess
  otu_table <- read.table(otu_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 2. Filter out Mito/Chloro
  otu_table <- otu_table %>%
    filter(!grepl("Mitochondria|Chloroplast", Class, ignore.case = TRUE))
  
  # 3. Pivot to long, define Media, remove zero abundances
  data_cols <- otu_table %>%
    dplyr::select(where(is.numeric)) %>%
    colnames()
  
  otu_table_long <- otu_table %>%
    pivot_longer(cols = all_of(data_cols), names_to = "Sample", values_to = "Abundance") %>%
    mutate(
      Media = case_when(
        str_detect(Sample, "HEART")   ~ "Heartwood",
        str_detect(Sample, "SAP")     ~ "Sapwood",
        str_detect(Sample, "LITTER")  ~ "Leaf Litter",
        str_detect(Sample, "BARK")    ~ "Bark",
        str_detect(Sample, "FOLIAGE") ~ "Foliage",
        str_detect(Sample, "MINERAL") ~ "Mineral Soil",
        str_detect(Sample, "ORGANIC") ~ "Organic Soil",
        str_detect(Sample, "BRANCH")  ~ "Branch",
        str_detect(Sample, "COARSE")  ~ "Coarse Root",
        str_detect(Sample, "FINE")    ~ "Fine Root",
        str_detect(Sample, "ROT")     ~ "Rot",
        TRUE                          ~ "Other"
      )
    ) %>%
    filter(Media != "Other") %>%
    filter(Abundance > 0) %>%
    mutate(OTU_ID = X)  # If needed; otherwise can omit
  
  # 4. Summarize (Sample+Media+OTU)
  otu_table_long <- otu_table_long %>%
    group_by(Sample, Media, OTU_ID) %>%
    summarize(
      Abundance = sum(Abundance),
      .groups   = "drop"
    )
  
  # 5. Calculate alpha diversity (Chao1 + Shannon) at sample level
  diversity_data_otu <- otu_table_long %>%
    group_by(Sample, Media) %>%
    summarize(
      Chao1   = estimateR(Abundance)[1],
      Shannon = diversity(Abundance, index = "shannon"),
      .groups = "drop"
    )
  
  # 6. Summarize means + SE
  summary_data <- diversity_data_otu %>%
    pivot_longer(cols = c(Chao1, Shannon),
                 names_to = "Metric", values_to = "Value") %>%
    group_by(Media, Metric) %>%
    summarize(
      Mean = mean(Value, na.rm = TRUE),
      SE   = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # 7. Dynamically reorder Media within each metric
  summary_data <- summary_data %>%
    group_by(Metric) %>%
    mutate(Media = fct_reorder(Media, Mean, .desc = TRUE)) %>%
    ungroup()
  
  # 8. Build the plot with consistent colors
  alpha_div_plot <- ggplot(summary_data, aes(x = Media, y = Mean, fill = Media)) +
    geom_bar(stat = "identity", alpha = 0.7, color = "black") +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  width = 0.2, size = 0.5) +
    facet_wrap(~ Metric, scales = "free_y") +
    scale_fill_manual(values = color_mapping) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(size = 10, angle = 45, hjust = 1),
      strip.text      = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(
      x     = "Media",
      y     = "Mean Diversity Metric ± SE",
      title = title_prefix
    )
  
  return(alpha_div_plot)
}

# ------------------------------------------------
# C. Create ITS and 16S plots with consistent colors
# ------------------------------------------------

# 1. ITS
its_plots <- process_data_and_create_plots(
  otu_table_path = "Black Oak/ITS/OTU_table.txt",
  title_prefix   = "ITS",
  color_mapping  = color_mapping
)

# 2. 16S
sixteen_s_plots <- process_data_and_create_plots(
  otu_table_path = "Black Oak/16S/OTU_table.txt",
  title_prefix   = "16S",
  color_mapping  = color_mapping
)

# ------------------------------------------------
# D. Combine and display
# ------------------------------------------------
combined_plot <- (sixteen_s_plots / its_plots)
combined_plot

# Save as A4 size .tif (210 x 297 mm = 8.27 x 11.69 inches)
ggsave("quve_alpha_diversity.tif", 
       plot = combined_plot,
       width = 8, 
       height = 8, 
       units = "in",
       dpi = 600)
