analyze_phyloseq_distances <- function(dist_matrix, phyloseq_obj, metric_name) {
  require(phyloseq)
  require(vegan)
  require(lme4)
  require(emmeans)
  require(effectsize)
  require(tidyverse)
  
  # Extract and clean metadata
  metadata_clean <- data.frame(sample_data(phyloseq_obj)) %>%
    filter(!is.na(Species) & !is.na(SampleType))
  
  # Get common samples
  dist_samples <- rownames(as.matrix(dist_matrix))
  metadata_samples <- rownames(metadata_clean)
  common_samples <- intersect(dist_samples, metadata_samples)
  
  print(paste("Total samples in distance matrix:", length(dist_samples)))
  print(paste("Total samples in metadata:", length(metadata_samples)))
  print(paste("Common samples:", length(common_samples)))
  
  # Subset distance matrix to common samples
  dist_matrix_subset <- as.dist(as.matrix(dist_matrix)[common_samples, common_samples])
  metadata_subset <- metadata_clean[common_samples, ]
  
  # PERMANOVA
  permanova_result <- adonis2(dist_matrix_subset ~ SampleType * Species, 
                              data = metadata_subset,
                              permutations = 999)
  print("PERMANOVA Results:")
  print(permanova_result)
  
  # Beta dispersion test
  beta_disp <- betadisper(dist_matrix_subset, metadata_subset$SampleType)
  beta_test <- permutest(beta_disp)
  print("Beta dispersion test results:")
  print(beta_test)
  
  # Convert distance matrix to long format
  dist_df <- as.data.frame(as.matrix(dist_matrix_subset))
  dist_df$Sample1 <- rownames(dist_df)
  
  dist_long <- dist_df %>%
    pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance") %>%
    filter(Sample1 != Sample2) %>%
    left_join(metadata_subset, by = c("Sample1" = "RowName")) %>%
    rename(Species1 = Species, SampleType1 = SampleType) %>%
    left_join(metadata_subset, by = c("Sample2" = "RowName")) %>%
    rename(Species2 = Species, SampleType2 = SampleType) %>%
    mutate(
      Sample1 = factor(Sample1),
      Sample2 = factor(Sample2),
      # Fixed recode syntax
      SampleType1 = case_when(
        SampleType1 == "Inner" ~ "Heartwood",
        SampleType1 == "Outer" ~ "Sapwood",
        TRUE ~ SampleType1
      ),
      SampleType2 = case_when(
        SampleType2 == "Inner" ~ "Heartwood",
        SampleType2 == "Outer" ~ "Sapwood",
        TRUE ~ SampleType2
      ),
      Comparison = case_when(
        Species1 == Species2 & SampleType1 == "Sapwood" & SampleType2 == "Heartwood" ~ "Sapwood_Heartwood_Within",
        Species1 == Species2 & SampleType1 == "Sapwood" & SampleType2 == "Sapwood" ~ "Sapwood_Sapwood_Within",
        Species1 == Species2 & SampleType1 == "Heartwood" & SampleType2 == "Heartwood" ~ "Heartwood_Heartwood_Within",
        Species1 != Species2 & SampleType1 == "Sapwood" & SampleType2 == "Heartwood" ~ "Sapwood_Heartwood_Between",
        Species1 != Species2 & SampleType1 == "Sapwood" & SampleType2 == "Sapwood" ~ "Sapwood_Sapwood_Between",
        Species1 != Species2 & SampleType1 == "Heartwood" & SampleType2 == "Heartwood" ~ "Heartwood_Heartwood_Between",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Comparison)) %>%
    mutate(Comparison = factor(Comparison, levels = c(
      "Sapwood_Heartwood_Within",
      "Sapwood_Sapwood_Within",
      "Heartwood_Heartwood_Within",
      "Sapwood_Heartwood_Between",
      "Sapwood_Sapwood_Between",
      "Heartwood_Heartwood_Between"
    )))
  
  # Rest of the function remains the same...
  # Fit mixed model
  mixed_model <- lmer(Distance ~ Comparison + 
                        (1|Species1) + (1|Species2) + 
                        (1|Sample1) + (1|Sample2),
                      data = dist_long)
  
  # Model diagnostics
  pdf(paste0("model_diagnostics_", metric_name, ".pdf"))
  plot(mixed_model, main = "Residual Plots")
  qqnorm(resid(mixed_model))
  qqline(resid(mixed_model))
  dev.off()
  
  # Calculate effect sizes
  effect_sizes <- eta_squared(mixed_model)
  print("Effect sizes (Eta-squared):")
  print(effect_sizes)
  
  # EMMs and post-hoc
  emm <- emmeans(mixed_model, specs = "Comparison")
  cld <- multcomp::cld(emm, 
                       alpha = 0.05,
                       adjust = "tukey",
                       Letters = letters[1:6])
  
  # Prepare plot data
  plot_data <- as.data.frame(cld) %>%
    dplyr::select(Comparison, emmean, SE, .group) %>%
    rename(Letters = .group) %>%
    left_join(
      dist_long %>%
        group_by(Comparison) %>%
        summarise(
          n = n_distinct(c(Sample1, Sample2))
        )
    )
  
  # Calculate plot dimensions
  y_max <- max(plot_data$emmean + plot_data$SE, na.rm = TRUE)
  bracket_height <- y_max * 0.1
  text_height <- bracket_height * 1.2
  
  # Generate plot
  metric_labels <- c(
    wunifrac = "Weighted UniFrac Distance",
    uunifrac = "Unweighted UniFrac Distance",
    braycurtis = "Bray-Curtis Distance"
  )
  
  comparison_labels <- c(
    "Sapwood_Heartwood_Within" = "Sapwood-\nHeartwood",
    "Sapwood_Sapwood_Within" = "Sapwood-\nSapwood",
    "Heartwood_Heartwood_Within" = "Heartwood-\nHeartwood",
    "Sapwood_Heartwood_Between" = "Sapwood-\nHeartwood",
    "Sapwood_Sapwood_Between" = "Sapwood-\nSapwood",
    "Heartwood_Heartwood_Between" = "Heartwood-\nHeartwood"
  )
  
  barplot_plot <- ggplot(plot_data, aes(x = Comparison, y = emmean, fill = Comparison)) +
    geom_bar(stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
    geom_text(aes(y = emmean + SE + y_max * 0.05, label = Letters), size = 5) +
    geom_text(aes(y = -y_max * 0.05, label = paste("n =", n)), size = 3) +
    annotate("text", x = 1, y = y_max * 1.2,
             label = paste("Eta² =", round(effect_sizes$Eta2[1], 3))) +
    
    # Brackets
    geom_segment(aes(x = 1, xend = 3, 
                     y = y_max + bracket_height, 
                     yend = y_max + bracket_height), size = 0.5) +
    annotate("text", x = 2, 
             y = y_max + text_height, 
             label = "Within Species", size = 4) +
    geom_segment(aes(x = 4, xend = 6, 
                     y = y_max + bracket_height, 
                     yend = y_max + bracket_height), size = 0.5) +
    annotate("text", x = 5, 
             y = y_max + text_height, 
             label = "Between Species", size = 4) +
    
    scale_x_discrete(labels = comparison_labels) +
    scale_fill_manual(values = c(
      "Sapwood_Heartwood_Within" = "grey50",
      "Sapwood_Sapwood_Within" = "#377EB8",
      "Heartwood_Heartwood_Within" = "#8B4513",
      "Sapwood_Heartwood_Between" = "grey50",
      "Sapwood_Sapwood_Between" = "#377EB8",
      "Heartwood_Heartwood_Between" = "#8B4513"
    )) +
    labs(x = "Comparison", 
         y = metric_labels[metric_name],
         title = paste(toupper(metric_name), "distances between sample types"),
         subtitle = paste("PERMANOVA p =", format.pval(permanova_result$"Pr(>F)"[1], digits = 3))) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save results
  results <- list(
    model = mixed_model,
    emmeans = emm,
    plot = barplot_plot,
    permanova = permanova_result,
    beta_dispersion = beta_test,
    effect_sizes = effect_sizes,
    data = dist_long,
    plot_data = plot_data
  )
  
  # Store in global environment
  assign(paste0("results_", metric_name), results, envir = .GlobalEnv)
  
  # Save plot
  ggsave(paste0("barplot_", metric_name, ".pdf"), plot = barplot_plot, 
         width = 8, height = 6)
  
  # Return results invisibly
  invisible(results)
}

# Process all distance matrices
results_list <- lapply(names(distance_matrices), function(metric) {
  analyze_phyloseq_distances(
    dist_matrix = as.dist(distance_matrices[[metric]]),
    phyloseq_obj = phyloseq_clean,
    metric_name = metric
  )
})

# Name the results
names(results_list) <- names(distance_matrices)