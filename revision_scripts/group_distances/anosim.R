# Additional multivariate analyses
analyze_multivariate_tests <- function(dist_matrices, metadata_clean) {
  require(vegan)
  require(tidyverse)
  
  # Initialize results list
  results_list <- list()
  
  # Function to run tests for each distance matrix
  run_tests <- function(dist_matrix, metric_name) {
    # Ensure samples match between distance matrix and metadata
    dist_samples <- rownames(as.matrix(dist_matrix))
    metadata_samples <- rownames(metadata_clean)
    common_samples <- intersect(dist_samples, metadata_samples)
    
    # Subset distance matrix and metadata
    dist_matrix_subset <- as.dist(as.matrix(dist_matrix)[common_samples, common_samples])
    metadata_subset <- metadata_clean[common_samples, ]
    
    # ANOSIM tests
    anosim_sample_type <- anosim(dist_matrix_subset, 
                                 grouping = metadata_subset$SampleType,
                                 permutations = 999)
    
    anosim_species <- anosim(dist_matrix_subset, 
                             grouping = metadata_subset$Species,
                             permutations = 999)
    
    # MRPP tests
    mrpp_sample_type <- mrpp(dist_matrix_subset, 
                             grouping = metadata_subset$SampleType,
                             permutations = 999)
    
    mrpp_species <- mrpp(dist_matrix_subset, 
                         grouping = metadata_subset$Species,
                         permutations = 999)
    
    # Create summary table
    results_df <- data.frame(
      Metric = metric_name,
      Test = c("ANOSIM", "ANOSIM", "MRPP", "MRPP"),
      Factor = c("Sample Type", "Species", "Sample Type", "Species"),
      Statistic = c(anosim_sample_type$statistic, 
                    anosim_species$statistic,
                    mrpp_sample_type$delta,
                    mrpp_species$delta),
      P_value = c(anosim_sample_type$signif,
                  anosim_species$signif,
                  mrpp_sample_type$Pvalue,
                  mrpp_species$Pvalue),
      Effect_size = c(anosim_sample_type$statistic,  # R statistic for ANOSIM
                      anosim_species$statistic,        # R statistic for ANOSIM
                      mrpp_sample_type$A,              # A statistic for MRPP
                      mrpp_species$A)                  # A statistic for MRPP
    )
    
    return(list(
      anosim_sample_type = anosim_sample_type,
      anosim_species = anosim_species,
      mrpp_sample_type = mrpp_sample_type,
      mrpp_species = mrpp_species,
      summary = results_df
    ))
  }
  
  # Run tests for each distance matrix
  for (metric in names(dist_matrices)) {
    results_list[[metric]] <- run_tests(dist_matrices[[metric]], metric)
  }
  
  # Combine all summaries
  all_summaries <- do.call(rbind, lapply(results_list, function(x) x$summary))
  
  # Save detailed results
  sink("multivariate_test_results.txt")
  cat("=== Multivariate Analysis Results ===\n\n")
  
  for (metric in names(results_list)) {
    cat(paste("\n\n=== Results for", metric, "===\n"))
    
    cat("\nANOSIM Results - Sample Type\n")
    cat("==========================\n")
    print(results_list[[metric]]$anosim_sample_type)
    
    cat("\nANOSIM Results - Species\n")
    cat("=======================\n")
    print(results_list[[metric]]$anosim_species)
    
    cat("\nMRPP Results - Sample Type\n")
    cat("========================\n")
    print(results_list[[metric]]$mrpp_sample_type)
    
    cat("\nMRPP Results - Species\n")
    cat("=====================\n")
    print(results_list[[metric]]$mrpp_species)
  }
  sink()
  
  # Create comparison plot
  plot_data <- all_summaries %>%
    mutate(
      Test_Factor = paste(Test, Factor),
      Significance = ifelse(P_value < 0.05, "*", "")
    )
  
  comparison_plot <- ggplot(plot_data, 
                            aes(x = Test_Factor, y = Effect_size, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = Significance), 
              position = position_dodge(width = 0.9),
              vjust = -0.5) +
    facet_wrap(~Test, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Effect Size", 
         title = "Comparison of Multivariate Tests",
         subtitle = "* indicates p < 0.05")
  
  # Save plot
  ggsave("multivariate_comparison.pdf", comparison_plot, width = 10, height = 6)
  
  return(list(
    results = results_list,
    summary = all_summaries,
    plot = comparison_plot
  ))
}

# Run the analysis with your existing data
multivariate_results <- analyze_multivariate_tests(distance_matrices, metadata_clean)

# Print summary table
print("Summary of all multivariate tests:")
print(multivariate_results$summary)

# Display comparison plot
print(multivariate_results$plot)