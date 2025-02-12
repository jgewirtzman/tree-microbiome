# Load necessary libraries
library(tidyverse)
library(boot)
library(car)

# Function to process and calibrate data
process_gc_data <- function() {
  # Load data files
  O2_standards <- read.csv("Downloads/GC_241010 - O2 Standards.csv")
  GHG_standards <- read.csv("Downloads/GC_241010 - GHG Standards.csv")
  GC_data <- read.csv("Downloads/GC_241010 - Run 1 Compiled (6).csv")
  
  # Handle missing values
  GC_data$CH4_Area[331] <- NA
  GC_data$O2_Area[which(GC_data$Sample.ID=="N2")] <- NA
  colnames(GC_data) <- make.names(colnames(GC_data))
  
  # Merge standards with GC data
  O2_data <- O2_standards %>%
    dplyr::left_join(GC_data, by = c("Sample" = "Sample.ID"))
  GHG_data <- GHG_standards %>%
    dplyr::left_join(GC_data, by = c("Sample" = "Sample.ID"))
  
  # Convert concentrations to numeric
  O2_data$O2_Area <- as.numeric(O2_data$O2_Area)
  O2_data$O2_ppm <- as.numeric(gsub(",", "", O2_data$X.O2...ppm.))
  GHG_data$N2O_Area <- as.numeric(GHG_data$N2O_Area)
  GHG_data$N2O_ppm <- as.numeric(gsub(",", "", GHG_data$X.N2O...ppm.))
  GHG_data$CH4_Area <- as.numeric(GHG_data$CH4_Area)
  GHG_data$CH4_ppm <- as.numeric(gsub(",", "", GHG_data$X.CH4...ppm.))
  GHG_data$CO2_Area <- as.numeric(GHG_data$CO2_Area)
  GHG_data$CO2_ppm <- as.numeric(gsub(",", "", GHG_data$X.CO2...ppm.))
  
  # Define standard ranges
  # O2 ranges
  low_range_o2_standards <- c("O1", "O1.5", "O2", "O2.5")
  high_range_o2_standards <- c("O3", "O4", "O5", "O6")
  
  # GHG ranges
  low_range_standards <- c("SB1", "SB3a", "SB3b", "SB3", "SB4a")
  high_range_standards <- c("SB4", "SB5a", "SB5", "S3a", "S3b", "S3c")
  
  # Filter data for O2 ranges
  low_range_o2_data <- O2_data %>% dplyr::filter(Sample %in% low_range_o2_standards)
  high_range_o2_data <- O2_data %>% dplyr::filter(Sample %in% high_range_o2_standards)
  all_high_range_o2_data <- O2_data %>% dplyr::filter(Sample %in% c(low_range_o2_standards, high_range_o2_standards))
  
  # Filter data for GHG ranges
  low_range_data <- GHG_data %>% dplyr::filter(Sample %in% low_range_standards)
  high_range_data <- GHG_data %>% dplyr::filter(Sample %in% high_range_standards)
  all_high_range_data <- GHG_data %>% dplyr::filter(Sample %in% c(low_range_standards, high_range_standards))
  
  # Create calibration curves
  O2_low_curve <- lm(O2_ppm ~ O2_Area, data = low_range_o2_data)
  O2_high_curve <- lm(O2_ppm ~ O2_Area, data = all_high_range_o2_data)
  N2O_low_curve <- lm(N2O_ppm ~ N2O_Area, data = low_range_data)
  N2O_high_curve <- lm(N2O_ppm ~ N2O_Area, data = all_high_range_data)
  CH4_low_curve <- lm(CH4_ppm ~ CH4_Area, data = low_range_data)
  CH4_high_curve <- lm(CH4_ppm ~ CH4_Area, data = all_high_range_data)
  CO2_low_curve <- lm(CO2_ppm ~ CO2_Area, data = low_range_data)
  CO2_high_curve <- lm(CO2_ppm ~ CO2_Area, data = all_high_range_data)
  
  # Get maximum peak areas
  max_O2_low_area <- max(low_range_o2_data$O2_Area, na.rm = TRUE)
  max_N2O_low_area <- max(low_range_data$N2O_Area, na.rm = TRUE)
  max_CH4_low_area <- max(low_range_data$CH4_Area, na.rm = TRUE)
  max_CO2_low_area <- max(low_range_data$CO2_Area, na.rm = TRUE)
  
  # Calculate both uncorrected and corrected concentrations
  GC_data <- GC_data %>%
    dplyr::mutate(
      # Uncorrected concentrations (for statistics)
      O2_concentration_uncorrected = dplyr::if_else(
        O2_Area <= max_O2_low_area,
        predict(O2_low_curve, newdata = data.frame(O2_Area = O2_Area)),
        predict(O2_high_curve, newdata = data.frame(O2_Area = O2_Area))
      ),
      N2O_concentration_uncorrected = dplyr::if_else(
        N2O_Area <= max_N2O_low_area,
        predict(N2O_low_curve, newdata = data.frame(N2O_Area = N2O_Area)),
        predict(N2O_high_curve, newdata = data.frame(N2O_Area = N2O_Area))
      ),
      CH4_concentration_uncorrected = dplyr::if_else(
        CH4_Area <= max_CH4_low_area,
        predict(CH4_low_curve, newdata = data.frame(CH4_Area = CH4_Area)),
        predict(CH4_high_curve, newdata = data.frame(CH4_Area = CH4_Area))
      ),
      CO2_concentration_uncorrected = dplyr::if_else(
        CO2_Area <= max_CO2_low_area,
        predict(CO2_low_curve, newdata = data.frame(CO2_Area = CO2_Area)),
        predict(CO2_high_curve, newdata = data.frame(CO2_Area = CO2_Area))
      ),
      
      # Corrected concentrations (for plotting)
      O2_concentration = pmax(0, O2_concentration_uncorrected * 3),
      O2_concentration = pmin(210000, O2_concentration),
      CH4_concentration = pmax(0, CH4_concentration_uncorrected * 3),
      CO2_concentration = pmax(0, CO2_concentration_uncorrected * 3),
      N2O_concentration = pmax(0, N2O_concentration_uncorrected * 3)
    )
  
  # Recode species and tissue types
  species_mapping <- c(
    "BB" = "Betula lenta",
    "H" = "Tsuga canadensis",
    "RM" = "Acer rubrum",
    "RO" = "Quercus rubra",
    "SM" = "Acer saccharum",
    "WP" = "Pinus strobus"
  )
  GC_data$Species <- dplyr::recode(GC_data$Species, !!!species_mapping)
  
  tissue_mapping <- c("S" = "sapwood", "H" = "heartwood")
  GC_data$Tissue <- dplyr::recode(GC_data$Tissue, !!!tissue_mapping)
  
  return(list(
    GC_data = GC_data,
    curves = list(
      O2_low = O2_low_curve,
      O2_high = O2_high_curve,
      N2O_low = N2O_low_curve,
      N2O_high = N2O_high_curve,
      CH4_low = CH4_low_curve,
      CH4_high = CH4_high_curve,
      CO2_low = CO2_low_curve,
      CO2_high = CO2_high_curve
    ),
    standard_data = list(
      O2_low = low_range_o2_data,
      O2_high = high_range_o2_data,
      GHG_low = low_range_data,
      GHG_all = GHG_data
    )
  ))
}

# Function to analyze gas using uncorrected concentrations
analyze_gas <- function(data, gas_prefix) {
  # Use uncorrected concentrations for statistical analysis
  gas_col <- paste0(gas_prefix, "_concentration_uncorrected")
  
  df_filtered <- data %>%
    dplyr::select(Species, Tree.No, Tissue, !!sym(gas_col)) %>%
    dplyr::filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(!!sym(gas_col))) %>%
    dplyr::filter(Tissue != "") %>%
    dplyr::mutate(value = as.numeric(!!sym(gas_col))) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::group_by(Species, Tree.No, Tissue) %>%
    dplyr::summarize(value = mean(value), .groups = 'drop')
  
  # Raw t-tests
  t_test_result <- t.test(value ~ Tissue, data = df_filtered)
  t_test_greater <- t.test(value ~ Tissue, data = df_filtered, alternative = "greater")
  t_test_less <- t.test(value ~ Tissue, data = df_filtered, alternative = "less")
  
  # Bootstrap analysis
  pivot_data <- tidyr::pivot_wider(df_filtered, names_from = Tissue, values_from = value)
  
  mean_diff <- function(data, indices) {
    sampled_data <- data[indices, ]
    sapwood_vals <- sampled_data$sapwood
    heartwood_vals <- sampled_data$heartwood
    return(mean(sapwood_vals - heartwood_vals))
  }
  
  set.seed(123)
  bootstrap_results <- boot(data = pivot_data, statistic = mean_diff, R = 10000)
  ci_results <- boot.ci(bootstrap_results, type = c("norm", "basic", "perc"))
  
  p_value_two_sided <- 2 * min(mean(bootstrap_results$t >= 0), mean(bootstrap_results$t <= 0))
  p_value_greater <- mean(bootstrap_results$t >= 0)
  p_value_less <- mean(bootstrap_results$t <= 0)
  
  # Variance analysis
  sapwood_data <- df_filtered$value[df_filtered$Tissue == "sapwood"]
  heartwood_data <- df_filtered$value[df_filtered$Tissue == "heartwood"]
  
  levene_test <- car::leveneTest(value ~ Tissue, data = df_filtered)
  fligner_test <- fligner.test(value ~ Tissue, data = df_filtered)
  
  observed_var_diff <- var(sapwood_data) - var(heartwood_data)
  
  # Permutation test for variance
  n_perm <- 10000
  perm_var_diffs <- numeric(n_perm)
  combined_data <- c(sapwood_data, heartwood_data)
  group_labels <- c(rep("sapwood", length(sapwood_data)), rep("heartwood", length(heartwood_data)))
  
  set.seed(123)
  for (i in 1:n_perm) {
    perm_labels <- sample(group_labels)
    perm_var_diffs[i] <- var(combined_data[perm_labels == "sapwood"]) - var(combined_data[perm_labels == "heartwood"])
  }
  
  var_p_value <- mean(abs(perm_var_diffs) >= abs(observed_var_diff))
  var_ci <- quantile(perm_var_diffs, c(0.025, 0.975))
  
  return(list(
    t_test = list(
      two_tailed = t_test_result,
      greater = t_test_greater,
      less = t_test_less
    ),
    bootstrap = list(
      estimate = mean(bootstrap_results$t),
      ci = ci_results,
      p_value_two_sided = p_value_two_sided,
      p_value_greater = p_value_greater,
      p_value_less = p_value_less
    ),
    variance = list(
      observed_diff = observed_var_diff,
      permutation_p_value = var_p_value,
      permutation_ci = var_ci,
      levene_test_p_value = levene_test$`Pr(>F)`[1],
      fligner_test_p_value = fligner_test$p.value
    )
  ))
}

# Function to plot standard curves
plot_standard_curve <- function(low_data, high_data, low_curve, high_curve, analyte, x, y, range_type) {
  low_summary <- summary(low_curve)
  high_summary <- summary(high_curve)
  
  low_equation <- paste("Low:", round(coef(low_curve)[2], 3), "*x +", round(coef(low_curve)[1], 3))
  high_equation <- paste("High:", round(coef(high_curve)[2], 3), "*x +", round(coef(high_curve)[1], 3))
  
  low_r2 <- paste("R² (low) =", round(low_summary$r.squared, 3))
  high_r2 <- paste("R² (high) =", round(high_summary$r.squared, 3))
  
  combined_high_data <- rbind(low_data, high_data)
  
  low_label <- paste("Low Range Standards\n", low_equation, "\n", low_r2)
  high_label <- paste("All Standards\n", high_equation, "\n", high_r2)
  
  ggplot() +
    geom_point(data = combined_high_data, aes_string(x = x, y = y, color = "'All Standards'"), alpha = 0.6) +
    geom_smooth(data = combined_high_data, aes_string(x = x, y = y, color = "'All Standards'"), 
                method = "lm", se = FALSE, linetype = "solid") +
    # Plot low range
    geom_point(data = low_data, aes_string(x = x, y = y, color = "'Low Range'")) +
    geom_smooth(data = low_data, aes_string(x = x, y = y, color = "'Low Range'"), 
                method = "lm", se = FALSE, linetype = "solid") +
    # Custom colors
    scale_color_manual(name = "Standard Curves",
                       values = c("Low Range" = "blue", "All Standards" = "red"),
                       labels = c("Low Range" = low_label, "All Standards" = high_label)) +
    # Labels and theme
    labs(title = paste(analyte, "Standard Curves"),
         subtitle = paste("Showing low range standards (blue) and all standards combined (red)"),
         x = "Peak Area",
         y = "Concentration (ppm)") +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10)
    )
}

# Function to plot low-range curves specifically
plot_low_range_curve <- function(low_data, low_curve, analyte, x, y) {
  low_summary <- summary(low_curve)
  low_equation <- paste("Low:", round(coef(low_curve)[2], 3), "*x +", round(coef(low_curve)[1], 3))
  low_r2 <- paste("R² (low) =", round(low_summary$r.squared, 3))
  
  x_min <- min(low_data[[x]], na.rm = TRUE)
  x_max <- max(low_data[[x]], na.rm = TRUE)
  
  ggplot() +
    geom_point(data = low_data, aes_string(x = x, y = y), color = "blue") +
    geom_smooth(data = low_data, aes_string(x = x, y = y), method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, label = paste(low_equation, low_r2, sep = "\n"), 
             hjust = 1.2, vjust = 2, size = 3.5, color = "blue", parse = FALSE) +
    scale_x_continuous(limits = c(x_min, x_max)) +
    labs(title = paste(analyte, "Standard Curve (Low-Range)"),
         x = "Peak Area",
         y = "Concentration (ppm)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Function to plot sample distributions
plot_sample_distributions <- function(GC_data) {
  # Add logging to check data structure
  cat("Columns in GC_data:", paste(names(GC_data), collapse=", "), "\n")
  
  # Filter for samples only and gather gas data
  sample_data <- GC_data %>%
    filter(!is.na(Species) & !is.na(Tissue)) %>%  # Filter for actual samples
    dplyr::select(Sample.ID, O2_concentration, CH4_concentration, 
           CO2_concentration, N2O_concentration) %>%
    gather(key = "Gas", value = "Concentration", -Sample.ID) %>%
    mutate(Gas = gsub("_concentration", "", Gas))
  
  # Add logging to check filtered data
  cat("Number of samples after filtering:", nrow(sample_data), "\n")
  cat("Unique gases:", paste(unique(sample_data$Gas), collapse=", "), "\n")
  
  # Create distribution plots
  ggplot(sample_data, aes(x = Concentration)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.2, fill = "blue") +
    facet_wrap(~Gas, scales = "free", ncol = 2) +
    theme_minimal() +
    labs(title = "Distribution of Gas Concentrations in Samples",
         x = "Concentration (ppm)",
         y = "Density") +
    theme(strip.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))
}

# Function to plot standards validation
plot_standards_validation <- function(O2_data, GHG_data) {
  # Process O2 standards
  o2_validation <- O2_data %>%
    dplyr::select(Sample, O2_Area, O2_ppm, X.O2...ppm.) %>%
    mutate(
      Expected = as.numeric(gsub(",", "", X.O2...ppm.)),
      Measured = O2_ppm,
      Gas = "O2"
    ) %>%
    dplyr::select(Gas, Expected, Measured)
  
  # Process other gases
  n2o_validation <- GHG_data %>%
    dplyr::select(Sample, N2O_Area, N2O_ppm, X.N2O...ppm.) %>%
    mutate(
      Expected = as.numeric(gsub(",", "", X.N2O...ppm.)),
      Measured = N2O_ppm,
      Gas = "N2O"
    ) %>%
    dplyr::select(Gas, Expected, Measured)
  
  ch4_validation <- GHG_data %>%
    dplyr::select(Sample, CH4_Area, CH4_ppm, X.CH4...ppm.) %>%
    mutate(
      Expected = as.numeric(gsub(",", "", X.CH4...ppm.)),
      Measured = CH4_ppm,
      Gas = "CH4"
    ) %>%
    dplyr::select(Gas, Expected, Measured)
  
  co2_validation <- GHG_data %>%
    dplyr::select(Sample, CO2_Area, CO2_ppm, X.CO2...ppm.) %>%
    mutate(
      Expected = as.numeric(gsub(",", "", X.CO2...ppm.)),
      Measured = CO2_ppm,
      Gas = "CO2"
    ) %>%
    dplyr::select(Gas, Expected, Measured)
  
  # Combine all gases
  all_validation <- bind_rows(
    o2_validation,
    n2o_validation,
    ch4_validation,
    co2_validation
  ) %>%
    filter(!is.na(Expected) & !is.na(Measured))
  
  # Create validation plot
  ggplot(all_validation, aes(x = Expected, y = Measured)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    facet_wrap(~Gas, scales = "free", ncol = 2) +
    theme_minimal() +
    labs(title = "Standard Values: Expected vs Measured",
         x = "Expected Concentration (ppm)",
         y = "Measured Concentration (ppm)") +
    theme(strip.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))
}

# Main analysis function
run_analysis <- function() {
  # Process data
  processed_data <- process_gc_data()
  GC_data <- processed_data$GC_data
  
  # Create standard curve plots
  plots <- list()
  
  # O2 curves
  plots$O2_low <- plot_standard_curve(
    processed_data$standard_data$O2_low,
    processed_data$standard_data$O2_high,
    processed_data$curves$O2_low,
    processed_data$curves$O2_high,
    "O2", "O2_Area", "O2_ppm", "Low-Range"
  )
  
  # GHG curves
  for (gas in c("N2O", "CH4", "CO2")) {
    area_col <- paste0(gas, "_Area")
    ppm_col <- paste0(gas, "_ppm")
    
    plots[[paste0(gas, "_low")]] <- plot_standard_curve(
      processed_data$standard_data$GHG_low,
      processed_data$standard_data$GHG_all,
      processed_data$curves[[paste0(gas, "_low")]],
      processed_data$curves[[paste0(gas, "_high")]],
      gas, area_col, ppm_col, "Low-Range"
    )
  }
  
  # Create diagnostic plots
  plots$sample_distributions <- plot_sample_distributions(GC_data)
  plots$standards_validation <- plot_standards_validation(
    processed_data$standard_data$O2_low,
    processed_data$standard_data$GHG_all
  )
  
  # Analyze each gas
  gas_prefixes <- c("O2", "CH4", "CO2", "N2O")
  results <- lapply(gas_prefixes, function(gas)
    analyze_gas(GC_data, gas))
  names(results) <- gas_prefixes
  
  # Print results for each gas
  for (gas in gas_prefixes) {
    cat("\n=== Results for", gas, "Concentration ===\n")
    
    cat("\nT-Tests:\n")
    cat("Two-tailed p-value:", results[[gas]]$t_test$two_tailed$p.value, "\n")
    cat("One-tailed (greater) p-value:", results[[gas]]$t_test$greater$p.value, "\n")
    cat("One-tailed (less) p-value:", results[[gas]]$t_test$less$p.value, "\n")
    
    cat("\nBootstrap Analysis:\n")
    cat("Mean Difference:", results[[gas]]$bootstrap$estimate, "\n")
    cat("Two-sided p-value:", results[[gas]]$bootstrap$p_value_two_sided, "\n")
    cat("One-sided p-value (greater):", results[[gas]]$bootstrap$p_value_greater, "\n")
    cat("One-sided p-value (less):", results[[gas]]$bootstrap$p_value_less, "\n")
    
    cat("\nVariance Analysis:\n")
    cat("Observed Variance Difference:", results[[gas]]$variance$observed_diff, "\n")
    cat("Permutation Test p-value:", results[[gas]]$variance$permutation_p_value, "\n")
    cat("95% CI (Permutation Test):", results[[gas]]$variance$permutation_ci[1], "to", 
        results[[gas]]$variance$permutation_ci[2], "\n")
    cat("\nLevene's Test p-value:", results[[gas]]$variance$levene_test_p_value, "\n")
    cat("Fligner-Killeen Test p-value:", results[[gas]]$variance$fligner_test_p_value, "\n")
  }
  
  # Save processed data to CSV
  write.csv(GC_data, "processed_GC_data.csv", row.names = FALSE)
  
  # Return all results
  return(list(
    data = GC_data,
    results = results,
    plots = plots
  ))
}

# Run the analysis
results <- run_analysis()

# Display plots
print(results$plots$sample_distributions)
print(results$plots$standards_validation)
print(results$plots$O2_low)
print(results$plots$CH4_low)
print(results$plots$CO2_low)
print(results$plots$N2O_low)

# Save the processed data object for use in other scripts
GC_data <- results$data