# Load necessary libraries
library(tidyverse)
library(boot)
library(car) # For Levene's test

# Function to process and calibrate data
process_gc_data <- function(use_exponential_O2_curve = TRUE) {
  # Load data files
  O2_standards <- read.csv("Downloads/GC_241010 - O2 Standards.csv")
  GHG_standards <- read.csv("Downloads/GC_241010 - GHG Standards.csv")
  GC_data <- read.csv("Downloads/GC_241010 - Run 1 Compiled (6).csv")
  
  # Handle missing values
  GC_data$CH4_Area[331] <- NA
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
  low_range_o2_standards <- c("N2", "O1", "O1.5", "O2", "O2.5")
  high_range_o2_standards <- c("O3", "O4", "O5", "O6")
  low_range_standards <- c("SB1", "SB3a", "SB3b", "SB3", "SB4a")
  high_range_standards <- c("SB4", "SB5a", "SB5", "S3a", "S3b", "S3c")
  
  # Filter data for ranges
  low_range_o2_data <- O2_data %>% dplyr::filter(Sample %in% low_range_o2_standards)
  high_range_o2_data <- O2_data %>% dplyr::filter(Sample %in% high_range_o2_standards)
  all_high_range_o2_data <- dplyr::bind_rows(low_range_o2_data, high_range_o2_data)
  low_range_data <- GHG_data %>% dplyr::filter(Sample %in% low_range_standards)
  
  # Create calibration curves
  if (use_exponential_O2_curve) {
    # Clean O2 data first
    clean_o2_data <- O2_data %>%
      dplyr::filter(!is.na(O2_Area) &
                      !is.na(O2_ppm) & O2_Area > 0 & O2_ppm > 0)
    
    # Calculate initial parameters for exponential fit
    initial_a <- min(clean_o2_data$O2_ppm)
    initial_b <- log(max(clean_o2_data$O2_ppm) / min(clean_o2_data$O2_ppm)) /
      (max(clean_o2_data$O2_Area) - min(clean_o2_data$O2_Area))
    
    # Fit exponential model
    O2_exp_fit <- nls(
      O2_ppm ~ a * exp(b * O2_Area),
      data = clean_o2_data,
      start = list(a = initial_a, b = initial_b),
      control = list(maxiter = 100)
    )
    
    # Extract coefficients
    exp_coeffs <- coef(O2_exp_fit)
    a_val <- exp_coeffs["a"]
    b_val <- exp_coeffs["b"]
    
    # Calculate O2 concentrations using exponential fit
    GC_data$O2_concentration <- a_val * exp(b_val * GC_data$O2_Area)
  } else {
    O2_low_curve <- lm(O2_ppm ~ O2_Area, data = low_range_o2_data)
    O2_high_curve <- lm(O2_ppm ~ O2_Area, data = all_high_range_o2_data)
    max_O2_low_area <- max(low_range_o2_data$O2_Area, na.rm = TRUE)
    GC_data$O2_concentration <- ifelse(
      GC_data$O2_Area <= max_O2_low_area,
      predict(O2_low_curve, newdata = data.frame(O2_Area = GC_data$O2_Area)),
      predict(O2_high_curve, newdata = data.frame(O2_Area = GC_data$O2_Area))
    )
  }
  
  # Process other gases
  N2O_low_curve <- lm(N2O_ppm ~ N2O_Area, data = low_range_data)
  N2O_high_curve <- lm(N2O_ppm ~ N2O_Area, data = GHG_data)
  CH4_low_curve <- lm(CH4_ppm ~ CH4_Area, data = low_range_data)
  CH4_high_curve <- lm(CH4_ppm ~ CH4_Area, data = GHG_data)
  CO2_low_curve <- lm(CO2_ppm ~ CO2_Area, data = low_range_data)
  CO2_high_curve <- lm(CO2_ppm ~ CO2_Area, data = GHG_data)
  
  max_N2O_low_area <- max(low_range_data$N2O_Area, na.rm = TRUE)
  max_CH4_low_area <- max(low_range_data$CH4_Area, na.rm = TRUE)
  max_CO2_low_area <- max(low_range_data$CO2_Area, na.rm = TRUE)
  
  # Calculate concentrations for other gases
  GC_data <- GC_data %>%
    dplyr::mutate(
      N2O_concentration = dplyr::if_else(
        N2O_Area <= max_N2O_low_area,
        predict(N2O_low_curve, newdata = data.frame(N2O_Area = N2O_Area)),
        predict(N2O_high_curve, newdata = data.frame(N2O_Area = N2O_Area))
      ),
      CH4_concentration = dplyr::if_else(
        CH4_Area <= max_CH4_low_area,
        predict(CH4_low_curve, newdata = data.frame(CH4_Area = CH4_Area)),
        predict(CH4_high_curve, newdata = data.frame(CH4_Area = CH4_Area))
      ),
      CO2_concentration = dplyr::if_else(
        CO2_Area <= max_CO2_low_area,
        predict(CO2_low_curve, newdata = data.frame(CO2_Area = CO2_Area)),
        predict(CO2_high_curve, newdata = data.frame(CO2_Area = CO2_Area))
      )
    )
  
  GC_data <- GC_data %>%
    dplyr::mutate(
      CH4_concentration = CH4_concentration * 3,
      O2_concentration = O2_concentration * 3,
      CO2_concentration = CO2_concentration * 3,
      N2O_concentration = N2O_concentration * 3
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
  
  return(list(GC_data = GC_data))  # Ensure GC_data is returned
}

# Function to perform comprehensive statistical analysis
analyze_gas <- function(data, gas_prefix) {
  # Construct the concentration column name
  gas_col <- paste0(gas_prefix, "_concentration")
  
  # Prepare data with more thorough filtering
  df_filtered <- data %>%
    dplyr::select(Species, Tree.No, Tissue, !!sym(gas_col)) %>%
    dplyr::filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(!!sym(gas_col))) %>%
    dplyr::filter(Tissue != "") %>%
    dplyr::mutate(value = as.numeric(!!sym(gas_col))) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::group_by(Species, Tree.No, Tissue) %>%
    dplyr::summarize(value = mean(value), .groups = 'drop')
  
  # Raw t-tests (two-tailed)
  t_test_result <- t.test(value ~ Tissue, data = df_filtered)
  
  # Raw t-tests (one-tailed)
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
  
  # Two-tailed bootstrap p-value
  p_value_two_sided <- 2 * min(mean(bootstrap_results$t >= 0), mean(bootstrap_results$t <= 0))
  
  # One-tailed bootstrap p-values
  p_value_greater <- mean(bootstrap_results$t >= 0)
  p_value_less <- mean(bootstrap_results$t <= 0)
  
  # Variance analysis
  sapwood_data <- df_filtered$value[df_filtered$Tissue == "sapwood"]
  heartwood_data <- df_filtered$value[df_filtered$Tissue == "heartwood"]
  
  # Levene’s test (test for homogeneity of variance)
  levene_test <- car::leveneTest(value ~ Tissue, data = df_filtered)
  
  # Fligner-Killeen test (robust test for homogeneity of variances)
  fligner_test <- fligner.test(value ~ Tissue, data = df_filtered)
  
  # Observed variance difference
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
  
  # Return results
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
      levene_test_p_value = levene_test$`Pr(>F)`[1], # Extract p-value from Levene's test
      fligner_test_p_value = fligner_test$p.value    # Extract p-value from Fligner-Killeen test
    )
  ))
}

# Main analysis
run_analysis <- function(use_exponential_O2_curve = FALSE) {
  # Process data
  processed_data <- process_gc_data(use_exponential_O2_curve)
  GC_data <- processed_data$GC_data  # Extract GC_data
  
  # Analyze each gas using final concentrations
  gas_prefixes <- c("O2", "CH4", "CO2", "N2O")  # Base names for gases
  results <- lapply(gas_prefixes, function(gas)
    analyze_gas(GC_data, gas))
  names(results) <- gas_prefixes
  
  # Print results for each gas using the proper concentration variable
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
    cat("95% CI (Permutation Test):", results[[gas]]$variance$permutation_ci[1], "to", results[[gas]]$variance$permutation_ci[2], "\n")
    
    cat("\nLevene's Test for Equal Variance:\n")
    cat("Levene's Test p-value:", results[[gas]]$variance$levene_test_p_value, "\n")
    
    cat("\nFligner-Killeen Test for Equal Variance:\n")
    cat("Fligner-Killeen Test p-value:", results[[gas]]$variance$fligner_test_p_value, "\n")
  }
  
  # Save processed data to CSV
  write.csv(GC_data, "processed_GC_data.csv", row.names = FALSE)
  
  # Return processed data and results
  return(list(data = GC_data, results = results))
}


# Run the analysis (set use_exponential_O2_curve = TRUE for exponential O2 curve)
results <- run_analysis(use_exponential_O2_curve = FALSE)

GC_data <- results$data
