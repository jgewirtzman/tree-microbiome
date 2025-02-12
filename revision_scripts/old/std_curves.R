# Function to create calibration curves and calculate concentrations with dilution correction
create_calibration_curves <- function(O2_standards, GC_data, 
                                    low_standards, high_standards,
                                    concentration_cutoff = 210000) {
  
  # Merge standards with GC data
  O2_data <- O2_standards %>%
    left_join(GC_data, by = c("Sample" = "Sample.ID")) %>%
    mutate(
      O2_Area = as.numeric(O2_Area),
      O2_ppm = as.numeric(gsub(",", "", X.O2...ppm.))
    )
  
  # Filter standards data based on input selections
  low_range_o2_data <- O2_data %>%
    filter(Sample %in% low_standards)
  
  high_range_o2_data <- O2_data %>%
    filter(Sample %in% high_standards)
  
  # Combine for high range curve
  all_high_range_o2_data <- bind_rows(low_range_o2_data, high_range_o2_data)
  
  # Create calibration curves
  O2_low_curve <- lm(O2_ppm ~ O2_Area, data = low_range_o2_data)
  O2_high_curve <- lm(O2_ppm ~ O2_Area, data = all_high_range_o2_data)
  
  # Get max area for low range
  max_O2_low_area <- max(low_range_o2_data$O2_Area, na.rm = TRUE)
  
  # Calculate concentrations with dilution correction and bounds
  results <- GC_data %>%
    mutate(
      O2_concentration = if_else(
        O2_Area <= max_O2_low_area,
        predict(O2_low_curve, newdata = data.frame(O2_Area = O2_Area)),
        predict(O2_high_curve, newdata = data.frame(O2_Area = O2_Area))
      )
    ) %>%
    mutate(
      O2_concentration = pmax(0, O2_concentration * 3), # Apply dilution factor
      O2_concentration = pmin(concentration_cutoff, O2_concentration) # Apply cutoff
    )
  
  return(list(
    data = results,
    low_curve = O2_low_curve,
    high_curve = O2_high_curve,
    r2_low = summary(O2_low_curve)$r.squared,
    r2_high = summary(O2_high_curve)$r.squared
  ))
}

# Function to analyze wood differences using bootstrap
analyze_wood_differences <- function(data) {
  # Prepare data
  wood_data <- data %>%
    select(Species, Tree.No, Tissue, O2_concentration) %>%
    filter(!is.na(O2_concentration)) %>%
    filter(Tissue %in% c("S", "H")) %>%
    group_by(Species, Tree.No, Tissue) %>%
    summarize(
      mean_O2 = mean(O2_concentration),
      var_O2 = var(O2_concentration),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      names_from = Tissue,
      values_from = c(mean_O2, var_O2),
      names_glue = "{.value}_{Tissue}"
    )
  
  # Bootstrap function for mean difference
  mean_diff <- function(data, indices) {
    sampled_data <- data[indices, ]
    return(mean(sampled_data$mean_O2_S - sampled_data$mean_O2_H))
  }
  
  # Perform bootstrap
  boot_results <- boot(data = wood_data, statistic = mean_diff, R = 10000)
  
  # Calculate confidence intervals
  ci <- boot.ci(boot_results, type = c("perc"))
  
  # One-sided p-value for testing if sapwood > heartwood
  p_value <- mean(boot_results$t <= 0)
  
  return(list(
    boot_results = boot_results,
    ci = ci,
    p_value = p_value,
    summary_stats = wood_data %>%
      summarise(
        mean_diff = mean(mean_O2_S - mean_O2_H),
        sd_diff = sd(mean_O2_S - mean_O2_H)
      )
  ))
}

# Define different calibration scenarios
calibration_scenarios <- list(
  base = list(
    low = c("O1", "O1.5", "O2", "O2.5"),
    high = c("O3", "O4", "O5", "O6")
  ),
  conservative = list(
    low = c("O1", "O1.5", "O2"),
    high = c("O3", "O4", "O5")
  ),
  extended = list(
    low = c("O1", "O1.5", "O2", "O2.5", "O3"),
    high = c("O3", "O4", "O5", "O6")
  )
)

cutoff_values <- c(210000, 180000, 240000) # Different O2 concentration cutoffs

# Run sensitivity analysis
sensitivity_results <- list()

for(scenario_name in names(calibration_scenarios)) {
  for(cutoff in cutoff_values) {
    key <- paste(scenario_name, cutoff, sep = "_")
    
    # Run calibration
    cal_result <- create_calibration_curves(
      O2_standards, 
      GC_data,
      calibration_scenarios[[scenario_name]]$low,
      calibration_scenarios[[scenario_name]]$high,
      cutoff
    )
    
    # Analyze differences
    diff_result <- analyze_wood_differences(cal_result$data)
    
    # Store results
    sensitivity_results[[key]] <- list(
      calibration = cal_result,
      differences = diff_result,
      scenario = scenario_name,
      cutoff = cutoff
    )
  }
}

# Create summary table
summary_df <- map_dfr(sensitivity_results, function(x) {
  tibble(
    scenario = x$scenario,
    cutoff = x$cutoff,
    mean_difference = x$differences$summary_stats$mean_diff,
    ci_lower = x$differences$ci$percent[4],
    ci_upper = x$differences$ci$percent[5],
    p_value = x$differences$p_value,
    r2_low = x$calibration$r2_low,
    r2_high = x$calibration$r2_high
  )
}, .id = "analysis_id")

# Plot sensitivity results
ggplot(summary_df, aes(x = cutoff, y = mean_difference, color = scenario)) +
  geom_point(size = 3) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = scenario), alpha = 0.2) +
  labs(title = "Sensitivity of Sapwood-Heartwood O2 Difference",
       subtitle = "Error bands show bootstrapped 95% confidence intervals",
       x = "O2 Concentration Cutoff (ppm)",
       y = "Mean Difference (Sapwood - Heartwood)") +
  theme_minimal()

# Print summary table
print(summary_df %>%
  arrange(abs(p_value - 0.05)) %>%
  mutate(across(where(is.numeric), ~round(., 3))))