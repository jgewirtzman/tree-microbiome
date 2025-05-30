# Comprehensive Statistical Diagnostics for All Three Datasets
# Gas concentrations, microbial functions, and ddPCR data

library(tidyverse)
library(moments)
library(nortest)
library(car)

# Function to run diagnostics for any dataset
run_comprehensive_diagnostics <- function(data, value_column, group_column, dataset_name) {
  cat("\n", rep("=", 70), "\n")
  cat("DIAGNOSTICS FOR", toupper(dataset_name), "\n")
  cat("Variable:", value_column, "\n")
  cat(rep("=", 70), "\n")
  
  # Prepare data
  df_clean <- data %>%
    select(all_of(c(value_column, group_column))) %>%
    filter(!is.na(.data[[value_column]]) & !is.na(.data[[group_column]])) %>%
    mutate(
      value = as.numeric(.data[[value_column]]),
      group = as.factor(.data[[group_column]])
    ) %>%
    filter(!is.na(value))
  
  # Split by group
  groups <- unique(df_clean$group)
  if (length(groups) != 2) {
    cat("ERROR: Expected exactly 2 groups, found:", length(groups), "\n")
    cat("Groups found:", paste(groups, collapse = ", "), "\n")
    return(NULL)
  }
  
  group1_data <- df_clean$value[df_clean$group == groups[1]]
  group2_data <- df_clean$value[df_clean$group == groups[2]]
  
  cat("Sample sizes:\n")
  cat(" ", groups[1], ": n =", length(group1_data), "\n")
  cat(" ", groups[2], ": n =", length(group2_data), "\n\n")
  
  # 1. NORMALITY TESTS
  cat("1. NORMALITY TESTS\n")
  cat(rep("-", 40), "\n")
  
  # Shapiro-Wilk test
  group1_shapiro <- shapiro.test(group1_data)
  group2_shapiro <- shapiro.test(group2_data)
  
  cat("Shapiro-Wilk test for normality:\n")
  cat(" ", groups[1], ": W =", round(group1_shapiro$statistic, 4), 
      ", p =", format(group1_shapiro$p.value, scientific = TRUE, digits = 3), "\n")
  cat(" ", groups[2], ": W =", round(group2_shapiro$statistic, 4), 
      ", p =", format(group2_shapiro$p.value, scientific = TRUE, digits = 3), "\n")
  
  # Anderson-Darling test (if samples large enough)
  if (length(group1_data) >= 7 && length(group2_data) >= 7) {
    group1_ad <- ad.test(group1_data)
    group2_ad <- ad.test(group2_data)
    
    cat("Anderson-Darling test for normality:\n")
    cat(" ", groups[1], ": A =", round(group1_ad$statistic, 4), 
        ", p =", format(group1_ad$p.value, scientific = TRUE, digits = 3), "\n")
    cat(" ", groups[2], ": A =", round(group2_ad$statistic, 4), 
        ", p =", format(group2_ad$p.value, scientific = TRUE, digits = 3), "\n")
  }
  
  # 2. EQUALITY OF VARIANCES
  cat("\n2. EQUALITY OF VARIANCES\n")
  cat(rep("-", 40), "\n")
  
  # F-test (parametric, assumes normality)
  f_test <- var.test(group1_data, group2_data)
  cat("F-test for equal variances:\n")
  cat("  F =", round(f_test$statistic, 4), ", p =", 
      format(f_test$p.value, scientific = TRUE, digits = 3), "\n")
  
  # Levene's test (less sensitive to non-normality)
  levene_result <- car::leveneTest(value ~ group, data = df_clean)
  cat("Levene's test for equal variances:\n")
  cat("  F =", round(levene_result$`F value`[1], 4), ", p =", 
      format(levene_result$`Pr(>F)`[1], scientific = TRUE, digits = 3), "\n")
  
  # Fligner-Killeen test (non-parametric)
  fligner_result <- fligner.test(value ~ group, data = df_clean)
  cat("Fligner-Killeen test for equal variances:\n")
  cat("  χ² =", round(fligner_result$statistic, 4), ", p =", 
      format(fligner_result$p.value, scientific = TRUE, digits = 3), "\n")
  
  # 3. DESCRIPTIVE STATISTICS
  cat("\n3. DESCRIPTIVE STATISTICS\n")
  cat(rep("-", 40), "\n")
  
  for (i in 1:2) {
    group_data <- if (i == 1) group1_data else group2_data
    cat(groups[i], ":\n")
    cat("  Mean:", round(mean(group_data), 6), "\n")
    cat("  Median:", round(median(group_data), 6), "\n")
    cat("  SD:", round(sd(group_data), 6), "\n")
    cat("  IQR:", round(IQR(group_data), 6), "\n")
    cat("  Skewness:", round(moments::skewness(group_data), 4), "\n")
    cat("  Kurtosis:", round(moments::kurtosis(group_data), 4), "\n")
    cat("  Min:", round(min(group_data), 6), "\n")
    cat("  Max:", round(max(group_data), 6), "\n\n")
  }
  
  # 4. OUTLIER DETECTION
  cat("4. OUTLIER DETECTION\n")
  cat(rep("-", 40), "\n")
  
  for (i in 1:2) {
    group_data <- if (i == 1) group1_data else group2_data
    q1 <- quantile(group_data, 0.25)
    q3 <- quantile(group_data, 0.75)
    iqr <- q3 - q1
    outliers <- sum(group_data < (q1 - 1.5 * iqr) | group_data > (q3 + 1.5 * iqr))
    cat(groups[i], ":", outliers, "outliers (IQR method)\n")
  }
  
  # 5. ZERO/NEGATIVE VALUE CHECK
  cat("\n5. DATA CHARACTERISTICS\n")
  cat(rep("-", 40), "\n")
  
  zero_count <- sum(df_clean$value == 0)
  negative_count <- sum(df_clean$value < 0)
  
  cat("Zero values:", zero_count, "\n")
  cat("Negative values:", negative_count, "\n")
  cat("Range:", round(min(df_clean$value), 6), "to", round(max(df_clean$value), 6), "\n")
  
  # Check if data might benefit from transformation
  if (min(df_clean$value) > 0) {
    cat("Log transformation possible (all values > 0)\n")
  } else if (min(df_clean$value) >= 0) {
    cat("Log+1 transformation possible (all values ≥ 0)\n")
  }
  
  # 6. RECOMMENDATIONS
  cat("\n6. STATISTICAL TEST RECOMMENDATIONS\n")
  cat(rep("-", 50), "\n")
  
  # Decision logic
  normal_assumption <- group1_shapiro$p.value > 0.05 & group2_shapiro$p.value > 0.05
  equal_var_assumption <- levene_result$`Pr(>F)`[1] > 0.05
  large_sample <- length(group1_data) >= 30 & length(group2_data) >= 30
  very_small_sample <- length(group1_data) < 10 | length(group2_data) < 10
  
  cat("Assumption checks:\n")
  cat("  Both groups normal (Shapiro-Wilk p > 0.05):", normal_assumption, "\n")
  cat("  Equal variances (Levene p > 0.05):", equal_var_assumption, "\n")
  cat("  Large samples (n ≥ 30 each):", large_sample, "\n")
  cat("  Very small samples (n < 10 either group):", very_small_sample, "\n\n")
  
  # Recommendation logic
  if (very_small_sample) {
    cat("RECOMMENDATION: NON-PARAMETRIC test (small sample size)\n")
    cat("  - Very small samples may not meet CLT assumptions\n")
    cat("  - Use Mann-Whitney U test\n")
    cat("  - Consider exact permutation tests\n")
  } else if (normal_assumption & equal_var_assumption) {
    cat("RECOMMENDATION: PARAMETRIC t-test\n")
    cat("  - Both normality and equal variance assumptions met\n")
    cat("  - Two-sample t-test with equal variances\n")
  } else if (normal_assumption & !equal_var_assumption) {
    cat("RECOMMENDATION: PARAMETRIC t-test with Welch correction\n")
    cat("  - Normality assumption met, but unequal variances\n")
    cat("  - Welch's t-test (unequal variances)\n")
  } else if (large_sample) {
    cat("RECOMMENDATION: Either PARAMETRIC (CLT) or NON-PARAMETRIC\n")
    cat("  - Large samples invoke Central Limit Theorem\n")
    cat("  - t-test should be robust, but Mann-Whitney U also appropriate\n")
    cat("  - Consider non-parametric if highly skewed\n")
  } else {
    cat("RECOMMENDATION: NON-PARAMETRIC test\n")
    cat("  - Normality assumptions violated with moderate samples\n")
    cat("  - Mann-Whitney U test recommended\n")
    cat("  - Bootstrap resampling also appropriate\n")
  }
  
  # Additional recommendations
  cat("\nAdditional considerations:\n")
  if (zero_count > 0.1 * nrow(df_clean)) {
    cat("  - High proportion of zeros - consider zero-inflated models\n")
  }
  if (min(df_clean$value) > 0 && (group1_shapiro$p.value < 0.05 | group2_shapiro$p.value < 0.05)) {
    cat("  - Consider log transformation if highly skewed\n")
  }
  
  cat("\nFor variance comparison:\n")
  if (normal_assumption) {
    cat("  - Use F-test for variance comparison\n")
  } else {
    cat("  - Use Fligner-Killeen test (non-parametric variance test)\n")
  }
  
  # Return summary
  return(list(
    dataset = dataset_name,
    variable = value_column,
    n1 = length(group1_data),
    n2 = length(group2_data),
    normal_group1 = group1_shapiro$p.value,
    normal_group2 = group2_shapiro$p.value,
    equal_variances = levene_result$`Pr(>F)`[1],
    recommendation = case_when(
      very_small_sample ~ "nonparametric_small",
      normal_assumption & equal_var_assumption ~ "parametric",
      normal_assumption & !equal_var_assumption ~ "parametric_welch",
      large_sample ~ "either_clt",
      TRUE ~ "nonparametric"
    ),
    groups = groups
  ))
}

# Load datasets
cat("COMPREHENSIVE STATISTICAL DIAGNOSTICS\n")
cat("=====================================\n\n")

# 1. Load ddPCR data
cat("Loading ddPCR data...\n")
ddpcr <- read.csv("/Users/jongewirtzman/Downloads/ddPCR_tree_transposed_data.csv")
ddpcr_filtered <- ddpcr %>%
  mutate(core_type = dplyr::recode(core_type, 
                                   "Inner" = "heartwood", 
                                   "Outer" = "sapwood"))

# 2. Load FAPROTAX data  
cat("Loading FAPROTAX data...\n")
faprotax <- read.csv("/Users/jongewirtzman/Downloads/metabolisms/16S_metabolisms_weighted.csv")
names(faprotax)[1] <- "SampleID"
faprotax_filtered <- faprotax %>%
  mutate(
    Sample_Type = case_when(
      grepl("Inner", SampleID) ~ "Inner",
      grepl("Outer", SampleID) ~ "Outer",
      grepl("Organic", SampleID) ~ "Organic",
      grepl("Mineral", SampleID) ~ "Mineral",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!(Sample_Type %in% c("Organic", "Mineral"))) %>%
  mutate(core_type = dplyr::recode(Sample_Type, 
                                   "Inner" = "heartwood", 
                                   "Outer" = "sapwood"))

# Store all diagnostic results
all_diagnostics <- list()

# 3. Run diagnostics for ddPCR data
cat("\n", rep("#", 80), "\n")
cat("ANALYZING ddPCR DATASETS\n")
cat(rep("#", 80), "\n")

all_diagnostics$pmoa <- run_comprehensive_diagnostics(
  ddpcr_filtered, "pmoa_loose", "core_type", "pmoA Gene Abundance"
)

all_diagnostics$mcra <- run_comprehensive_diagnostics(
  ddpcr_filtered, "mcra_probe_loose", "core_type", "mcrA Gene Abundance"
)

# 4. Run diagnostics for FAPROTAX data
cat("\n", rep("#", 80), "\n")
cat("ANALYZING MICROBIAL FUNCTIONAL DATASETS\n")
cat(rep("#", 80), "\n")

all_diagnostics$methanotrophy <- run_comprehensive_diagnostics(
  faprotax_filtered, "methanotrophy", "core_type", "Methanotrophy Function"
)

all_diagnostics$methanogenesis <- run_comprehensive_diagnostics(
  faprotax_filtered, "hydrogenotrophic_methanogenesis", "core_type", "Hydrogenotrophic Methanogenesis"
)

all_diagnostics$chemoheterotrophy <- run_comprehensive_diagnostics(
  faprotax_filtered, "aerobic_chemoheterotrophy", "core_type", "Aerobic Chemoheterotrophy"
)

all_diagnostics$fermentation <- run_comprehensive_diagnostics(
  faprotax_filtered, "fermentation", "core_type", "Fermentation"
)

# 5. Summary table
cat("\n", rep("#", 80), "\n")
cat("SUMMARY OF ALL RECOMMENDATIONS\n")
cat(rep("#", 80), "\n")

summary_df <- data.frame(
  Dataset = character(),
  Variable = character(),
  N1 = integer(),
  N2 = integer(),
  Normal_Group1 = numeric(),
  Normal_Group2 = numeric(),
  Equal_Var_P = numeric(),
  Recommendation = character(),
  stringsAsFactors = FALSE
)

for (name in names(all_diagnostics)) {
  if (!is.null(all_diagnostics[[name]])) {
    summary_df <- rbind(summary_df, data.frame(
      Dataset = all_diagnostics[[name]]$dataset,
      Variable = all_diagnostics[[name]]$variable,
      N1 = all_diagnostics[[name]]$n1,
      N2 = all_diagnostics[[name]]$n2,
      Normal_Group1 = round(all_diagnostics[[name]]$normal_group1, 4),
      Normal_Group2 = round(all_diagnostics[[name]]$normal_group2, 4),
      Equal_Var_P = round(all_diagnostics[[name]]$equal_variances, 4),
      Recommendation = all_diagnostics[[name]]$recommendation,
      stringsAsFactors = FALSE
    ))
  }
}

print(summary_df)

cat("\n\nKEY:\n")
cat("parametric = Standard t-test (equal variances)\n")
cat("parametric_welch = Welch's t-test (unequal variances)\n") 
cat("nonparametric = Mann-Whitney U test\n")
cat("nonparametric_small = Mann-Whitney U (small samples)\n")
cat("either_clt = Either parametric or non-parametric (large samples)\n")

# Save results
write.csv(summary_df, "statistical_diagnostics_summary.csv", row.names = FALSE)







###


# Non-parametric Statistical Tests for All Datasets
# Based on diagnostic results showing non-normal distributions

library(tidyverse)
library(broom)

# Function to run Mann-Whitney U test and format results
run_mann_whitney <- function(data, value_col, group_col, dataset_name, variable_name) {
  cat("\n", rep("=", 60), "\n")
  cat("Mann-Whitney U Test:", toupper(dataset_name), "-", variable_name, "\n")
  cat(rep("=", 60), "\n")
  
  # Prepare data
  df_clean <- data %>%
    select(all_of(c(value_col, group_col))) %>%
    filter(!is.na(.data[[value_col]]) & !is.na(.data[[group_col]])) %>%
    mutate(
      value = as.numeric(.data[[value_col]]),
      group = as.factor(.data[[group_col]])
    ) %>%
    filter(!is.na(value))
  
  # Get group names
  groups <- sort(unique(df_clean$group))
  group1_data <- df_clean$value[df_clean$group == groups[1]]
  group2_data <- df_clean$value[df_clean$group == groups[2]]
  
  # Sample sizes
  n1 <- length(group1_data)
  n2 <- length(group2_data)
  
  # Descriptive statistics
  cat("Sample sizes:\n")
  cat("  ", groups[1], ": n =", n1, "\n")
  cat("  ", groups[2], ": n =", n2, "\n\n")
  
  cat("Descriptive statistics:\n")
  cat("  ", groups[1], ": Median =", round(median(group1_data), 4), 
      ", IQR =", round(IQR(group1_data), 4), "\n")
  cat("  ", groups[2], ": Median =", round(median(group2_data), 4), 
      ", IQR =", round(IQR(group2_data), 4), "\n\n")
  
  # Mann-Whitney U test (two-sided)
  mw_test <- wilcox.test(value ~ group, data = df_clean, exact = FALSE)
  
  # One-sided tests
  mw_greater <- wilcox.test(group1_data, group2_data, alternative = "greater", exact = FALSE)
  mw_less <- wilcox.test(group1_data, group2_data, alternative = "less", exact = FALSE)
  
  cat("Mann-Whitney U Test Results:\n")
  cat("  Two-sided test: W =", mw_test$statistic, ", p =", format(mw_test$p.value, scientific = TRUE, digits = 3), "\n")
  cat("  One-sided (", groups[1], " > ", groups[2], "): p =", format(mw_greater$p.value, scientific = TRUE, digits = 3), "\n")
  cat("  One-sided (", groups[1], " < ", groups[2], "): p =", format(mw_less$p.value, scientific = TRUE, digits = 3), "\n")
  
  # Effect size (rank-biserial correlation)
  r_effect <- 1 - (2 * mw_test$statistic) / (n1 * n2)
  cat("  Effect size (rank-biserial r) =", round(r_effect, 4), "\n")
  
  # Fligner-Killeen test for variance differences
  fk_test <- fligner.test(value ~ group, data = df_clean)
  cat("\nFligner-Killeen Test for Equal Variances:\n")
  cat("  χ² =", round(fk_test$statistic, 4), ", p =", format(fk_test$p.value, scientific = TRUE, digits = 3), "\n")
  
  # Return results for summary
  return(list(
    dataset = dataset_name,
    variable = variable_name,
    n1 = n1,
    n2 = n2,
    group1 = groups[1],
    group2 = groups[2],
    median1 = median(group1_data),
    median2 = median(group2_data),
    iqr1 = IQR(group1_data),
    iqr2 = IQR(group2_data),
    w_statistic = mw_test$statistic,
    p_two_sided = mw_test$p.value,
    p_greater = mw_greater$p.value,
    p_less = mw_less$p.value,
    effect_size = r_effect,
    variance_test_stat = fk_test$statistic,
    variance_p_value = fk_test$p.value
  ))
}

# =============================================================================
# 1. GAS CONCENTRATION DATA
# =============================================================================
cat("ANALYZING GAS CONCENTRATION DATA (NON-PARAMETRIC)\n")
cat(rep("#", 80), "\n")

# Prepare gas data (same filtering as original analysis)
gas_results <- list()

for (gas in c("O2", "CH4", "CO2", "N2O")) {
  gas_col <- paste0(gas, "_concentration_uncorrected")
  
  # Filter and aggregate data (same as original analysis)
  df_filtered <- GC_data %>%
    select(Species, Tree.No, Tissue, all_of(gas_col)) %>%
    filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(.data[[gas_col]])) %>%
    filter(Tissue != "") %>%
    mutate(value = as.numeric(.data[[gas_col]])) %>%
    filter(!is.na(value)) %>%
    group_by(Species, Tree.No, Tissue) %>%
    summarize(value = mean(value), .groups = 'drop')
  
  gas_results[[gas]] <- run_mann_whitney(
    df_filtered, "value", "Tissue", "Gas Concentrations", paste(gas, "Concentration")
  )
}

# =============================================================================
# 2. ddPCR GENE ABUNDANCE DATA
# =============================================================================
cat("\n\nANALYZING ddPCR GENE ABUNDANCE DATA (NON-PARAMETRIC)\n")
cat(rep("#", 80), "\n")

# Load ddPCR data
ddpcr <- read.csv("/Users/jongewirtzman/Downloads/ddPCR_tree_transposed_data.csv")
ddpcr_filtered <- ddpcr %>%
  mutate(core_type = dplyr::recode(core_type, 
                                   "Inner" = "heartwood", 
                                   "Outer" = "sapwood"))

ddpcr_results <- list()

# pmoA gene
ddpcr_results$pmoa <- run_mann_whitney(
  ddpcr_filtered, "pmoa_loose", "core_type", "ddPCR", "pmoA Gene Abundance"
)

# mcrA gene
ddpcr_results$mcra <- run_mann_whitney(
  ddpcr_filtered, "mcra_probe_loose", "core_type", "ddPCR", "mcrA Gene Abundance"
)

# =============================================================================
# 3. MICROBIAL FUNCTIONAL DATA (FAPROTAX)
# =============================================================================
cat("\n\nANALYZING MICROBIAL FUNCTIONAL DATA (NON-PARAMETRIC)\n")
cat(rep("#", 80), "\n")

# Load FAPROTAX data
faprotax <- read.csv("/Users/jongewirtzman/Downloads/metabolisms/16S_metabolisms_weighted.csv")
names(faprotax)[1] <- "SampleID"
faprotax_filtered <- faprotax %>%
  mutate(
    Sample_Type = case_when(
      grepl("Inner", SampleID) ~ "Inner",
      grepl("Outer", SampleID) ~ "Outer",
      grepl("Organic", SampleID) ~ "Organic",
      grepl("Mineral", SampleID) ~ "Mineral",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!(Sample_Type %in% c("Organic", "Mineral"))) %>%
  mutate(core_type = dplyr::recode(Sample_Type, 
                                   "Inner" = "heartwood", 
                                   "Outer" = "sapwood"))

functional_results <- list()

# Methanotrophy
functional_results$methanotrophy <- run_mann_whitney(
  faprotax_filtered, "methanotrophy", "core_type", "Microbial Functions", "Methanotrophy"
)

# Hydrogenotrophic methanogenesis
functional_results$methanogenesis <- run_mann_whitney(
  faprotax_filtered, "hydrogenotrophic_methanogenesis", "core_type", "Microbial Functions", "Hydrogenotrophic Methanogenesis"
)

# Aerobic chemoheterotrophy
functional_results$chemoheterotrophy <- run_mann_whitney(
  faprotax_filtered, "aerobic_chemoheterotrophy", "core_type", "Microbial Functions", "Aerobic Chemoheterotrophy"
)

# Fermentation
functional_results$fermentation <- run_mann_whitney(
  faprotax_filtered, "fermentation", "core_type", "Microbial Functions", "Fermentation"
)

# =============================================================================
# 4. COMPREHENSIVE SUMMARY TABLE
# =============================================================================
cat("\n\n", rep("#", 80), "\n")
cat("COMPREHENSIVE SUMMARY OF ALL NON-PARAMETRIC TESTS\n")
cat(rep("#", 80), "\n")

# Combine all results
all_results <- c(gas_results, ddpcr_results, functional_results)

# Create summary dataframe
summary_df <- data.frame(
  Dataset = character(),
  Variable = character(),
  Group1 = character(),
  Group2 = character(),
  N1 = integer(),
  N2 = integer(),
  Median1 = numeric(),
  Median2 = numeric(),
  W_Statistic = numeric(),
  P_Value_Two_Sided = numeric(),
  P_Value_Directional = numeric(),
  Effect_Size = numeric(),
  Variance_P_Value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (name in names(all_results)) {
  result <- all_results[[name]]
  if (!is.null(result)) {
    # Determine which directional p-value to report based on medians
    if (result$median1 > result$median2) {
      directional_p <- result$p_greater
      direction <- paste(result$group1, ">", result$group2)
    } else {
      directional_p <- result$p_less
      direction <- paste(result$group1, "<", result$group2)
    }
    
    # Interpretation
    interpretation <- case_when(
      result$p_two_sided < 0.001 ~ "Highly significant difference",
      result$p_two_sided < 0.01 ~ "Significant difference", 
      result$p_two_sided < 0.05 ~ "Significant difference",
      result$p_two_sided < 0.1 ~ "Marginally significant",
      TRUE ~ "No significant difference"
    )
    
    summary_df <- rbind(summary_df, data.frame(
      Dataset = result$dataset,
      Variable = result$variable,
      Group1 = result$group1,
      Group2 = result$group2,
      N1 = result$n1,
      N2 = result$n2,
      Median1 = round(result$median1, 4),
      Median2 = round(result$median2, 4),
      W_Statistic = result$w_statistic,
      P_Value_Two_Sided = result$p_two_sided,
      P_Value_Directional = directional_p,
      Effect_Size = round(result$effect_size, 4),
      Variance_P_Value = result$variance_p_value,
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    ))
  }
}

print(summary_df)

# Save results
write.csv(summary_df, "nonparametric_test_results.csv", row.names = FALSE)

# =============================================================================
# 5. MANUSCRIPT-READY RESULTS
# =============================================================================
cat("\n\n", rep("#", 80), "\n")
cat("MANUSCRIPT-READY STATISTICAL REPORTING\n")
cat(rep("#", 80), "\n")

cat("\nFor gas concentrations:\n")
cat("Statistical comparisons used Mann-Whitney U tests due to non-normal distributions.\n\n")

for (gas in names(gas_results)) {
  result <- gas_results[[gas]]
  significance <- if (result$p_two_sided < 0.001) "p < 0.001" else paste("p =", round(result$p_two_sided, 3))
  
  if (result$p_two_sided < 0.05) {
    direction <- if (result$median1 > result$median2) "higher in sapwood" else "higher in heartwood"
    cat(gas, "concentrations were significantly", direction, "(Mann-Whitney U,", significance, ")\n")
  } else {
    cat(gas, "concentrations showed no significant difference between tissues (Mann-Whitney U,", significance, ")\n")
  }
}

cat("\nFor ddPCR gene abundances:\n")
for (gene in names(ddpcr_results)) {
  result <- ddpcr_results[[gene]]
  gene_name <- if (gene == "pmoa") "pmoA" else "mcrA"
  significance <- if (result$p_two_sided < 0.001) "p < 0.001" else paste("p =", round(result$p_two_sided, 3))
  
  if (result$p_two_sided < 0.05) {
    direction <- if (result$median1 > result$median2) "higher in sapwood" else "higher in heartwood" 
    cat(gene_name, "gene abundance was significantly", direction, "(Mann-Whitney U,", significance, ")\n")
  } else {
    cat(gene_name, "gene abundance showed no significant difference (Mann-Whitney U,", significance, ")\n")
  }
}

cat("\nFor microbial functions:\n")
function_names <- list(
  methanotrophy = "Methanotrophy",
  methanogenesis = "Hydrogenotrophic methanogenesis", 
  chemoheterotrophy = "Aerobic chemoheterotrophy",
  fermentation = "Fermentation"
)

for (func in names(functional_results)) {
  result <- functional_results[[func]]
  func_name <- function_names[[func]]
  significance <- if (result$p_two_sided < 0.001) "p < 0.001" else paste("p =", round(result$p_two_sided, 3))
  
  if (result$p_two_sided < 0.05) {
    direction <- if (result$median1 > result$median2) "higher in sapwood" else "higher in heartwood"
    cat(func_name, "was significantly", direction, "(Mann-Whitney U,", significance, ")\n")
  } else {
    cat(func_name, "showed no significant difference (Mann-Whitney U,", significance, ")\n")
  }
}

cat("\nVariance differences (Fligner-Killeen test):\n")
for (name in names(all_results)) {
  result <- all_results[[name]]
  if (result$variance_p_value < 0.05) {
    significance <- if (result$variance_p_value < 0.001) "p < 0.001" else paste("p =", round(result$variance_p_value, 3))
    cat(result$variable, "showed significantly different variances between tissues (", significance, ")\n")
  }
}