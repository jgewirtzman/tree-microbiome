library(tidyverse)

# Code to determine exact sample sizes for Figure 3 legend

# 1. For ddPCR data (panels 3a and 3b - pmoA and mcrA)
ddpcr <- read.csv("/Users/jongewirtzman/Downloads/ddPCR_tree_transposed_data.csv")
ddpcr_filtered <- ddpcr %>%
  mutate(core_type = dplyr::recode(core_type, 
                                   "Inner" = "heartwood", 
                                   "Outer" = "sapwood"))

# Count samples for pmoA analysis
n_pmoa <- ddpcr_filtered %>%
  filter(!is.na(pmoa_loose)) %>%
  nrow()

# Count samples for mcrA analysis  
n_mcra <- ddpcr_filtered %>%
  filter(!is.na(mcra_probe_loose)) %>%
  nrow()

cat("pmoA samples (panel 3a):", n_pmoa, "\n")
cat("mcrA samples (panel 3b):", n_mcra, "\n")

# 2. For FAPROTAX data (panels 3c, 3d, 3e, 3f)
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

# Count samples for each FAPROTAX analysis
n_methanotrophy <- faprotax_filtered %>%
  filter(!is.na(methanotrophy)) %>%
  nrow()

n_hydrogenotrophic <- faprotax_filtered %>%
  filter(!is.na(hydrogenotrophic_methanogenesis)) %>%
  nrow()

n_aerobic_chemo <- faprotax_filtered %>%
  filter(!is.na(aerobic_chemoheterotrophy)) %>%
  nrow()

n_fermentation <- faprotax_filtered %>%
  filter(!is.na(fermentation)) %>%
  nrow()

cat("Methanotrophy samples (panel 3c):", n_methanotrophy, "\n")
cat("Hydrogenotrophic methanogenesis samples (panel 3d):", n_hydrogenotrophic, "\n")
cat("Aerobic chemoheterotrophy samples (panel 3e):", n_aerobic_chemo, "\n")
cat("Fermentation samples (panel 3f):", n_fermentation, "\n")

# 3. For dumbbell plot (panel 3j)
# This uses the GC_data from your main analysis
df_filtered <- GC_data %>%
  dplyr::select(Species, Tree.No, Tissue, CH4_concentration) %>%
  filter(!is.na(Species) & !is.na(Tree.No) & !is.na(Tissue) & !is.na(CH4_concentration)) %>%
  filter(Tissue != "") %>%
  mutate(CH4_concentration = as.numeric(CH4_concentration)) %>%
  filter(!is.na(CH4_concentration)) %>%
  group_by(Species, Tree.No, Tissue) %>%
  summarize(CH4_concentration = mean(CH4_concentration), .groups = 'drop')

df_wide <- df_filtered %>%
  pivot_wider(names_from = Tissue, values_from = CH4_concentration) %>%
  filter(!is.na(heartwood) & !is.na(sapwood))

n_paired_trees <- nrow(df_wide)

cat("Paired tree measurements (panel 3j):", n_paired_trees, "\n")

# 4. For rainfall plots panels 3g and 3h (already correctly stated as n=60)
# Verify this number
gc_samples <- GC_data %>%
  filter(Tissue != "") %>%
  distinct(Species, Tree.No) %>%
  nrow()

cat("Total trees in GC analysis (panels 3g, 3h):", gc_samples, "\n")

# Summary output for easy copy-paste into legend
cat("\n=== SAMPLE SIZES FOR FIGURE LEGEND ===\n")
cat("Panel 3a (pmoA): n=", n_pmoa, " biologically independent wood samples\n")
cat("Panel 3b (mcrA): n=", n_mcra, " biologically independent wood samples\n")
cat("Panel 3c (methanotrophy): n=", n_methanotrophy, " biologically independent samples\n")
cat("Panel 3d (methanogens): n=", n_hydrogenotrophic, " biologically independent samples\n")
cat("Panel 3e (aerobic chemo): n=", n_aerobic_chemo, " biologically independent samples\n")
cat("Panel 3f (fermentation): n=", n_fermentation, " biologically independent samples\n")
cat("Panel 3g,3h (GC analysis): n=", gc_samples, " trees\n")
cat("Panel 3j (dumbbell): n=", n_paired_trees, " individual trees with paired measurements\n")