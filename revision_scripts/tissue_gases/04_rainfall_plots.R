# Load required libraries
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

library(tidyverse)
library(cowplot)
library(ggplot2)
library(dplyr)

# Process GC data for O2 and CH4
df_filtered <- GC_data %>%
  filter(Tissue != "") %>%
  mutate(
    O2_Area = O2_concentration/1000000*100,  # Convert to percentage
    CH4_Area = CH4_concentration,
    log10ch4 = log10(CH4_Area)
  )

# Import ddPCR data
ddpcr <- read.csv("/Users/jongewirtzman/Downloads/ddPCR_tree_transposed_data.csv")
ddpcr <- ddpcr %>%
  mutate(core_type = dplyr::recode(core_type, 
                                   "Inner" = "heartwood", 
                                   "Outer" = "sapwood"))

# Import FAPROTAX data
faprotax <- read.csv("/Users/jongewirtzman/Downloads/metabolisms/16S_metabolisms_weighted.csv")
names(faprotax)[1] <- "SampleID"

# Process FAPROTAX data
faprotax <- faprotax %>%
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

# Create O2 plots
density_plot <- ggplot(df_filtered, aes(x = O2_Area, fill = Tissue)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(df_filtered, aes(x = O2_Area, y = Tissue, fill = Tissue)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab(expression(O[2] ~ "(%)"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
PA <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

# Create CH4 plots
density_plot <- ggplot(df_filtered, aes(x = log10ch4, fill = Tissue)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(df_filtered, aes(x = log10ch4, y = Tissue, fill = Tissue)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab(expression(log[10] ~ "(CH"[4] ~ "ppm)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
PB <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))



#calculate pmoA per g wood

ddpcr$pmoa_loose<-ddpcr$pmoa_loose/ddpcr$Sample.Mass.Added.to.Tube..mg.*(75/2)*1000

# Create pmoA plots
density_plot <- ggplot(ddpcr, aes(x = log10(pmoa_loose), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(ddpcr, aes(x = log10(pmoa_loose), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab(expression(log[10] ~ "(pmoA)"))

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
combined_plot1 <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

ddpcr$mcra_probe_loose<-ddpcr$mcra_probe_loose/ddpcr$Sample.Mass.Added.to.Tube..mg.*(75/2)*1000

# Create mcrA plots
density_plot <- ggplot(ddpcr, aes(x = log10(mcra_probe_loose), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(ddpcr, aes(x = log10(mcra_probe_loose), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab(expression(log[10] ~ "(mcrA)"))

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
combined_plot2 <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

# Create aerobic chemoheterotrophy plots
density_plot <- ggplot(faprotax, aes(x = aerobic_chemoheterotrophy, fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(faprotax, aes(x = aerobic_chemoheterotrophy, y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab("% aerobic chemoheterotrophy")

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
combined_plot3 <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

# Create fermentation plots
density_plot <- ggplot(faprotax, aes(x = fermentation, fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity', bw = 5) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(faprotax, aes(x = fermentation, y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab("% fermentation")

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
combined_plot4 <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

# Create hydrogenotrophic methanogenesis plots
density_plot <- ggplot(faprotax, aes(x = log10(hydrogenotrophic_methanogenesis), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(faprotax, aes(x = log10(hydrogenotrophic_methanogenesis), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab(expression(log[10] ~ "(% " * H[2] * " methanogenesis)"))

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
combined_plot5 <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

# Create methanotrophy plots
density_plot <- ggplot(faprotax, aes(x = log10(methanotrophy), fill = core_type)) +
  geom_density(alpha = 0.7, adjust = 1.5, position = 'identity') +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  labs(fill = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

box_jitter_plot <- ggplot(faprotax, aes(x = log10(methanotrophy), y = core_type, fill = core_type)) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c("heartwood" = "#a6611a", "sapwood" = "#1f78b4")) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  ) +
  xlab(expression(log[10] ~ "(% methanotrophy)"))

aligned <- align_plots(density_plot, box_jitter_plot, align = "v", axis = "tblr")
p1 <- ggdraw(aligned[[1]])
p2 <- ggdraw(aligned[[2]])
combined_plot6 <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

# Create final rainfall plot combining all plots
rainfall <- plot_grid(
  combined_plot1, combined_plot2,
  combined_plot6, combined_plot5,
  combined_plot3, combined_plot4,
  PA, PB,
  nrow = 4
)


rainfall

# Save the plot
ggsave("rainfall_plot.png", rainfall, width = 5.33, height = 8)



# Create final rainfall plot combining all plots with labels
rainfall_labeled <- plot_grid(
  combined_plot1, combined_plot2,
  combined_plot6, combined_plot5,
  combined_plot3, combined_plot4,
  PA, PB,
  nrow = 4,
  labels = c("a", "b", "c", "d", "e", "f", "g", "h"),  # Add labels
  label_size = 14,  # Adjust size as needed
  label_x = 0.00,   # Adjust horizontal position
  label_y = 0.95,   # Adjust vertical position
  hjust = 0,        # Left align the labels
  vjust = 1         # Top align the labels
)

rainfall_labeled

# Save the plot
ggsave("rainfall_plot.png", rainfall, width = 5.33, height = 8)




# Perform Welch's t-tests for mean differences between heartwood and sapwood

# O2 concentration
t_test_O2 <- t.test(O2_Area ~ Tissue, data = df_filtered, var.equal = FALSE)

# CH4 concentration (log-transformed)
t_test_CH4 <- t.test(log10ch4 ~ Tissue, data = df_filtered, var.equal = FALSE)

# pmoA abundance (log-transformed) - Remove NAs before the test
t_test_pmoA <- t.test(log10(pmoa_loose) ~ core_type, 
                      data = ddpcr %>% filter(!is.na(pmoa_loose) & pmoa_loose > 0), 
                      var.equal = FALSE)

# mcrA abundance (log-transformed) - Remove NAs before the test
t_test_mcrA <- t.test(log10(mcra_probe_loose) ~ core_type, 
                      data = ddpcr %>% filter(!is.na(mcra_probe_loose) & mcra_probe_loose > 0), 
                      var.equal = FALSE)

# Aerobic chemoheterotrophy - Remove NAs before the test
t_test_aerobic_chemo <- t.test(aerobic_chemoheterotrophy ~ core_type, 
                               data = faprotax %>% filter(!is.na(aerobic_chemoheterotrophy)), 
                               var.equal = FALSE)

# Fermentation - Remove NAs before the test
t_test_fermentation <- t.test(fermentation ~ core_type, 
                              data = faprotax %>% filter(!is.na(fermentation)), 
                              var.equal = FALSE)

# Hydrogenotrophic methanogenesis (log-transformed) - Remove NAs and ensure values > 0
t_test_hydrogen_methano <- t.test(log10(hydrogenotrophic_methanogenesis) ~ core_type, 
                                  data = faprotax %>% filter(!is.na(hydrogenotrophic_methanogenesis) & hydrogenotrophic_methanogenesis > 0), 
                                  var.equal = FALSE)

# Methanotrophy (log-transformed) - Remove NAs and ensure values > 0
t_test_methanotrophy <- t.test(log10(methanotrophy) ~ core_type, 
                               data = faprotax %>% filter(!is.na(methanotrophy) & methanotrophy > 0), 
                               var.equal = FALSE)


# Print results
list(
  O2_t_test = t_test_O2,
  CH4_t_test = t_test_CH4,
  pmoA_t_test = t_test_pmoA,
  mcrA_t_test = t_test_mcrA,
  Aerobic_Chemo_t_test = t_test_aerobic_chemo,
  Fermentation_t_test = t_test_fermentation,
  Hydrogen_Methano_t_test = t_test_hydrogen_methano,
  Methanotrophy_t_test = t_test_methanotrophy
)


