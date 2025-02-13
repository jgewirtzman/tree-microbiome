# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)


################
# Process 16S  #
################

# Convert OTU table to matrix
otu_table_matrix <- as.matrix(otu_table(ps_rarefied))
otu_table_matrix <- t(otu_table_matrix)  # Transpose so samples are rows

# Compute Chao1
chao1_matrix <- estimateR(otu_table_matrix)
chao1_values <- chao1_matrix["S.chao1", ]
names(chao1_values) <- colnames(chao1_matrix)

# Compute Shannon diversity
shannon_values <- diversity(otu_table_matrix, index = "shannon")

# Extract metadata and ensure alignment
sample_names <- rownames(sample_data(ps_rarefied))

# Create dataframe with matched sample IDs
alpha_div <- data.frame(
  SampleID = sample_names,
  Chao1 = chao1_values[sample_names],
  Shannon = shannon_values[sample_names],
  core_type = metadata[sample_names, "core_type"]
)

# Rename core types in the data
alpha_div <- alpha_div %>%
  mutate(core_type = recode(core_type, "Inner" = "Heartwood", "Outer" = "Sapwood")) %>%
  filter(core_type %in% c("Heartwood", "Sapwood", "Organic", "Mineral"))

# Compute mean and standard deviation (SD) for each core type
summary_df <- alpha_div %>%
  group_by(core_type) %>%
  summarise(
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    SD_Chao1 = sd(Chao1, na.rm = TRUE), 
    mean_Shannon = mean(Shannon, na.rm = TRUE),
    SD_Shannon = sd(Shannon, na.rm = TRUE) 
  )%>%
  mutate(core_type = factor(core_type, levels = c("Sapwood", "Heartwood", "Mineral", "Organic")))

# Define colors for core types
CORE_COLORS <- c("Heartwood" = "#a6611a", "Sapwood" = "#dfc27d", 
                 "Mineral" = "#80cdc1", "Organic" = "#018571")

# Violin plot for 16S Chao1
p3 <- ggplot(alpha_div, aes(x = core_type, y = Chao1, 
                            color = core_type, fill = core_type)) +
  geom_violin(alpha = 0.2, position = position_dodge(width = 0.8), color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 3, alpha = 0.5) +
  geom_point(data = summary_df, aes(y = mean_Chao1, group = core_type),
             size = 5, shape = 21, color = "black",
             position = position_dodge(width = 0.8)) +
  geom_errorbar(data = summary_df, aes(y = mean_Chao1, 
                                       ymin = mean_Chao1 - SD_Chao1, 
                                       ymax = mean_Chao1 + SD_Chao1,
                                       group = core_type),
                position = position_dodge(width = 0.8),
                width = 0.2, color = "black") +
  scale_color_manual(values = CORE_COLORS) +
  scale_fill_manual(values = CORE_COLORS) +
  labs(y = "Chao1", x=NULL) +
  theme_classic() +
  theme(legend.position = "none")

# Violin plot for 16S Shannon
p4 <- ggplot(alpha_div, aes(x = core_type, y = Shannon, 
                            color = core_type, fill = core_type)) +
  geom_violin(alpha = 0.2, position = position_dodge(width = 0.8), color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 3, alpha = 0.5) +
  geom_point(data = summary_df, aes(y = mean_Shannon, group = core_type),
             size = 5, shape = 21, color = "black",
             position = position_dodge(width = 0.8)) +
  geom_errorbar(data = summary_df, aes(y = mean_Shannon, 
                                       ymin = mean_Shannon - SD_Shannon, 
                                       ymax = mean_Shannon + SD_Shannon,
                                       group = core_type),
                position = position_dodge(width = 0.8),
                width = 0.2, color = "black") +
  scale_color_manual(values = CORE_COLORS) +
  scale_fill_manual(values = CORE_COLORS) +
  labs(y = "Shannon", x=NULL) +
  theme_classic() +
  theme(legend.position = "none")

################
# Process ITS  #
################

# Calculate alpha diversity for ITS
otu_table_matrix_its <- as.matrix(otu_table(ps_rarefied_its))
otu_table_matrix_its <- t(otu_table_matrix_its)

# Compute Chao1
chao1_matrix_its <- estimateR(otu_table_matrix_its)
chao1_values_its <- chao1_matrix_its["S.chao1", ]
names(chao1_values_its) <- colnames(chao1_matrix_its)

# Compute Shannon diversity
shannon_values_its <- diversity(otu_table_matrix_its, index = "shannon")

# Extract metadata and ensure alignment
sample_names_its <- rownames(sample_data(ps_rarefied_its))

# Create dataframe with matched sample IDs
alpha_div_its <- data.frame(
  SampleID = sample_names_its,
  Chao1 = chao1_values_its[sample_names_its],
  Shannon = shannon_values_its[sample_names_its],
  core_type = metadata_its[sample_names_its, "core_type"]
)

# Rename core types in the data and filter
alpha_div_its <- alpha_div_its %>%
  mutate(core_type = recode(core_type, "Inner" = "Heartwood", "Outer" = "Sapwood")) %>%
  filter(core_type %in% c("Heartwood", "Sapwood", "Organic", "Mineral"))

# Compute summary statistics - add factor() with levels here
summary_df_its <- alpha_div_its %>%
  group_by(core_type) %>%
  summarise(
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    SD_Chao1 = sd(Chao1, na.rm = TRUE), 
    mean_Shannon = mean(Shannon, na.rm = TRUE),
    SD_Shannon = sd(Shannon, na.rm = TRUE)
  ) %>%
  mutate(core_type = factor(core_type, levels = c("Sapwood", "Heartwood", "Mineral", "Organic")))

# Violin plot for ITS Chao1
p5 <- ggplot(alpha_div_its, aes(x = core_type, y = Chao1, 
                                color = core_type, fill = core_type)) +
  geom_violin(alpha = 0.2, position = position_dodge(width = 0.8), color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 3, alpha = 0.5) +
  geom_point(data = summary_df_its, aes(y = mean_Chao1, group = core_type),
             size = 5, shape = 21, color = "black",
             position = position_dodge(width = 0.8)) +
  geom_errorbar(data = summary_df_its, aes(y = mean_Chao1, 
                                           ymin = mean_Chao1 - SD_Chao1, 
                                           ymax = mean_Chao1 + SD_Chao1,
                                           group = core_type),
                position = position_dodge(width = 0.8),
                width = 0.2, color = "black") +
  scale_color_manual(values = CORE_COLORS) +
  scale_fill_manual(values = CORE_COLORS) +
  labs(y = "Chao1", x=NULL) +
  theme_classic() +
  theme(legend.position = "none")

# Violin plot for ITS Shannon
p6 <- ggplot(alpha_div_its, aes(x = core_type, y = Shannon, 
                                color = core_type, fill = core_type)) +
  geom_violin(alpha = 0.2, position = position_dodge(width = 0.8), color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 3, alpha = 0.5) +
  geom_point(data = summary_df_its, aes(y = mean_Shannon, group = core_type),
             size = 5, shape = 21, color = "black",
             position = position_dodge(width = 0.8)) +
  geom_errorbar(data = summary_df_its, aes(y = mean_Shannon, 
                                           ymin = mean_Shannon - SD_Shannon, 
                                           ymax = mean_Shannon + SD_Shannon,
                                           group = core_type),
                position = position_dodge(width = 0.8),
                width = 0.2, color = "black") +
  scale_color_manual(values = CORE_COLORS) +
  scale_fill_manual(values = CORE_COLORS) +
  labs(y = "Shannon", x=NULL) +
  theme_classic() +
  theme(legend.position = "none")

# Add titles to the plots
p3 <- p3 + ggtitle("16S Alpha Diversity") + 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
p5 <- p5 + ggtitle("ITS Alpha Diversity") + 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))

# Create final plot combinations
pa <- p3 + p4
pb <- p5 + p6

# Stack plots and display
final_plot <- pa / pb + plot_layout(heights = c(1, 1))
print(final_plot)

# Save workspace
#save(alpha_div, alpha_div_its, summary_df, summary_df_its,      file = "3_alpha_diversity_workspace.RData")