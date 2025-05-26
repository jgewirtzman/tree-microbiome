# Required packages
library(ggplot2)
library(ggridges)
library(dplyr)
library(scales)
library(viridis)

# Read the data
data <- read.csv("Other-Metadata/annotated_metadata/16S_tree_sample_table_with_meta.csv")

# Convert copies/µl to copies/g
elution_volume <- 100  # µl
data$copies_per_g <- (data$X16S_per_ul * elution_volume) / 
  (data$Sample.Mass.Added.to.Tube..mg..x / 1000)

# Set species order
species_order <- c(
  "BELE", "BEAL", "BEPA", "CAOV", "QURU", "QUVE", "QUAL", 
  "FAGR", "PRSE", "ACSA", "ACRU", "FRAM", "KALA", "SAAL", 
  "PIST", "TSCA"
)

# Transform to log10 and remove infinite/NA values
data$log10_copies <- log10(data$copies_per_g)
data <- data %>% 
  filter(!is.infinite(log10_copies), !is.na(log10_copies)) %>%
  # Filter to only include Inner/Outer and recode them
  filter(core_type %in% c("Inner", "Outer")) %>%
  mutate(
    core_type = recode(core_type,
                       "Inner" = "Heartwood",
                       "Outer" = "Sapwood"),
    species.x = factor(species.x, levels = species_order)  # Set factor order
  )

# Create species name mapping
species_names <- c(
  "FAGR" = "Fagus grandifolia",
  "BEAL" = "Betula alleghaniensis",
  "ACRU" = "Acer rubrum",
  "ACSA" = "Acer saccharum",
  "BELE" = "Betula lenta",
  "TSCA" = "Tsuga canadensis",
  "QURU" = "Quercus rubra",
  "FRAM" = "Fraxinus americana",
  "PRSE" = "Prunus serotina",
  "KALA" = "Kalmia latifolia",
  "SAAL" = "Sassafras albidum",
  "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina",
  "BEPA" = "Betula papyrifera",
  "CAOV" = "Carya ovata",
  "PIST" = "Pinus strobus"
)

# Let's print out the counts before creating the labels
print("Species counts:")
species_counts <- data %>%
  group_by(species.x) %>%
  summarise(
    total_n = n(),
    latin_name = species_names[as.character(first(species.x))]
  ) %>%
  arrange(match(species.x, species_order))

print(species_counts)

# Create the labels with counts
species_labels <- species_counts %>%
  mutate(full_label = paste0(latin_name, " (n=", total_n, ")")) %>%
  pull(full_label, species.x)

# Process the data for means
processed_data <- data %>%
  filter(!is.na(species.x)) %>%  # Remove NA species
  group_by(species.x, core_type) %>%
  summarise(
    mean_log_copies = mean(log10_copies, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

# Create base plot
p <- ggplot() +
  # Reference lines
  geom_vline(
    xintercept = c(4, 6, 8),
    color = "red",
    linetype = "dashed",
    alpha = 0.9
  ) +
  # Labels
  annotate(
    "text",
    x = c(4, 6, 8),
    y = 0,
    label = c("Human Stomach", "Seawater", "Soil"),
    angle = 0,
    vjust = -0.5,
    hjust=1.1,
    color = "red"
  ) +
  # Axis labels and theme
  labs(
    x = expression("16S" ~ Log[10](Copies~g^-1)),
    y = NULL,
    color = NULL
  ) +
  scale_y_discrete(
    labels = species_labels,
    limits = rev(species_order)  # This reverses the order
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "italic")  # Italicize species names
  ) +
  scale_x_continuous(limits = c(2, 10))


# Add points
p <- p +
  geom_point(
    data = processed_data,
    aes(x = mean_log_copies, 
        y = species.x,
        color = core_type),
    size = 3,
    shape = 21,
    stroke = 1,
    fill = NA
  ) +
  scale_color_manual(
    values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4")
  )

# Add ridges
p <- p + 
  geom_density_ridges_gradient(
    data = data %>% filter(!is.na(species.x)),
    aes(x = log10_copies, y = species.x, fill = after_stat(x)),
    scale = 0.8,
    rel_min_height = 0.01,
    bandwidth = 0.6,
    alpha = 0.4  # Lower transparency for better visibility
  )+
  scale_fill_viridis(option="D", alpha=0.6, guide="none")+theme(legend.position = 'top')

print(p)

p+ggsave("16S.tif", width = 7, height = 7, dpi = 300, units = "in")
