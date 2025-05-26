# Load required libraries
library(phyloseq)
library(ggvenn)
library(patchwork)

################
# Process 16S  #
################

# Step 1: Remove mitochondria and chloroplasts
remove_taxa <- function(physeq) {
  mitochondria <- rownames(tax_table(physeq))[which(tax_table(physeq)[,5] == "Mitochondria")]
  chloroplast <- rownames(tax_table(physeq))[which(tax_table(physeq)[,4] == "Chloroplast")]
  
  bad_taxa <- c(mitochondria, chloroplast)
  return(prune_taxa(!taxa_names(physeq) %in% bad_taxa, physeq))
}

ps_clean <- remove_taxa(raw_ps)

# Step 2: Create FILTERED (but NOT rarefied) dataset
ps_filtered <- prune_taxa(taxa_sums(ps_clean) > 10, ps_clean)

# Step 3: Create FILTERED AND RAREFIED dataset
set.seed(1234)  # Ensure reproducibility
ps_rarefied <- rarefy_even_depth(ps_filtered, sample.size = 3500)

# Step 4: Extract sample metadata & confirm core_type values
metadata <- as.data.frame(sample_data(ps_clean))
print(unique(metadata$core_type))  # Confirm core_type values

# Function to extract ASVs per core type
extract_asvs <- function(ps_obj, core_type) {
  sample_ids <- rownames(metadata[metadata$core_type == core_type, ])
  ps_subset <- prune_samples(sample_ids, ps_obj)
  ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)  # Remove absent taxa
  return(unique(taxa_names(ps_subset)))
}

# Define environment names based on `core_type` in metadata
core_types <- c("Outer", "Inner", "Mineral", "Organic")

# Get ASVs per environment for filtered dataset
asvs_filtered <- list(
  Sapwood = extract_asvs(ps_filtered, "Outer"),
  Heartwood = extract_asvs(ps_filtered, "Inner"),
  Mineral = extract_asvs(ps_filtered, "Mineral"),
  Organic = extract_asvs(ps_filtered, "Organic")
)

# Get ASVs per environment for filtered AND rarefied dataset
asvs_rarefied <- list(
  Sapwood = extract_asvs(ps_rarefied, "Outer"),
  Heartwood = extract_asvs(ps_rarefied, "Inner"),
  Mineral = extract_asvs(ps_rarefied, "Mineral"),
  Organic = extract_asvs(ps_rarefied, "Organic")
)

# Calculate ASV counts for each environment
asv_counts_filtered <- sapply(asvs_filtered, length)
asv_counts_rarefied <- sapply(asvs_rarefied, length)

# Modify the titles to include ASV counts
filtered_title_16s <- paste0(
  "Filtered but NOT Rarefied 16S ASVs\n", 
  "Sapwood: ", asv_counts_filtered["Sapwood"], " | ",
  "Heartwood: ", asv_counts_filtered["Heartwood"], " | ",
  "Mineral: ", asv_counts_filtered["Mineral"], " | ",
  "Organic: ", asv_counts_filtered["Organic"]
)

rarefied_title_16s <- paste0(
  "Filtered AND Rarefied 16S ASVs\n", 
  "Sapwood: ", asv_counts_rarefied["Sapwood"], " | ",
  "Heartwood: ", asv_counts_rarefied["Heartwood"], " | ",
  "Mineral: ", asv_counts_rarefied["Mineral"], " | ",
  "Organic: ", asv_counts_rarefied["Organic"]
)

p1 <- ggvenn(asvs_filtered, 
             fill_color = c("#dfc27d", "#a6611a", "#80cdc1", "#018571"),
             stroke_size = 0.5, 
             set_name_size = 5) +
  ggtitle(filtered_title_16s)

p2 <- ggvenn(asvs_rarefied, 
             fill_color = c("#dfc27d", "#a6611a", "#80cdc1", "#018571"),
             stroke_size = 0.5, 
             set_name_size = 5) +
  ggtitle(rarefied_title_16s)

################
# Process ITS  #
################

# Step 1: Directly use the raw ITS dataset (no need to remove mitochondria/chloroplasts)
ps_clean_its <- raw_ps_its  

# Step 2: Create FILTERED (but NOT rarefied) dataset
ps_filtered_its <- prune_taxa(taxa_sums(ps_clean_its) > 10, ps_clean_its)

# Step 3: Create FILTERED AND RAREFIED dataset
set.seed(1234)  # Ensure reproducibility
ps_rarefied_its <- rarefy_even_depth(ps_filtered_its, sample.size = 3500)

# Step 4: Extract sample metadata & confirm core_type values
metadata_its <- as.data.frame(sample_data(ps_clean_its))
print(unique(metadata_its$core_type))  # Confirm core_type values

# Function to extract ASVs per core type
extract_asvs_its <- function(ps_obj, core_type) {
  sample_ids <- rownames(metadata_its[metadata_its$core_type == core_type, ])
  ps_subset <- prune_samples(sample_ids, ps_obj)
  ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)  # Remove absent taxa
  return(unique(taxa_names(ps_subset)))
}

# Get ASVs per environment for filtered dataset
asvs_filtered_its <- list(
  Sapwood = extract_asvs_its(ps_filtered_its, "Outer"),
  Heartwood = extract_asvs_its(ps_filtered_its, "Inner"),
  Mineral = extract_asvs_its(ps_filtered_its, "Mineral"),
  Organic = extract_asvs_its(ps_filtered_its, "Organic")
)

# Get ASVs per environment for filtered AND rarefied dataset
asvs_rarefied_its <- list(
  Sapwood = extract_asvs_its(ps_rarefied_its, "Outer"),
  Heartwood = extract_asvs_its(ps_rarefied_its, "Inner"),
  Mineral = extract_asvs_its(ps_rarefied_its, "Mineral"),
  Organic = extract_asvs_its(ps_rarefied_its, "Organic")
)

# Calculate ASV counts for each environment
asv_counts_filtered_its <- sapply(asvs_filtered_its, length)
asv_counts_rarefied_its <- sapply(asvs_rarefied_its, length)

# Modify the titles to include ASV counts
filtered_title_its <- paste0(
  "Filtered but NOT Rarefied ITS ASVs\n", 
  "Sapwood: ", asv_counts_filtered_its["Sapwood"], " | ",
  "Heartwood: ", asv_counts_filtered_its["Heartwood"], " | ",
  "Mineral: ", asv_counts_filtered_its["Mineral"], " | ",
  "Organic: ", asv_counts_filtered_its["Organic"]
)

rarefied_title_its <- paste0(
  "Filtered AND Rarefied ITS ASVs\n", 
  "Sapwood: ", asv_counts_rarefied_its["Sapwood"], " | ",
  "Heartwood: ", asv_counts_rarefied_its["Heartwood"], " | ",
  "Mineral: ", asv_counts_rarefied_its["Mineral"], " | ",
  "Organic: ", asv_counts_rarefied_its["Organic"]
)

p3 <- ggvenn(asvs_filtered_its, 
             fill_color = c("#dfc27d", "#a6611a", "#80cdc1", "#018571"),
             stroke_size = 0.5, 
             set_name_size = 5) +
  ggtitle(filtered_title_its)

p4 <- ggvenn(asvs_rarefied_its, 
             fill_color = c("#dfc27d", "#a6611a", "#80cdc1", "#018571"),
             stroke_size = 0.5, 
             set_name_size = 5) +
  ggtitle(rarefied_title_its)

# Combine all plots
final_venn_plot <- (p1 + p2) / (p3 + p4)
print(final_venn_plot)

# Save as A4 size .tif (210 x 297 mm = 8.27 x 11.69 inches)
ggsave("venn_diagram.tif", 
       plot = final_venn_plot,
       width = 12, 
       height = 12, 
       units = "in",
       dpi = 600)

# Save processed data for subsequent scripts
#save(ps_filtered, ps_rarefied, ps_filtered_its, ps_rarefied_its,      metadata, metadata_its, file = "2_venn_workspace.RData")