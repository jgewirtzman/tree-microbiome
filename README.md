# Tree Microbiome Analysis

## Overview
This repository contains scripts and datasets for analyzing microbial communities within tree samples. The primary script, Tree_microbiome_all.R, and various additional scripts perform high-throughput sequencing data processing, statistical analysis, and visualization of ecological and microbiological insights.

## Data Availability

**Important:** The large sequencing data files are hosted on Figshare due to their size. You'll need to download them separately before running the analysis scripts.

### Required Data Download

1. **Visit the Figshare dataset:**
   [A diverse and distinct microbiome inside living trees](https://figshare.com/s/b781d3c5c24f1089ed5d?)

2. **Download the following files and extract them to your local repository:**
   - `Forest_Microbiome.zip` (14.1 MB) → Extract to create `16S/` and `ITS/` directories
   - `Black_Oak.zip` (1.38 MB) → Extract to create `Black Oak/16S/` and `Black Oak/ITS/` directories
   - `GC_Data.zip` (84.02 kB) → Extract to `Downloads/` directory

3. **Your final directory structure should look like:**
   ```
   tree-microbiome/
   ├── 16S/                     # From Forest_Microbiome.zip
   │   ├── OTU_table.txt        # OTU abundance table
   │   ├── rooted_tree.nwk      # Phylogenetic tree (rooted)
   │   ├── unrooted_tree.nwk    # Phylogenetic tree (unrooted)
   │   ├── seq.fasta            # FASTA sequence file
   │   ├── taxa.biom            # Taxonomic data in BIOM format
   │   ├── taxonomy_table.txt   # Taxonomy assignments
   │   └── tree_16s_mapping_dada2_corrected.txt.no_gz
   ├── ITS/                     # From Forest_Microbiome.zip
   │   └── (Same structure as 16S)
   ├── Black Oak/               # From Black_Oak.zip
   │   ├── 16S/                 # 16S data for Black Oak
   │   └── ITS/                 # ITS data for Black Oak
   ├── Downloads/               # From GC_Data.zip (gas chromatography data)
   └── [rest of repository files]
   ```

## Repository Structure

```
tree-microbiome/
├── 16S/                     # 16S rRNA gene sequencing data (DOWNLOAD FROM FIGSHARE)
├── ITS/                     # Internal Transcribed Spacer (ITS) sequencing data (DOWNLOAD FROM FIGSHARE)
├── Black Oak/               # Black Oak-specific analyses (DOWNLOAD FROM FIGSHARE)
├── Downloads/               # Gas chromatography data (DOWNLOAD FROM FIGSHARE)
├── Other-Metadata/          # Metadata and supplementary files
│   ├── 16s_w_metadata.csv   # Metadata for 16S
│   ├── ITS_w_metadata.csv   # Metadata for ITS
│   ├── annotated_metadata/  # Processed metadata files
│   └── ddPCR_meta_all_data.csv
├── revision_scripts/        # Various analysis scripts
│   ├── gas_metabolisms/     # Gas metabolism analyses
│   ├── group_distances/     # Community distance calculations
│   ├── incubation/          # Incubation gas processing scripts
│   ├── methane_cycling/     # Methanotrophic/methanogenic abundance
│   ├── myco_phylo_differences/ # Fungal guild and phylogenetic heatmaps
│   ├── reads_asvs/          # Alpha diversity and OTU analyses
│   ├── tissue_gases/        # Tissue gas composition analyses
│   ├── old/                 # Legacy scripts
│   └── outputs/             # Figures and processed data outputs
├── S.PhyloMaker-master/     # Phylogenetic tree reconstruction tools
├── Tree_microbiome_all.R    # Main analysis script
├── metabolisms/             # Functional predictions (weighted/unweighted)
├── LICENSE                  # License file
└── README.md                # This document
```

## Features

### Data Processing and Preprocessing
- Import and format microbial sequencing data (16S, ITS)
- Quality filtering and OTU/ASV table generation
- Metadata integration for ecological interpretation

### Community Composition and Diversity
- **Alpha Diversity:** Shannon, Simpson, and Faith's PD indices (alpha_diversity.R)
- **Beta Diversity:** Bray-Curtis and UniFrac distance calculations (group_distances/)
- **Ordinations:** PCA, NMDS, and UMAP for microbial community visualization

### Statistical and Functional Analyses
- **Differential Abundance Testing:** Identifies significant microbial taxa (barplot_lme_optimized.R)
- **Microbial Source Tracking:** Infers microbial origins using FEAST
- **Metabolic Function Prediction:**
  - FAPROTAX-based functional assignments (gas_metabolisms/faprotax.R)
  - Methanotroph and methanogen abundance analysis (methane_cycling/black_oak_methanome_abundance.R)
- **Tissue Gas Composition:**
  - Process and visualize tissue-level gas data (tissue_gases/)
  - Methane and oxygen concentration analysis

### Visualization and Outputs
- **Heatmaps:**
  - Taxonomic composition (genus_heatmap.R, class_heatmap.R)
  - Fungal guild composition (myco_phylo_differences/guild_heatmaps_its.R)
- **Publication-ready Figures:**
  - Generated in outputs/ (e.g., funguild_heatmap.pdf, rainfall_plot.png)

## Installation

1. **Clone this repository:**
   ```bash
   git clone https://github.com/jgewirtzman/tree-microbiome.git
   cd tree-microbiome
   ```

2. **Download required data from Figshare** (see Data Availability section above)

3. **Install R and required packages:**
   ```r
   # Install required packages
   install.packages(c("phyloseq", "RColorBrewer", "ggplot2", "microeco", "pheatmap",
                      "ggpubr", "MicrobiotaProcess", "tidyr", "GUniFrac", "FEAST",
                      "dplyr", "umap", "phytools", "ape", "vegan"))
   
   # Install Bioconductor packages
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install("phyloseq")
   ```

## Usage

1. **Prepare Data:** 
   - Download and extract all Figshare files as described above
   - Update paths to sequencing files and metadata in scripts if needed

2. **Run the Analysis:**
   - Execute `Tree_microbiome_all.R` for full pipeline processing
   - Use sub-scripts in `revision_scripts/` for specific analyses

3. **Visualize and Interpret Results:**
   - Check `outputs/` for plots, heatmaps, and statistical summaries

## Demo

Example dataset processing times vary depending on system specifications. For questions or support, contact the repository maintainers.

## Citation

If you use this code or data, please cite:

**Manuscript:** "A diverse and distinct microbiome inside living trees" (Manuscript #2024-05-11045B)

**Data:** Available at [Figshare dataset](https://figshare.com/s/b781d3c5c24f1089ed5d?)

## License

This project is open-source and available under the MIT License.
