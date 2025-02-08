# Tree Microbiome Analysis

## Overview
This repository contains the R script `Tree_microbiome_all.R` designed for analyzing microbial communities within tree samples. The script utilizes a comprehensive suite of statistical and graphical packages in R to process, analyze, and visualize high-throughput sequencing data, focusing on ecological and microbiological insights.

## Features
- **Data Import and Preprocessing**: Import of data files and preprocessing to ensure data quality and compatibility for analysis.
- **Microbial Community Analysis**: Includes community composition analysis, diversity assessments, and microbial abundances.
- **Microbial Source Tracking**: Uses FEAST to infer microbial community source/sink dynamics.
- **Bioinformatic Analysis**: Processes sequencing data using tools tailored for microbial ecology.
- **Ordinations**: Implements ordination methods to explore complex microbial community structures and relationships, such as PCA (Principal Component Analysis) and NMDS (Non-metric Multidimensional Scaling).
- **Statistical Testing**: Applies statistical tests to compare microbial communities across different conditions or treatments, leveraging packages for ecological statistics.

## Installation
Ensure that R is installed on your system, along with the following R packages:
```R
# Install R packages
install.packages(c("phyloseq", "RColorBrewer", "ggplot2", "microeco", "pheatmap",
                   "ggpubr", "MicrobiotaProcess", "tidyr", "GUniFrac", "FEAST",
                   "dplyr", "umap", "phytools", "ape"))

```

## Usage
Update the paths to the data files in the script to match your local setup. After updating the file paths, you can run the script in an R environment to perform the analysis.

## Demo
The code uses previously published packages for all analyses. Install times and run times vary by package and analysis, but range from seconds to hours on a standard desktop computer. For assistance running code or for more information on specific analyses, please feel free to contact the authors.

## License
This project is open source and available under the MIT License.
