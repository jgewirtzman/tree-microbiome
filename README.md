# Tree 16s Analysis

## Overview
This repository contains the R script `Tree_16s_all.R`, which is designed for analyzing microbial communities in tree samples. The script uses several statistical and graphical packages in R to process and visualize data derived from different sources such as ddPCR and qPCR.

## Features
- **Data Import**: Automated import of CSV data files.
- **Data Merging**: Integrates multiple data sources to create comprehensive datasets.
- **Visualization**: Uses ggplot2 and other visualization tools to create insightful plots of the processed data.

## Installation
To run this script, you will need R and the following R packages:
- phyloseq
- RColorBrewer
- ggplot2
- microeco
- pheatmap
- ggpubr
- MicrobiotaProcess
- tidyr
- GUniFrac
- FEAST
- dplyr
- umap
- phytools
- ape

You can install these packages using R commands like:
```R
install.packages("ggplot2")
install.packages("RColorBrewer")
# Add other packages similarly
