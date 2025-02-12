#procrustes analysis
#run after Tree_microbiome_all.R

#micro_unifrac=microbiome distance matrix
#tree_dist=tree species distance matrix

z_score_normalize <- function(mat) {
  mat <- as.matrix(mat)  # Ensure it's a matrix
  (mat - mean(mat)) / sd(mat)
}

# Apply Z-score normalization to each matrix
micro_unifrac_z <- z_score_normalize(micro_unifrac)
tree_dist_z <- z_score_normalize(tree_dist)

# Run the procrustes test with normalized matrices
proc_result <- vegan::protest(as.dist(micro_unifrac_z), as.dist(tree_dist_z))
print(proc_result)
