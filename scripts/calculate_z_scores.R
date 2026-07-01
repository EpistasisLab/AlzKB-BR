###########
# z-scores 
###########

# the z-score for each gene is computed across all the 231 brain regions

expr_df <- read.csv("output/Method4_allSbjs_Brainregions_noduplicatedGenes.csv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_df) <- expr_df$Gene
expr_matrix <- as.matrix(expr_df[, -1])

# calculate z-scores
row_mean <- rowMeans(expr_matrix, na.rm = TRUE)
row_sd <- apply(expr_matrix, 1, sd, na.rm = TRUE)

#calculate mean-centered values
mean_centered <- sweep(expr_matrix, 1, row_mean, "-")

#calculate z-scores
z_scores <- sweep(mean_centered, 1, row_sd, "/")

write.csv(z_scores, "231brainregions_zscores.csv")
