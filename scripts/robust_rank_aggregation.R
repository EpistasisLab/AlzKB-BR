#############################################################
# Prioritization and Selection of Potential Drug Candidates
############################################################

## combine the output of RotatE link predictoon 
seed_42 <- read.csv("100trials_seed42_Alzkb-br_rotate_head_drug_prediction.csv")
seed42_rank <- seed_42[,c("rank_seed42", "head_label")]
seed_43 <- read.csv("100trials_seed43_Alzkb-br_rotate_head_drug_prediction.csv")
seed43_rank <- seed_43[,c("rank_seed43", "head_label")]
seed_44 <- read.csv("100trials_seed44_Alzkb-br_rotate_head_drug_prediction.csv")
seed44_rank <- seed_44[,c("rank_seed44", "head_label")]
seed_45 <- read.csv("100trials_seed45_Alzkb-br_rotate_head_drug_prediction.csv")
Seed45_rank <- seed_45[,c("rank_seed45", "head_label")]
seed_46 <- read.csv("100trials_seed46_Alzkb-br_rotate_head_drug_prediction.csv")
seed46_rank <- seed_46[,c("rank_seed46", "head_label")]
seed_47 <- read.csv("100trials_seed47_Alzkb-br_rotate_head_drug_prediction.csv")
seed47_rank <- seed_47[,c("rank_seed47", "head_label")]
seed_48 <- read.csv("100trials_seed48_Alzkb-br_rotate_head_drug_prediction.csv")
seed48_rank <- seed_48[,c("rank_seed48", "head_label")]
seed_49 <- read.csv("100trials_seed49_Alzkb-br_rotate_head_drug_prediction.csv")
seed49_rank <- seed_49[,c("rank_seed49", "head_label")]
seed_50 <- read.csv("100trials_seed50_Alzkb-br_rotate_head_drug_prediction.csv")
seed50_rank <- seed_50[,c("rank_seed50", "head_label")]
seed_51 <- read.csv("100trials_seed51_Alzkb-br_rotate_head_drug_prediction.csv")
seed51_rank <- seed_51[,c("rank_seed51", "head_label")]

merged <- merge(seed42_rank, seed43_rank, by = "head_label", all = TRUE)
merged_2 <- merge(merged, seed44_rank, by = "head_label", all = TRUE)
merged_3 <- merge(merged_2, Seed45_rank, by = "head_label", all = TRUE)
merged_4 <- merge(merged_3, seed46_rank, by = "head_label", all = TRUE)
merged_5 <- merge(merged_4, seed47_rank, by = "head_label", all = TRUE)
merged_6 <- merge(merged_5, seed48_rank, by = "head_label", all = TRUE)
merged_7 <- merge(merged_6, seed49_rank, by = "head_label", all = TRUE)
merged_8 <- merge(merged_7, seed50_rank, by = "head_label", all = TRUE)
merged_9 <- merge(merged_8, seed51_rank, by = "head_label", all = TRUE)

## Robust Rank Aggregation (RRA) method
install.packages("RobustRankAggreg")  
library(RobustRankAggreg)

#extract rank columns 
only_ranks_10seeds <- merged_9
only_ranks_10seeds$head_label <- paste0("Drug_", only_ranks_10seeds$head_label) # the drug IDs are numeric, change them to names by adding Drug_ to their front, this helps with the creating list for the next step

# convert columns into lists of ranked items, remove the first column (head_label) that has drug IDs and use that column to create a ranking list of those drugs across 10 seeds (remainig 10 columns)
ranked_lists <- lapply(only_ranks_10seeds[, -1], function(column) {  
  ordered_drugs <- only_ranks_10seeds$head_label[order(column)]  # Order drug names by rank  
  return(ordered_drugs)  
})
print(ranked_lists)
str(ranked_lists) 

# ensure ranked_lists are correctly formatted  
glist <- lapply(ranked_lists, as.character)  

# create rank matrix explicitly  
rmat <- rankMatrix(glist = glist, N = NA, full = FALSE)  

# perform Robust Rank Aggregation explicitly, and create Score: RRA scores    
result <- aggregateRanks(  
  glist = glist,  
  rmat = rmat,  
  N = NA,  
  method = "RRA",  
  full = FALSE,  
  exact = FALSE,  
  topCutoff = NA  
)  

## apply Bonferroni correction and calculate the corrected/adjusted p-value
n_tests <- nrow(result)  
result$Bonferroni_pvalue <- pmin(result$Score * n_tests, 1)

# apply benjamini-hochberg to caculate the corrected/adjusted p-value, this step is optional
result$BH_pvalue <- p.adjust(result$Score, method = "BH")  



