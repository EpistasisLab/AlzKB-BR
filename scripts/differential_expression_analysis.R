#############################################################################################################
#Integrating Brain Region Data into AlzKB: Analysis of differential gene expression in Brain Regions (limma)
#############################################################################################################

## create differentially expressed genes (DEGs)

# used liraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
library(limma)
library(edgeR)

# we used the output file from the step Data Acquisition and Curation to create DEGs. Notes: The gene expression data were organized into a matrix in which rows represented genes and columns represented brain regions.
expr_df <- read.csv("output/Method4_allSbjs_Brainregions_noduplicatedGenes.csv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_df) <- make.unique(expr_df$Genes)
expr_df$Genes <- NULL

# create expression matrix, design matrix, use lmFit function (limma) to define a contrast comparing the Target region with the average expression across all Other brain regions.
regions.
expr <- data.matrix(expr_df)
for (i in colnames(expr_df)){
  target <- i
  
  grp <- factor(ifelse(colnames(expr) == target, "Target", "Other"),
                levels = c("Other","Target"))
  design <- model.matrix(~0 + grp)
  colnames(design) <- c("Other","Target")
  
  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, makeContrasts(TargetVsOther = Target - Other, levels = design))
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = "TargetVsOther", number = Inf, sort.by = "P", adjust.method = "BH")
  res <- cbind(Gene = rownames(res), res)
  head(res)
  write.csv(res, paste0("limma_", target, "_vs_others.csv"), row.names = FALSE)
}

## filter DEGs based on adj-p-value and LogFC (FDR < 0.05, |logFC| > 1)
fdr_cut <- 0.05  
lfc_cut <- 1  

# Folder containing the CSV files  
input_folder <- "DEGs/allsbjs_231brainRegions_DEGsandVolacanoPlot"  
output_folder <- "allsbjs_Filtered_DEGs_Output"  # Folder where output files will be saved  

# create output folder if it doesn't exist  
if (!dir.exists(output_folder)) {  
  dir.create(output_folder)  
}  

# get list of all CSV files in the input folder  
file_list <- list.files(path = input_folder, pattern = "*.csv", full.names = TRUE)  

# loop through each file and process  
for (file_path in file_list) {  
  # Read the CSV file  
  limma_df <- read.csv(file_path)  
  
  # extract the region name from the file name  
  file_name <- basename(file_path)  # Get the file name (e.g., "limma_3_vs_others.csv")  
  region_name <- gsub("limma_|_vs_others.csv", "", file_name)  # Extract the region name (e.g., "3")  
  
  # select relevant columns  
  limma_df <- limma_df[, c("Gene", "logFC", "adj.P.Val")]  
  
  # add Brain_region column with the region name  
  limma_df$Brain_region <- region_name  
  
  # apply cutoff filtering  
  limma_df_cutoff <- subset(limma_df, adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)  
  
  # save filtered data to a new CSV file  
  output_file <- paste0(output_folder, "/limma_DEGs_cutoff_", region_name, "_vs_others.csv")  
  write.csv(limma_df_cutoff, output_file, row.names = FALSE)  
  
  # print progress  
  cat("Processed and saved:", output_file, "\n")  
}  

cat("All files have been processed!\n")  

## combine filtered DEG results from all brain regions into a single table

# folder containing the filtered CSV files  
output_folder <- "allsbjs_Filtered_DEGs_Output"  

# output file name for the combined data  
combined_output_file <- "allsbjs_filtred_combined_DEGs.csv"  

# get list of all CSV files in the output folder  
filtered_files <- list.files(path = output_folder, pattern = "*.csv", full.names = TRUE)  

# initialize an empty list to store data frames  
data_list <- list()  

# loop through each file and read the data  
for (file_path in filtered_files) {  
  # Read the CSV file  
  data <- read.csv(file_path)  
  
  # Append the data frame to the list  
  data_list[[length(data_list) + 1]] <- data  
}  

# combine all data frames into one  
combined_data <- do.call(rbind, data_list)  

# save the combined data to a new CSV file  
write.csv(combined_data, combined_output_file, row.names = FALSE)  

cat("All files have been combined into:", combined_output_file, "\n")  

#notes: the output file has 24,516 significant logFC values, 12,611 DEGs, and 122 brain regions 


