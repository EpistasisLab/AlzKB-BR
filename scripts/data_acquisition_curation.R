##########################################################################
#Integrating Brain Region Data into AlzKB: Data Acquisition and Curation
#########################################################################

# Used liraries
library(tidyr)
library(dplyr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("preprocessCore")
library(preprocessCore)

## dataset from each brain

# sbj1: Read the expression from MicroarrayExpression.csv  
maexp_178236545 <- read.csv(file = "data/allenbrain/178236545_ds/MicroarrayExpression.csv", header = FALSE)
# Read in sample annotations from SampleAnnot.csv
sampleann_178236545 <- read.csv("data/allenbrain/178236545_ds/SampleAnnot.csv")
rownames(maexp_178236545) <- maexp_178236545[,1]
maexp_178236545 <- maexp_178236545[,-1]
colnames(maexp_178236545) <- sampleann_178236545$structure_acronym

# sbj2: Read the expression from MicroarrayExpression.csv  
maexp_178238266 <- read.csv(file = "data/allenbrain/178238266_ds/MicroarrayExpression.csv", header = FALSE)
# Read in sample annotations from SampleAnnot.csv
sampleann_178238266 <- read.csv("data/allenbrain/178238266_ds/SampleAnnot.csv")
rownames(maexp_178238266) <- maexp_178238266[,1]
maexp_178238266 <- maexp_178238266[,-1]
colnames(maexp_178238266) <- sampleann_178238266$structure_acrony 

# sbj3: Read the expression from MicroarrayExpression.csv  
maexp_178238316 <- read.csv(file = "data/allenbrain/178238316_ds/MicroarrayExpression.csv", header = FALSE)
# Read in sample annotations from SampleAnnot.csv
sampleann_178238316 <- read.csv("data/allenbrain/178238316_ds/SampleAnnot.csv")
rownames(maexp_178238316) <- maexp_178238316[,1]
maexp_178238316 <- maexp_178238316[,-1]
colnames(maexp_178238316) <- sampleann_178238316$structure_acronym

# sbj4: Read the expression from MicroarrayExpression.csv  
maexp_178238359 <- read.csv(file = "data/allenbrain/178238359_ds/MicroarrayExpression.csv", header = FALSE)
# Read in sample annotations from SampleAnnot.csv
sampleann_178238359 <- read.csv("data/allenbrain/178238359_ds/SampleAnnot.csv")
rownames(maexp_178238359) <- maexp_178238359[,1]
maexp_178238359 <- maexp_178238359[,-1]
colnames(maexp_178238359) <- sampleann_178238359$structure_acronym

# sbj5: Read the expression from MicroarrayExpression.csv 
maexp_178238373<- read.csv(file = "data/allenbrain/178238373_ds/MicroarrayExpression.csv", header = FALSE)
# Read in sample annotations from SampleAnnot.csv
sampleann_178238373<- read.csv("data/allenbrain/178238373_ds/SampleAnnot.csv")
rownames(maexp_178238373) <- maexp_178238373[,1]
maexp_178238373<- maexp_178238373[,-1]
colnames(maexp_178238373) <- sampleann_178238373$structure_acronym

# sbj6: Read the expression from MicroarrayExpression.csv  
maexp_178238387<- read.csv(file = "data/allenbrain/178238387_ds/MicroarrayExpression.csv", header = FALSE)
# Read in sample annotations from SampleAnnot.csv
sampleann_178238387<- read.csv("data/allenbrain/178238387_ds/SampleAnnot.csv")
rownames(maexp_178238387) <- maexp_178238387[,1]
maexp_178238387<- maexp_178238387[,-1]
colnames(maexp_178238387) <- sampleann_178238387$structure_acronym

## merge all brains (maexp_...)

ds1_2 <- merge(maexp_178236545[,1:501], maexp_178238266[,1:470], by="row.names", all= TRUE)
row.names(ds1_2) <- ds1_2$Row.names
ds1_2 <- ds1_2[, -1]

# Merge result with maexp_178238316
ds1_2_3 <- merge(ds1_2[,1:971], maexp_178238316[,1:529], by="row.names", all= TRUE)
row.names(ds1_2_3) <- ds1_2_3$Row.names
ds1_2_3 <- ds1_2_3[, -1]

# Merge result with maexp_178238359 
ds1_2_3_4 <- merge(ds1_2_3[,1:1500], maexp_178238359[,1:363], by="row.names", all= TRUE)
row.names(ds1_2_3_4) <- ds1_2_3_4$Row.names
ds1_2_3_4 <- ds1_2_3_4[, -1]

# Merge result with maexp_178238373 
ds1_2_3_4_5 <- merge(ds1_2_3_4[,1:1863], maexp_178238373[,1:893], by="row.names", all= TRUE)
row.names(ds1_2_3_4_5) <- ds1_2_3_4_5$Row.names
ds1_2_3_4_5 <- ds1_2_3_4_5[, -1]

# Merge with maexp_178238387
ds1_2_3_4_5_6 <- merge(ds1_2_3_4_5[,1:2756], maexp_178238387[,1:946], by="row.names", all= TRUE)
row.names(ds1_2_3_4_5_6) <- ds1_2_3_4_5_6$Row.names
ds1_2_3_4_5_6 <- ds1_2_3_4_5_6[, -1]

rm(ds1_2, ds1_2_3, ds1_2_3_4, ds1_2_3_4_5, onto_178236545, onto_178238266, onto_178238316, onto_178238359, onto_178238373, onto_178238387, probes_178236545, probes_178238316, probes_178238266, probes_178238359, probes_178238373, sampleann_178236545, sampleann_178238266, sampleann_178238316, sampleann_178238359, sampleann_178238373, sampleann_178238387)


## clean and preprocess the dataset (ds1_2_3_4_5_6) and prepare it for analysis

#add rownames as a column and name it probe_id
ds1_2_3_4_5_6_probeID <- cbind(probe_id=row.names(ds1_2_3_4_5_6), ds1_2_3_4_5_6)
rownames(ds1_2_3_4_5_6_probeID) <- NULL

# Clean Brain Region Names: Remove .x, .y, etc from the column names.note example "Pr5.1" -> "Pr5"
colnames(ds1_2_3_4_5_6_probeID) <- gsub("\\..*", "", colnames(ds1_2_3_4_5_6_probeID))
ds1_2_3_4_5_6_probeID <- apply(ds1_2_3_4_5_6_probeID,2, as.numeric) # make all value numeric 
tdata_nb <- t(ds1_2_3_4_5_6_probeID) #transpose data by switching rows and columns

#create a wide format of the data
library(reshape2)
#add rownnames as a column 
tdata_nb <- cbind(rownames(tdata_nb), tdata_nb) #Add Row Names as a Column
rownames(tdata_nb) <- NULL
colnames(tdata_nb) <- tdata_nb[1,] #set column names 
tdata_nb <- tdata_nb[-1, ] 
colnames(tdata_nb)[1] <- "brain_region"
tdata_nb <- as.data.frame(tdata_nb) #Convert to Data Frame
str(tdata_nb)
sapply(tdata_nb, class)
sbj_melt <- melt(tdata_nb, id.vars = "brain_region") #Reshape the Data with melt(), notes:Keep the brain_region column as an identifier. Melt every other column.
#remove digits were added to each brain region for duplications/replication
sbj_melt$brain_region <- gsub("\\.\\d+$", "", sbj_melt$brain_region) #Clean Brain Region Names

# mean expression value  
colnames(sbj_melt) <- c("Brain_region", "probe_id", "Value")
sbj_melt$Brain_region <- as.character(sbj_melt$Brain_region)
sbj_melt$probe_id <- as.character(sbj_melt$probe_id)
sbj_melt$Value <- as.numeric(as.character(sbj_melt$Value))
sbj_melt_mean <- aggregate(Value ~ Brain_region + probe_id, sbj_melt, FUN=mean, na.rm=TRUE, na.action="na.pass")

# mapped probes to genes using the Ensembl Genome Browser
load(file = "output/p2ensg_all.RData")
sbj_melt_mean_p2ensg <- merge(sbj_melt_mean, p2ensg_orig, by.x="probe_id", by.y="probe_id", all.x=T)
sbj_melt_mean_p2ensg_mean_noDuplicatedGene <- aggregate(Value ~ Brain_region + name, sbj_melt_mean_p2ensg, FUN=mean, na.rm=TRUE, na.action="na.pass")
cleaned_sbj_melt <- sbj_melt_mean_p2ensg_mean_noDuplicatedGene %>% filter(!is.na(name), !name %in% c("None", "NA"))

#to convert a long-format data frame (cleaned_sbj_melt) back into a wide-format data frame
Sbj_Brain_genes <- cleaned_sbj_melt %>% 
  pivot_wider(
    id_cols = name, #Each unique gene becomes a row
    names_from = Brain_region, # Brain regions become columns
    values_from = Value
  )
colnames(Sbj_Brain_genes)[1] <- "Genes" 
write.csv(Sbj_Brain_genes, "Method4_allSbjs_Brainregions_noduplicatedGenes.csv", row.names = FALSE)
