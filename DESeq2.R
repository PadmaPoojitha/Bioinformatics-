#DESeq2:
#AIM: gene expression analysis using DESeq2
BiocManager::install("DESeq2")
BiocManager::install("airway")

library(DESeq2)
library(tidyverse)
library(airway)

getwd()
setwd("/Users/poojithaalla/Desktop/Bioinformatics/DESeq2")

#Getting data: 
data(airway)
airway
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

counts_data <- read.csv('counts_data.csv')
col_data <- read.csv('sample_info.csv')

#Check that sample names match in both files: 
all(colnames(counts_data) %in% rownames(col_data))
all(colnames(counts_data) == rownames(col_data))

#Object: 
data_mat <- DESeqDataSetFromMatrix(countData = counts_data,
                                   colData = col_data,
                                   design = ~ dexamethasone)
data_mat

new_dat <- rowSums(counts(data_mat)) >= 10
data_mat <- data_mat[new_dat,]

#Setting factor level: 
data_mat$dexamathasone <- relevel(data_mat$dexamethasone, ref = 'untreated')
data_mat <- DESeq(data_mat)
data_results <- results(data_mat)
data_results

summary(data_results)
data_results1 <- results(data_mat, alpha =0.01)
summary(data_results1)

#Visualize:
plotMA(data_results)
