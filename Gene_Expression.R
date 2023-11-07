#Gene Expression: 

#Manipulate gene expression: 
getwd()
setwd("~/Desktop/Bioinformatics/Gene_expression/")
getwd()

#install packages: 
install.packages("dplyr")
install.packages("tidyverse")
install.packages('BiocManager')
BiocManager::install("GEOquery")

#loading libraries: 
library(dplyr)
library(tidyverse)
library(BiocManager)
library(GEOquery)

#read the data: 
data <- read.csv("GSE183947_fpkm.csv") #gene expression data 
dim(data)

#obtaining metadata file: 
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]])) 
head(metadata)

metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>% 
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ","",tissue)) %>% #remove word tissue
  mutate(metastasis = gsub("metastasis: ","",metastasis)) #removes metastasis

head(data)

#reshaping data: 
data.long <- data %>%
  rename(gene = X) %>%
  gather(key ='samples', value = 'FPKM', -gene)

#joining dataframes gene expression & metadata: 
data.long <- data.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

#data exploration: 
#Extract expression of BRCA1 Vs BRCA2 in Tumor Vs Normal samples: 

data.long %>%
  filter (gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM), 
            median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM) #ascending order, type - before the column name for column name for descending order 

#Visualize using ggplot2: 

library(ggplot2)

#barplot:comparing expression of BRCA1 between tissues
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = samples,y = FPKM, fill = tissue)) +
  geom_col()

#density: Compare distribution of expression between tumor and normal tissues
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.4)

#Boxplot: Compare expression between samples based on metastasis status
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis,y = FPKM)) +
  geom_boxplot()

#Scatterplot: compare expression between BRCA1 & BRCA2: 
data.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) + #compare normal and brest cancer tissue
  geom_point() +
  geom_smooth(method ='lm', se = FALSE) #Look at trend: positive slope - might have correlation 

#Heatmap: compare multiple genes across all samples
gene_interest <- c('BRCA1', 'BRCA2', 'TP53','ALK','MYCN')

hm <- data.long %>%
  filter(gene %in% gene_interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() + 
  scale_fill_gradient(low = 'yellow', high = 'black')

ggsave(hm, filename = 'heatmap1.png', width = 10, height = 8)


