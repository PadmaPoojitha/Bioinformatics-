#Single-cell RNA data analysis using Harmony: 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("harmony", version = "3.18")
install.packages("devtools")
devtools::install_github('satijalab/seurat-data', force=TRUE)
devtools::install_github('satijalab/seurat-wrappers', force=TRUE)

set.seed(1234)

#Load Libraries: 

library(harmony)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(tidyverse)
library(ggplot2)

#Load Data: 
AvailableData()

#We will use ifnb.SeuratData: 
install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source") 
library(ifnb.SeuratData)

#InstallData("ifnb")
LoadData("ifnb")
str(ifnb)
UpdateSeuratObject("ifnb.SeuratData")

#QC: 
ifnb$mitopercent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
ifnb_filter <- subset(ifnb, subset = nCount_RNA > 800 &
                      nFeature_RNA >200 &
                        mitopercent < 5)

#Normalization: 
ifnb_filter <- NormalizeData(ifnb_filter)

#Variables Features: 
ifnb_filter <- FindVariableFeatures(ifnb_filter)

#Scaling: 
ifnb_filter <- ScaleData(ifnb_filter)

#Linear dimensionality reduction: 
ifnb_filter <- RunPCA(ifnb_filter)
ElbowPlot(ifnb_filter)

#Non-liner dimensionality reduction:
ifnb_filter <- RunUMAP(ifnb_filter, dims = 1:20, reduction = 'pca')
PLot1 <- DimPlot(ifnb_filter, reduction = 'umap', group.by = 'stim')

#Run Harmony: 
ifnb_harmony <- ifnb_filter %>% 
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb_harmony@reductions

ifnb_harem <- Embeddings(ifnb_harmony, 'harmony')

#UMAP & clustering for harmony: 
ifnb_harmony <- ifnb_harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

Plot2 <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'stim')
PLot1 | Plot2

#Downstream Analysis: 

