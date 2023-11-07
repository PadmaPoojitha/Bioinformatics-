#Single cell RNA Seq: 

install.packages('Seurat')
install.packages('hdf5r')

#libraries: 
library(Seurat)
library(tidyverse)
library(patchwork)

#Load Data:
getwd()
setwd("/Users/poojithaalla/Desktop/Bioinformatics/Single_cell")
nsclc.data <- Read10X_h5('20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5') 
cnts <- nsclc.data$`Gene Expression`

#Create a Seurat Object: 
nsclc.seobj <- CreateSeuratObject(counts = cnts, project = "NSCLC", min.cells = 3, min.features = 200)
nsclc.seobj

#Quality Control: 
View(nsclc.seobj@meta.data)
#Mitochondrial reads %:
nsclc.seobj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seobj, pattern = "^MT-")
View(nsclc.seobj@meta.data)

#Visualize QC metrics: 
VlnPlot(nsclc.seobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships:
FeatureScatter(nsclc.seobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

FeatureScatter(nsclc.seobj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_smooth(method = 'lm')

#Filter cells that have unique feature counts over 2,500 or less than 200 
#Low-quality cells or empty droplets will often have very few genes, Cell doublets or multiplets may exhibit an aberrantly high gene count
#Filter cells that have >5% mitochondrial counts - Low-quality / dying cells often exhibit extensive mitochondrial contamination
nsclc.seobj <- subset(nsclc.seobj, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                             percent.mt < 5)

nsclc.seobj

#Normalization:
nsclc.seobj <- NormalizeData(nsclc.seobj)

#Feature selection: 
nsclc.seobj <- FindVariableFeatures(nsclc.seobj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(nsclc.seobj), 10) #top10 highly variable genes 
plot1 <- VariableFeaturePlot(nsclc.seobj) 
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Scaling: 
all.genes <- rownames(nsclc.seobj)
nsclc.seobj <- ScaleData(nsclc.seobj, features = all.genes)

#Linear Dimensionality: 
nsclc.seobj <- RunPCA(nsclc.seobj, 
                      features = VariableFeatures(object = nsclc.seobj))

#visualize top 5 PCA results: 
print(nsclc.seobj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seobj, dims = 1, cells = 500, balanced = TRUE)

#Determine dimensionality: 
ElbowPlot(nsclc.seobj)

#Cluster: 
nsclc.seobj <- FindNeighbors(nsclc.seobj, dims = 1:15)
nsclc.seobj <- FindClusters(nsclc.seobj, resolution = c(0.1,0.3, 0.5, 0.7, 1))

DimPlot(nsclc.seobj, group.by = "RNA_snn_res.0.3", label = TRUE)

#setting identity of clusters
Idents(nsclc.seobj) <- "RNA_snn_res.0.3"
Idents(nsclc.seobj)

#Non-lineraity dimensional reduction: 
reticulate::py_install(packages ='umap-learn') #UMAP
nsclc.seobj <- RunUMAP(nsclc.seobj, dims = 1:15)
DimPlot(nsclc.seobj, reduction = "umap")

#Differentially expressed markers: 
nsclc.markers <- FindAllMarkers(nsclc.seobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nsclc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

VlnPlot(nsclc.seobj, features = c("APOE", "SPP1"))

nsclc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(nsclc.seobj, features = top10$gene) + NoLegend()


