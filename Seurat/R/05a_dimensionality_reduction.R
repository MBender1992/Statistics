#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

## clear workspace and memory
rm(list = ls())
gc()
set.seed(2521)

#*****************************************
# 2. Install and load packages 

## load packages
library(tidyverse)
library(Seurat)
library(scCustomize)
library(ggpubr)

#*****************************************
# 3. Create custom functions


#*****************************************
# 4. Load data

## load seurat object
scaled_seurat <- readRDS("../../scRNAseq_data/GSE171524/rds/04_GSE171524_seurat_scaled.rds")
colors <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")

#*****************************************
# 5. Perform dimensionality reduction and clustering

## reduce dimensionality
scaled_seurat <- RunPCA(scaled_seurat)
ElbowPlot(scaled_seurat, ndims = 30, reduction = "pca")

## find neighbours and clusterss
scaled_seurat <- FindNeighbors(scaled_seurat, dims = 1:30, reduction = "pca")
scaled_seurat <- FindClusters(scaled_seurat, resolution = 1, cluster.name = "unintegrated_clusters", algorithm = 4)
scaled_seurat <- RunUMAP(scaled_seurat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

## visualize possible confounders
## Best practice: Run UMAP before integration and plot possible confounders. Re asses after integration
DimPlot_scCustom(scaled_seurat, reduction = "umap.unintegrated", group.by = c("Sample", "unintegrated_clusters", "Condition", "Phase"),
                 raster = FALSE, num_columns =  = 1)

FeaturePlot_scCustom(scaled_seurat, c("MitoRatio", "DoubletScore"), raster = FALSE, combine = TRUE, num_columns = 1) & 
  scale_color_viridis_c(option = "inferno", direction = 1)

## save data
saveRDS(scaled_seurat, file = "../../scRNAseq_data/GSE171524/rds/05a_GSE171524_seurat_reduced.rds")
