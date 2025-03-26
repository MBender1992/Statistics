#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/

## clear workspace and memory
rm(list = ls())
gc()


#*****************************************
# 2. Install and load packages 

## load packages
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(reticulate)
library(scCustomize)
library(ggpubr)

# Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/bin/python")
Sys.setenv(RETICULATE_PYTHON = "~/miniforge3/envs/scvi-env/bin/python")

#*****************************************
# 3. Create custom functions


#*****************************************
# 4. Load data

## load seurat object
scaled_seurat <- readRDS("~/scRNAseq/GSE171524/rds/04_GSE171524_seurat_scaled.rds")

## reduce dimensionality
scaled_seurat <- RunPCA(scaled_seurat)
ElbowPlot(scaled_seurat, ndims = 30, reduction = "pca")

## find neighbours and clusterss
scaled_seurat <- FindNeighbors(scaled_seurat, dims = 1:30, reduction = "pca")
scaled_seurat <- FindClusters(scaled_seurat, resolution = 0.8, cluster.name = "unintegrated_clusters")
scaled_seurat <- RunUMAP(scaled_seurat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

## visualize possible confounders
## Best practice: Run UMAP before integration and plot possible confounders. Re asses after integration
colors <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
p1 <- DimPlot(scaled_seurat, reduction = "umap.unintegrated", group.by = c("orig.ident", "condition", "Phase"), cols = colors, raster = FALSE, ncol = 1)
p2 <- FeaturePlot(scaled_seurat, "mitoRatio", raster = FALSE) +  ggtitle("UMAP colored by percentage of mtRNA") + xlab("UMAP1") + ylab("UMAP2")
## Doublet score plotten? 

#*****************************************
# 5. Data integration

## integrate data with scvi and harmony
# "conda activate scvi-env" to run scVI
integrated_seurat <- IntegrateLayers(
  object = scaled_seurat, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "~/miniforge3/envs/scvi-env", verbose = TRUE
)

integrated_seurat <- IntegrateLayers(
  object = integrated_seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = TRUE
)


## process integrated scvi data
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.scvi", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 2, cluster.name = "scvi_clusters")
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")

## process integrated harmony data
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.harmony", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 2, cluster.name = "harmony_clusters")
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

## plot UMAPS for integrated data
p3 <- DimPlot(integrated_seurat,  reduction = "umap.scvi",  group.by = c("orig.ident", "condition", "Phase"),  combine = TRUE,
              label.size = 2, ncol = 1, cols = colors, raster = FALSE) 
p4 <- DimPlot(integrated_seurat,  reduction = "umap.harmony",  group.by = c("orig.ident", "condition", "Phase"),  combine = TRUE,
              label.size = 2, ncol = 1, cols = colors, raster = FALSE)

## combine plots
p_umaps <- ggarrange(p1, p3, p4, ncol = 3, labels = LETTERS[1:3])

# ## save plots
ggsave("Seurat/results/UMAP_groups_integration.png", p_umaps, width = 500, height = 297, unit = "mm")
# ggsave("Seurat/results/UMAP_gradient_unintegrated.png", p2, width = 210, height = 100, unit = "mm")
