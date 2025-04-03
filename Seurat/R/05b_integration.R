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
library(SeuratWrappers)
library(reticulate)
library(scCustomize)
library(ggpubr)

# set python environment
Sys.setenv(RETICULATE_PYTHON = "~/miniforge3/envs/scvi-env/bin/python")

#*****************************************
# 3. Create custom functions


#*****************************************
# 4. Load data

## load seurat object
reduced_seurat <- readRDS("../../scRNAseq_data/GSE171524/rds/05a_GSE171524_seurat_reduced.rds")
colors <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")

#*****************************************
# 5. Data integration

## integrate data with scvi and harmony
# "conda activate scvi-env" to run scVI
integrated_seurat <- IntegrateLayers(
  object = reduced_seurat, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "~/miniforge3/envs/scvi-env", verbose = TRUE
)

integrated_seurat <- IntegrateLayers(
  object = integrated_seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = TRUE
)

## reorder clusters from scvi clustering 
integrated_seurat$scvi_clusters <- factor(integrated_seurat$scvi_clusters, levels = as.character(sort(as.numeric(levels(integrated_seurat$scvi_clusters))))) ## reorder levels

## process integrated scvi data
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.scvi", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 1, cluster.name = "scvi_clusters", algorithm = 4)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")

## process integrated harmony data
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.harmony", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 1, cluster.name = "harmony_clusters", algorithm = 4)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

## plot UMAPS for integrated data
p1 <- DimPlot_scCustom(integrated_seurat, reduction = "umap.unintegrated", group.by = c("Sample", "unintegrated_clusters", "Condition", "Phase"), combine = TRUE,
              raster = FALSE, num_columns = 1)
p2 <- DimPlot_scCustom(integrated_seurat,  reduction = "umap.scvi",  group.by = c("Sample", "scvi_clusters", "Condition", "Phase"),  combine = TRUE,
              num_columns = 1, raster = FALSE) 
p3 <- DimPlot_scCustom(integrated_seurat,  reduction = "umap.harmony",  group.by = c("Sample", "harmony_clusters", "Condition", "Phase"),  combine = TRUE,
              num_columns = 1, raster = FALSE)

## plot feature plots for integrated data
p4 <- FeaturePlot_scCustom(integrated_seurat, c("MitoRatio", "DoubletScore"), raster = FALSE, combine = TRUE, num_columns = 1) 
p5 <- FeaturePlot_scCustom(integrated_seurat, c("MitoRatio", "DoubletScore"), reduction = "umap.scvi", raster = FALSE, combine = TRUE, num_columns  = 1)
p6 <- FeaturePlot_scCustom(integrated_seurat, c("MitoRatio", "DoubletScore"), reduction = "umap.harmony", raster = FALSE, combine = TRUE, num_columns = 1) 

## combine plots
group_umaps <- ggarrange(p1, p2, p3, ncol = 3, labels = LETTERS[1:3])
gradient_umaps <- ggarrange(p4, p5, p6, ncol = 3, labels = LETTERS[1:3])

# ## save plots
ggsave("Seurat/results/UMAP_groups_integration.png", group_umaps, width = 500, height = 297, unit = "mm")
ggsave("Seurat/results/UMAP_gradient_integration.png", gradient_umaps, width = 500, height = 400, unit = "mm")

## save data
saveRDS(integrated_seurat, file = "../../scRNAseq_data/GSE171524/rds/05b_GSE171524_seurat_integrated.rds")
