#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

## clear workspace and memory
rm(list = ls())
gc()

#*****************************************
# 2. Install and load packages 

## load packages
library(tidyverse)
library(Seurat)

#*****************************************
# 3. Create custom functions


#*****************************************
# 4. Load data

## load seurat object
filtered_seurat <- readRDS("~/scRNAseq/GSE171524/rds/03_GSE171524_seurat_filtered.rds")

## split object for separate normalization
filtered_seurat[["RNA"]] <- split(filtered_seurat[["RNA"]], f = filtered_seurat$orig.ident)

## log-normalize data
filtered_seurat <- NormalizeData(filtered_seurat)

## calculate cell cycle scores
filtered_seurat <- JoinLayers(filtered_seurat)
filtered_seurat <- CellCycleScoring(filtered_seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
filtered_seurat[["RNA"]] <- split(filtered_seurat[["RNA"]], f = filtered_seurat$orig.ident)

## find variable features
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)

## scale data
filtered_seurat <- ScaleData(filtered_seurat)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(filtered_seurat, layer = "scale.data"), 10)

# # plot variable features with and without labels
p1 <- VariableFeaturePlot(filtered_seurat)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
p3 <- ggpubr::ggarrange(p1,p2, nrow = 2)
ggsave("Seurat/results/variable_features.png", p3, width = 210, height =297, unit = "mm")

## save data
saveRDS(filtered_seurat, file = "~/scRNAseq/GSE171524/rds/04_GSE171524_seurat_scaled.rds")


## compare with a scTransformed model (in a second step)

## Perform scTransform (without accounting for batch. This will be done in a non-linear fashion by integration)
## Dimensionality reduction and selection of highly variable genes (HVG genes will be selected within scTransform)
## Integrate data


