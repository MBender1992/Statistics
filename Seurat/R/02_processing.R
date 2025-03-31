#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

## clear workspace and memory
rm(list = ls())
gc()

#*****************************************
# 2. Install and load packages 

# load packages
library(tidyverse)
library(Seurat)
library(scales)
library(cowplot)
library(RCurl)
library(AnnotationHub)
library(ensembldb)

#*****************************************
# 3. Create custom functions

#*****************************************
# 4. Load data
list_seurat <- readRDS("~/scRNAseq/GSE171524/rds/01_GSE171524_seurat_list.rds")
doublets <- readRDS("Seurat/data/01_GSE171524_doublets.rds")

## merge seurat dataset
merged_seurat <- merge(x = list_seurat[[1]],
                       y = list_seurat[-1], 
                       add.cell.id = names(list_seurat))

## Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat$doublet <- ifelse(doublets$scDblFinder.class == "doublet", TRUE, FALSE)
merged_seurat$DoubletScore <- doublets$scDblFinder.score

## remove doublets
merged_seurat <- merged_seurat[,which(!merged_seurat$doublet)]

##******
## add relevant metadata

## calculate novelty score
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

## transform and process metadata
metadata <- merged_seurat@meta.data
metadata <- dplyr::rename(metadata, nUMI = nCount_RNA,nGene = nFeature_RNA)
metadata$Condition <- ifelse(str_detect(metadata$orig.ident, "ctr"), "Control", "Covid19")
metadata$Sample <- metadata$orig.ident
metadata$MitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")/100
metadata$Cells <- rownames(metadata)

## add metadta back to Seurat object
merged_seurat@meta.data <- metadata

## Create .Rdata object to load seurat object
saveRDS(merged_seurat, file = "~/scRNAseq/GSE171524/rds/02_GSE171524_merged_seurat.rds")

