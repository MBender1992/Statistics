#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

# Have a look at supplementary from the following paper to adjust workflow
# Excerpt taken from "Single-cell RNA sequencing reveals melanoma cell state dependent
# heterogeneity of response to MAPK inhibitors" 
# Lim et al. 2024

## clear workspace and memory
rm(list = ls())
gc()

#*****************************************
# 2. Install and load packages 

## load packages
library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(scDblFinder)
library(ggprism)
library(ggsci)

#*****************************************
# 3. Create custom functions (add into ekbSeq package)

theme_prism2 <- function(base_size = 9, base_family = "", ...){
  list(theme_prism(...),  
       theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1), 
             plot.title = element_text(hjust = 0.5, face = "bold")),
       guides(x = "prism_offset", y = "prism_offset"),
       scale_color_npg(),
       scale_fill_npg())
}

#*****************************************
# 4. Load data
list_seurat <- readRDS("~/scRNAseq/GSE171524/rds/01_GSE171524_seurat_list.rds") 
## this object was generated with "make_seurat_object" stored on the workstation with locally
## stored large files which would not fit into the github folder (in the ~/scRNAseq folder)

## merge seurat dataset
merged_seurat <- merge(x = list_seurat[[1]],
                       y = list_seurat[-1], 
                       add.cell.id = names(list_seurat))

## Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat)

#*****************************************
# 5. Identification of doublets with scDblFinder
## for droplet based methods: no empty droplets should be present

## scDblFinder safeguards when examining transitional cells
## 1. Artificial doublets likely do not mimic true transitional states
## 2. Feature selection: KNN derived metrics are used rather than pure transcriptional overlap
## 3. Empirical results: <5% false positve rate for intermediate populations in differentiation datasets 

## convert seurat to sce object
sce <- as.SingleCellExperiment(merged_seurat)

## run scDblFinder for each sample separately (yields better results)
set.seed(123)
dbl <- scDblFinder(sce, samples = "orig.ident") ## specification of the samples argument performs doublet detection per sample which is recommended

## tabularize results of number of singlets and doublets
dbl_tbl <-table(dbl$scDblFinder.class)

## visualize Doublet score distribution before filtering with default threshold
p1 <- data.frame(class = dbl$scDblFinder.class, score = dbl$scDblFinder.score) %>%
  ggplot(aes(x = class, y = score, fill = class)) + 
  geom_violin() +
  geom_jitter(size = 0.5, alpha = 0.1) +
  annotate("text", x = 1, y = 0.9, label = paste0("N = ", dbl_tbl[1])) +
  annotate("text", x = 2, y = 0.2, label = paste0("N = ", dbl_tbl[2])) +
  theme_prism2() +  ylab("scDblFinder score") +
  ggtitle("Doublet score distribution")

## extract and format metadata from sce
meta_dbl <- as.data.frame(colData(dbl)) %>%
  dplyr::select(starts_with("scDblFinder")) %>%
  `rownames<-`(colnames(dbl))

## add doublet classifications to seurat object
merged_seurat <- AddMetaData(object = merged_seurat, metadata = meta_dbl[, c("scDblFinder.score","scDblFinder.class"), drop = FALSE])

#*****************************************
# 6. Visualisation of results

## split object for separate normalization
merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$orig.ident)

## run seurat pipeline
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.8)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)

## plot UMAP
p2 <- DimPlot(merged_seurat, reduction = "umap") + ggtitle("Default UMAP plot") + xlab("UMAP1") + ylab("UMAP2")
p3 <- FeaturePlot(merged_seurat, "scDblFinder.score", raster = FALSE) + scale_color_viridis_c(begin = 0, end = 1, option = "inferno") + ggtitle("UMAP colored by doublet score") + xlab("UMAP1") + ylab("UMAP2")
p4 <- ggpubr::ggarrange(p1, p2, p3, nrow = 1)
ggsave("Seurat/results/doublet_detection.png", p4, width = 300, height =100, unit = "mm")

## export doublet indices
doublets_meta <- merged_seurat@meta.data %>%
  dplyr::select(scDblFinder.score, scDblFinder.class) %>%
  rownames_to_column("UMI")

## save index indicating doublets
saveRDS(doublets_meta, "Seurat/data/01_GSE171524_doublets.rds")  


