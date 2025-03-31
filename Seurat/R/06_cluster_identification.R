#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

## workflow
## 1. Manual annotation of cell types known to be present within the sample based on marker gene annotation in PanglaoDB
## 2. Remaining clusters can be cross referenced with PanglaoDB based on their marker gene expression (findMarkers) 
## 3. Check clusters based on their functional annotation (feed top1000 genes to DAVID: https://davidbioinformatics.nih.gov/) and check whether enriched terms make sense for assigned cell type
## 4. Automated annotation with singleR 
## 5. Comparison of manual and automated approach and refinement of clusters
## 6. Cluster validation based on new analysis with findMarkers for cell identities and cross-referencing with known markers in PanglaoDB

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
library(singleR)
library(celldex)
library(ggpubr)
library(pheatmap)

BiocManager::install("celldex")
BiocManager::install("Rhdf5lib")



#*****************************************
# 3. Create custom functions

#*****************************************
# 4. Load data

## load seurat object
# integrated_seurat <- readRDS("../../scRNAseq_data/GSE171524/rds/05b_GSE171524_seurat_integrated.rds")
integrated_seurat <- readRDS("~/scRNAseq/GSE171524/rds/05b_GSE171524_seurat_integrated.rds")
integrated_seurat <- JoinLayers(integrated_seurat)
Idents(object = integrated_seurat) <- "scvi_clusters"

#*****************************************
# 5.a Find markers (manually)

## calculate cluster markers
markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)

## filter top 3 markers for each cluster and remove lncRNA as well as pseudogenes
top3_markers <- markers %>% 
  dplyr::filter(avg_log2FC > 1 & !str_detect(gene,  "LINC|AC\\d{3,}|AL\\d{3,}|AP\\d{3,}")) %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

## save top 1000 marker genes for every cluster for GO analysis
for(i in unique(markers$cluster)){
  markers %>% 
    dplyr::filter(cluster == i & avg_log2FC > 1 & !str_detect(gene,  "LINC|AC\\d{3,}|AL\\d{3,}|AP\\d{3,}")) %>% 
    arrange(desc(avg_log2FC)) %>% 
    head(1000) %>%
    .$gene %>% write.csv(paste0("Seurat/results/marker_genes/Cluster_", i, ".csv"))
}

## Dotplot showing top3 Genes for every cell identity
p <- DotPlot(integrated_seurat, features = unique(top3_markers$gene)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("Seurat/results/Dotplot_marker_genes.png", p, width = 400, height = 200, unit = "mm", bg = "white")

## identification of fibroblasts and smooth muscle cells based on the top3 markers and research within PanglaoDB
features <- c("CDH1", "AQP1", "SEC14L3", "ADH7") ## check marker genes
p_clusters <- DimPlot(integrated_seurat,  reduction = "umap.scvi", combine = TRUE, label = TRUE,  label.size = 4, cols = colors, raster = FALSE)
p_features <- FeaturePlot(integrated_seurat, features = features, reduction = "umap.scvi", max.cutoff = 0.5, raster = FALSE)
ggarrange(p_clusters, p_features)

## show whether features identified with PanglaoDB are included in markers
markers %>% dplyr::filter(gene %in% features)

## find top 3 markers for a given cluster to check with PanglaoDB
top3_markers %>% 
  dplyr::filter(cluster %in% c(20)) %>% ## enter cluster which should be examined
  .$gene %>%
  paste(collapse = ",")

#*****************************************
# 5.b Find markers (automatically)

integrated_seurat <- readRDS("~/scRNAseq/GSE171524/rds/06_GSE171524_seurat_annotated.rds") ## delete afterwards

## get reference 
ref <- HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))
## expression values are log counts

## extract counts
counts <- GetAssayData(integrated_seurat, slot = "counts")


#*****************************************
# 6. Annnotate and validate markers
predicted_cell_types <- list(
  "1" = "Macrophages",
  "2" = "Fibroblasts",
  "3" = "CD4+ T-cells", 
  "4" = "Macrophages",
  "5" = "AT1",
  "6" = "AT2",
  "7" = "AT2",
  "8" = "Plasma cells",
  "9" = "Macrophages", 
  "10" = "Monocytes",
  "11" = "Fibroblasts",
  "12" = "Fibroblasts",
  "13" = "Natural killer cells", 
  "14" = "Endothelial cells",
  "15" = "15", # airway?
  "16" = "Neuronal cells", # DAVID supports this notion
  "17" = "Endothelial cells",
  "18" = "18", ## Endothelial cells/Pericytes/Neuronal cells
  "19" = "Airway epithelial", 
  "20" = "20", # airway epithelial?
  "21" = "B-cells",
  "22" = "CD8+ T-cells", 
  "23" = "Mast cells", # unclear or neuronal cells
  "24" = "Fibroblasts",
  "25" = "Pericytes",
  "26" = "Smooth muscle cells",
  "27" = "Dendritic cells",
  "28" = "Erythroid-like",
  "29" = "Fibroblasts",
  "30" = "B-cells", 
  "31" = "31" # ?
)

# set Cell label as ident
integrated_seurat@meta.data$Cell_label <- factor(integrated_seurat@meta.data$scvi_clusters, labels = unlist(predicted_cell_types))
Idents(object = integrated_seurat) <- "Cell_label"

p_clusters <- DimPlot(integrated_seurat,  reduction = "umap.scvi", combine = TRUE, label = TRUE,  label.size = 4, cols = colors, raster = FALSE)

## save results
saveRDS(integrated_seurat, file = "~/scRNAseq/GSE171524/rds/06_GSE171524_seurat_annotated.rds")


