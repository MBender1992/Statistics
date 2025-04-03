#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

## workflow
## 1. Manual annotation of cell types known to be present within the sample based on marker gene annotation in PanglaoDB
## 2. Remaining clusters can be cross referenced with PanglaoDB based on their marker gene expression (findMarkers) 
## 3. Check clusters based on their functional annotation (feed top1000 genes to DAVID: https://davidbioinformatics.nih.gov/) and check whether enriched terms make sense for assigned cell type
## 4. Automated annotation with singleR 
## 5. Comparison of manual and automated approach and refinement of clusters
## 6. Cluster validation based on new analysis with findMarkers for cell identities and cross-referencing with known markers in PanglaoDB

## for further reference see the publication that generated this data:
## https://www.nature.com/articles/s41586-021-03569-1?fromPaywallRec=false
## they propose a structure into cell_type_main, cell_type_intermediate and cell_type_fine
## first cluster coarsely, then inspect each subcluster separately and annotate cells accordingly
## to exactly match the cell types as assigned in the paper the cluster resolution would need to be increased. Also subclusters (e.g. T-cells)
## could be more closely examined

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
library(SingleR)
library(ggpubr)
library(pheatmap)
library(BiocParallel)

#*****************************************
# 3. Load data

## load seurat object
integrated_seurat <- readRDS("../../scRNAseq_data/GSE171524/rds/05b_GSE171524_seurat_integrated.rds")
integrated_seurat <- JoinLayers(integrated_seurat)
Idents(object = integrated_seurat) <- "scvi_clusters"

#*****************************************
# 4.a Find markers (manually)

## calculate cluster markers
markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1) ## based on identss
markers <- markers %>%
  dplyr::filter(!str_detect(gene,  "LINC|AC\\d{3,}|AL\\d{3,}|AP\\d{3,}")) %>%
  mutate(score = (avg_log2FC * pct.1) / pct.2)

## filter top 3 markers for each cluster and remove lncRNA as well as pseudogenes
top3_markers <- markers %>% 
  group_by(cluster) %>%
  top_n(n = 3, wt = score)

## save top 1000 marker genes for every cluster for GO analysis
for(i in unique(markers$cluster)){
  markers %>% 
    arrange(desc(score)) %>% 
    head(1000) %>%
    .$gene %>% 
    write.csv(paste0("Seurat/results/marker_genes/Cluster_", i, ".csv"))
}

## Dotplot showing top3 Genes for every cell identity
p <- DotPlot(integrated_seurat, features = unique(top3_markers$gene)) + theme(axis.text.x = element_text(angle = 90))
ggsave("Seurat/results/Dotplot_marker_genes.png", p, width = 400, height = 200, unit = "mm", bg = "white")

## identification of fibroblasts and smooth muscle cells based on the top3 markers and research within PanglaoDB
features <- c("MS4A7", "OLR1", "C1QA") ## check marker genes
p_clusters <- DimPlot_scCustom(integrated_seurat, group.by = ("scvi_clusters"),  reduction = "umap.scvi", combine = TRUE, label = TRUE,  label.size = 4, raster = FALSE)
p_features <- FeaturePlot_scCustom(integrated_seurat, features = features, reduction = "umap.scvi",  raster = FALSE, , alpha_exp = 0.75)
ggarrange(p_clusters, p_features, align = "hv")



## show whether features identified with PanglaoDB are included in markers
markers %>% dplyr::filter(gene %in% features)
markers %>%
  arrange(desc(score)) %>%
  dplyr::filter(cluster == c(31)) %>%
  head(10) %>%
  .$gene %>%
  paste(collapse = ",")

## assign manual cell types
manual_cell_types <- list(
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

# add manual labels
integrated_seurat$manual_labels <- factor(integrated_seurat$scvi_clusters, labels = unlist(manual_cell_types))
integrated_seurat$manual_labels <- factor(integrated_seurat$manual_labels, levels = sort(levels(integrated_seurat$manual_labels)))

#*****************************************
# 4.b Find markers (automatically)

## get reference 
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))
## expression values are log counts

## extract counts
counts <- GetAssayData(integrated_seurat, layer = "counts")

## get predictions
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main, BPPARAM = MulticoreParam(12)) # use 12 cores
keep_labels <- names(table(pred$labels)[table(pred$labels) > 500]) ## remove labels with less than 500 cells for easier comparison
pred <- pred[pred$labels %in% keep_labels,]

## assign SingleR labels
integrated_seurat$singleR_labels <- pred$labels[match(rownames(integrated_seurat@meta.data), rownames(pred))]

## plot annotation results
p_manual <- DimPlot_scCustom(integrated_seurat, reduction = "umap.scvi", group.by = ("manual_labels"), label = TRUE, repel = TRUE, raster = FALSE) 
p_singleR <- DimPlot_scCustom(integrated_seurat, reduction = "umap.scvi", group.by = ("singleR_labels"), raster = FALSE)
p_annotation <- ggarrange(p_clusters, p_manual, p_singleR, ncol = 1, align = "hv")

## run diagnostics
p_deltas <- plotDeltaDistribution(pred, labels.use = keep_labels) ## see which cells were confidently called from the beginning and which cells took more refinement

## compare to unsupervised clustering
tbl <- table(Assigned = integrated_seurat$singleR_labels, Clusters = integrated_seurat$scvi_clusters[])
p_h1 <- plotScoreHeatmap(pred, cells.use = which(pred$labels %in% keep_labels), labels.use = keep_labels) ## careful interpretation when data has been pruned
p_h2 <- pheatmap(log10(tbl+10), color = colorRampPalette(c("white", "blue"))(10))
p_heatmaps <- ggarrange(ggplotify::as.ggplot(p_h1), ggplotify::as.ggplot(p_h2), ncol = 1, labels = LETTERS[1:2]) + bgcolor("white")        

# save diagnostics plots
ggsave("Seurat/results/SingleR_diag_UMAP_annotation.png", p_annotation, width = 300, height = 400, unit = "mm")
ggsave("Seurat/results/SingleR_diag_Delta_Distribution.png", p_deltas, width = 210, height = 297, unit = "mm")
ggsave("Seurat/results/SingleR_diag_Heatmaps.png", p_heatmaps, width = 210, height = 297, unit = "mm")

#*****************************************
# 4.c Refine labels
## assign manual cell types
refined_cell_types <- list(
  "1" = "Macrophages",
  "2" = "Fibroblasts",
  "3" = "T-cells", 
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
  "15" = "Airway epithelial", 
  "16" = "Airway epithelial", 
  "17" = "Endothelial cells",
  "18" = "Neuronal cells", 
  "19" = "Airway epithelial", 
  "20" = "Fibroblasts", 
  "21" = "B-cells",
  "22" = "Cycling NK_T-cells", 
  "23" = "Mast cells", 
  "24" = "Fibroblasts",
  "25" = "Smooth muscle cells",
  "26" = "Fibroblasts",
  "27" = "Dendritic cells",
  "28" = "Erythroid-like",
  "29" = "Fibroblasts",
  "30" = "B-cells", 
  "31" = "Fibroblasts" # ?
)

## assign cell labels labels
integrated_seurat$Cell_label <- factor(integrated_seurat$scvi_clusters, labels = unlist(refined_cell_types))
integrated_seurat$Cell_label <- factor(integrated_seurat$Cell_label, levels = sort(levels(integrated_seurat$Cell_label)))

#*****************************************
# 5. Validate markers
Idents(object = integrated_seurat) <- "Cell_label"
markers_validation <- FindAllMarkers(integrated_seurat, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25) ## based on identss
# markers_validation_orig <- markers_validation
# markers_validation <- markers_validation_orig

markers_validation <- markers_validation %>%
  dplyr::filter(!str_detect(gene,  "LINC|AC\\d{3,}|AL\\d{3,}|AP\\d{3,}")) %>%
  mutate(score = (avg_log2FC * pct.1) / pct.2)

## filter top 3 marker genes for each cluster
top3_markers_validation <- markers_validation %>% 
  group_by(cluster) %>%
  top_n(n = 3, wt = score)

## plot UMAP and extract cluster names
p_umap <- DimPlot_scCustom(integrated_seurat, group.by = ("Cell_label"),  reduction = "umap.scvi", combine = TRUE, 
                           label = TRUE, repel = TRUE, label.size = 4, raster = FALSE)
clusters <- unique(top3_markers_validation$cluster)

## loop over clusters to save top3 marker genes for each assigned cell identity
for(clst in clusters){
  features <- top3_markers_validation %>% filter(cluster == clst) %>%  .$gene
  p <- FeaturePlot_scCustom(integrated_seurat, features = features, reduction = "umap.scvi", raster = FALSE, ncol = 3)
  p <- ggarrange(p_umap, p, nrow = 1, align = "h", widths = c(0.3, 0.7))
  p <- annotate_figure(p, top = text_grob(paste0(clst, " marker genes"), face = "bold", size = 16))
  ggsave(paste0("Seurat/results/marker_genes/", clst,"_marker_genes.png"), p, width =600, height = 200, unit = "mm", bg = "white")
}

## check whether assigned cluster express markers of their assigned cell dientitdy andw hether they show markers genes as shown in PanglaoDB
p <- DotPlot(integrated_seurat, features = unique(top3_markers_validation$gene)) +   theme(axis.text.x = element_text(angle = 90))
ggsave("Seurat/results/Dotplot_marker_genes_validation.png", p, width = 400, height = 200, unit = "mm", bg = "white")

## save results
saveRDS(integrated_seurat, file = "../../scRNAseq_data/GSE171524/rds/06_GSE171524_seurat_annotated.rds")
write.csv(markers, "Seurat/results/marker_genes/original_markers.csv")
write.csv(markers_validation, "Seurat/results/marker_genes/validated_markers.csv")


p1 <- DimPlot_scCustom(integrated_seurat, group.by = ("Cell_label"),  reduction = "umap.harmony", combine = TRUE, 
                 label = TRUE, repel = TRUE, label.size = 4, raster = FALSE)
p2 <- DimPlot_scCustom(integrated_seurat, group.by = ("harmony_clusters"),  reduction = "umap.harmony", combine = TRUE, 
                 label = TRUE, repel = TRUE, label.size = 4, raster = FALSE)
ggsave("Seurat/results/test.png", ggarrange(p1, p2, nrow = 2), width = 210, height = 297, unit = "mm")

# integrated_seurat$Name <- unlist(str_remove(rownames(integrated_seurat@meta.data), "^[^_]*_"))
# integrated_seurat$Name <-str_replace(integrated_seurat$Name, "\\.", "-")
# 
# metadata_study <- read_delim("../../scRNAseq_data/GSE171524/metadata/GSE171524_lung_metaData.txt")[-1,]
# metadata_study <- metadata_study %>% select(NAME, contains("cell_type"))
# 
# integrated_seurat$cell_type_intermediate <- metadata_study$cell_type_intermediate[match(integrated_seurat@meta.data$Name, metadata_study$NAME)]