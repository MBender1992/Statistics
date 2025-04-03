#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE


## https://www.biostars.org/p/454292/
# I find that the scRNA-seq dim plots will sometimes melt illustrator as well. 
# I found a handy library, ggraster, that will let you rasterize the points but keep the labels and such as vectors. 
# I figured I would pass this along in case you haven't seen it.


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
library(BiocParallel)
library(RColorBrewer)
library(ComplexHeatmap)

#*****************************************
# 3. Load data

## load seurat object
integrated_seurat <- readRDS("../../scRNAseq_data/GSE171524/rds/05b_GSE171524_seurat_annotated.rds")

## apply in the QC filter step (should be changed to another script)
integrated_seurat <- integrated_seurat[!str_detect(rownames(integrated_seurat),  "LINC|AC\\d{3,}|AL\\d{3,}|AP\\d{3,}") ,]

#*****************************************
# 4. Differential expression analysis

## save metadata object
metadata <- integrated_seurat@meta.data

## plot cell type composition
p <- metadata %>% 
  group_by(Sample, Condition, Cell_label)  %>%
  summarize(n_cells = n()) %>% 
  mutate(percentage = n_cells/sum(n_cells)) %>%
  ggplot(aes(Cell_label, percentage, color = Condition)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(0.1)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  ylab("Frequency among all cells") + 
  scale_color_manual(values = c("darkblue", "darkred"))

ggsave("Seurat/results/cell_composition.png", p, height = 100, width = 210, unit = "mm")



## identify top markers for AT1 vs. AT2 and AT2 vs. AT1
topmarkers1 <- FindMarkers(integrated_seurat, ident.1 = "AT1", ident.2 = "AT2", 
                           only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA") 
topmarkers2 <- FindMarkers(integrated_seurat, ident.1 = "AT2", ident.2 = "AT1", 
                           only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
de_markers <-rbind(top_n(topmarkers1, 25, avg_log2FC), top_n(topmarkers2, 25, avg_log2FC))

##
keep_labels <- integrated_seurat$Cell_label %in% c("AT1", "AT2")

## add known markers, filter normalized data object and convert to matrix. Scale matrix
genes <- c(rownames(de_markers), "SFTPB", "SFTPC", "SFTPD", "ETV5", "AGER", "CLIC5")
mat <- GetAssayData(integrated_seurat, layer = "data")[genes, keep_labels] %>% as.matrix()
mat<- t(scale(t(mat)))

## define annotation bar
cluster_anno <- integrated_seurat@meta.data$Cell_label[keep_labels]

scales::viridis_pal()(50)
## make the black color map to 0. the yellow map to highest and the purle map to the lowest
Seurat::PurpleAndYellow()
col_fun <- circlize::colorRamp2(c(-2, 0, 2),   c("#FF00FF", "black", "#FFFF00"))

## plot Heatmap
ht <- Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        clustering_distance_rows = "euclidean", 
        clustering_distance_columns = "euclidean", 
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 4),
        column_title_rot = 90,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 5)

## save Heatmap
png("Seurat/results/AT_differential_expression_Heatmap.png", width = 210, height = 297, unit = "mm", res = 300)
draw(ht)
dev.off()

## violin plots
p1 <- VlnPlot_scCustom(integrated_seurat[, integrated_seurat$Cell_label == "AT1"], features = "CAV1", colors_use = c("darkblue", "darkred"), 
                       group.by = "Cell_label", alpha = 0.4, raster = F, split.by = "Condition")
p2 <- VlnPlot_scCustom(integrated_seurat[, integrated_seurat$Cell_label == "AT2"], features = "ETV5", colors_use = c("darkblue", "darkred"), 
                       group.by = "Cell_label", alpha = 0.4, raster = F, split.by = "Condition")
p3 <- ggarrange(p1, p2, nrow = 1, align = "hv")
ggsave("Seurat/results/marker_expression_by_condition.png", p3, height = 100, width = 210, unit = "mm")

## differential expression with Pseudobulk Analysis and DESeq2

## Go analysis similar to bulk RNASeq

## Score gene signature
## https://satijalab.org/seurat/reference/addmodulescore
## for plotting set na_cutoff = NA in featurePlot_scCustom
## --> can also be used for cell identification (higher score corersponds to higher likelihood of specific cell identity)
