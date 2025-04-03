#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

## clear workspace and memory
rm(list = ls())
gc()

#*****************************************
# 2. Install and load packages 

## load packages
library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(scales)
library(cowplot)
library(RCurl)
library(ggprism)
library(ggsci)
library(ggpubr)
library(ggridges)
library(scCustomize)

#*****************************************
# 3. Create custom functions

## adjusted prism theme
theme_prism2 <- function(base_size = 9, base_family = "", ...){
  list(theme_prism(...),  
       theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1), 
             plot.title = element_text(hjust = 0.5, face = "bold")),
       guides(x = "prism_offset", y = "prism_offset"))
}


## calculate qc thresholds for plotting
calc_thresholds <- function(data, metric, nmad = 5, log.transform = FALSE){
  if(log.transform == TRUE){
    lower <- median(log(data[[metric]])) - nmad*mad(log(data[[metric]]))
    upper <- median(log(data[[metric]])) + nmad*mad(log(data[[metric]]))
  } else {
    lower <- median(data[[metric]] - nmad*mad(data[[metric]]))
    upper <- median(data[[metric]] + nmad*mad(data[[metric]]))
  }

  c(lower, upper)
}

## wrapper around RidgePlot
custom_RidgePlot <- function(seurat.obj, metric, upper.xlim = NULL){
  metadata <- seurat.obj@meta.data
  p <- RidgePlot(seurat.obj, metric, cols = colors) +
    ggtitle(paste0(metric," per cell")) +
    theme_prism2() +
    theme(legend.position = "none",
          axis.title.y = element_blank()) +
    geom_vline(xintercept = median(metadata[[metric]]), size = 0.9, lty = 1, color = "darkred") 
  if(!is.null(upper.xlim)) p + xlim(0,upper.xlim) else p
}

#*****************************************
# 4. Load data

## load seurat object
merged_seurat <- readRDS("../../scRNAseq_data/GSE171524/rds/02_GSE171524_merged_seurat.rds")
metadata <- merged_seurat@meta.data

#*****************************************
# 5 Visualization of low quality cells and noise

## define colors from scCustomize palette
colors <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")

## Visualize the number of cell counts per sample
p1 <- metadata %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) + 
  geom_bar() +
  theme_prism2() +
  theme(legend.position = "right", axis.title.x = element_blank()) +
  scale_fill_manual(values = colors) +
  xlab("Identity") +
  ylab("Counts") +
  ggtitle("Number of cells")

## Visualize metrics
p2 <- custom_RidgePlot(seurat.obj = merged_seurat, metric = "nUMI", upper.xlim = 10000)
p3 <- custom_RidgePlot(seurat.obj = merged_seurat,metric = "nGene", upper.xlim = 5000)
p4 <- custom_RidgePlot(seurat.obj = merged_seurat,metric = "log10GenesPerUMI")
p5 <- custom_RidgePlot(seurat.obj = merged_seurat,metric = "MitoRatio")

## arrange plots
p6 <- ggarrange(p2,p3,p4, p5, ncol = 2, nrow = 2, labels = LETTERS[2:5])
p7 <- ggarrange(p1, p6, ncol = 1, nrow = 2, labels = c(LETTERS[1], ""), heights = c(0.37, 0.63))

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p8 <- ggplot(metadata) +
  geom_point(aes(x=nUMI,y=nGene,fill=MitoRatio > 0.05),shape=21,alpha=0.4) + 
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text())  +
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~orig.ident, ncol = 9) +
  geom_vline(xintercept = 500,color="red",linetype="dotted")+
  geom_hline(yintercept=250,color="red", linetype="dotted")

## violin plot of different metrics 
p9a <- VlnPlot_scCustom(merged_seurat, features = "nUMI", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9b <- VlnPlot_scCustom(merged_seurat, features = "nGene", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9c <- VlnPlot_scCustom(merged_seurat, features = "log10GenesPerUMI", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9d <- VlnPlot_scCustom(merged_seurat, features = "MitoRatio", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9  <- ggarrange(p9a, p9b, p9c, p9d, ncol = 2, nrow = 2, labels = LETTERS[1:4])

## save plots
ggsave("Seurat/results/QC_RidgePlot.png", p7, width = 240, height =350, unit = "mm")
ggsave("Seurat/results/combined_qc_scatter.png", p8, width = 300, height =250, unit = "mm")
ggsave("Seurat/results/QC_ViolinPlot.png", p9, width = 300, height =250, unit = "mm")

#*****************************************
#* 6 Application of filters and reexamination of quality criteria

## calculate boolean to specify which cells to keep (the number of MADs is taken from Germain et. al 2020 who suggest 5 MADs for all QC metrics
## except mtRNA for which they suggest 3 MADs)
nUMI_keep  <- !isOutlier(metadata$nUMI, nmads=5, log=TRUE) 
nGene_keep <- !isOutlier(metadata$nGene, nmads=5, log=TRUE)
log10GenesPerUMI_keep <- !isOutlier(metadata$log10GenesPerUMI, nmads=5) 
mito_keep  <- metadata$MitoRatio < 0.1
# mito.keep  <- !(isOutlier(metadata$mitoRatio, nmads=3, type="higher")) 

## calculate thresholds for plotting 
nUMI_thresholds <- calc_thresholds(metadata, "nUMI", nmad = 5, log.transform = TRUE) ## good default would be higher than 500 
nGene_thresholds <- calc_thresholds(metadata, "nGene", nmad = 5, log.transform = TRUE)  ## good default is higher than 200-250
log10GenesPerUMI_thresholds <- calc_thresholds(metadata, "log10GenesPerUMI", nmad = 5, log.transform = FALSE) ## good default is higher than 0.8
MitoRatio_thresholds <- calc_thresholds(metadata, "MitoRatio", nmad = 3, log.transform = FALSE) ## good default is lower than 20% or 10%

## apply all fitlers simultaneously
qc_pass <- nUMI_keep & nGene_keep & log10GenesPerUMI_keep & mito_keep

## filter seurat object based on median absolute deviations (as suggested by Sanbomics and Germain et. al 2020)
filtered_seurat <- merged_seurat[,which(qc_pass)]

## remove genes with zero counts
## extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

## Output logical matrix specifying for each gene whether or not there are more than 0 counts per cell
nonzero <- counts > 0

## sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

## Only keep those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

## Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

## extract filtered metadata
metadata_filtered <- filtered_seurat@meta.data

#*****************************************
# 7. Visualization after filtering

## Visualize the number of cell counts per sample
p1 <- metadata_filtered %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) + 
  geom_bar() +
  theme_prism2() +
  theme(legend.position = "right", axis.title.x = element_blank()) +
  scale_fill_manual(values = colors) +
  xlab("Identity") +
  ylab("Counts") +
  ggtitle("Number of cells")

## Visualize metrics
p2 <- custom_RidgePlot(seurat.obj = filtered_seurat, metric = "nUMI", upper.xlim = 10000)
p3 <- custom_RidgePlot(seurat.obj = filtered_seurat,metric = "nGene", upper.xlim = 5000)
p4 <- custom_RidgePlot(seurat.obj = filtered_seurat,metric = "log10GenesPerUMI")
p5 <- custom_RidgePlot(seurat.obj = filtered_seurat,metric = "MitoRatio")

## arrange plots
p6 <- ggarrange(p2,p3,p4, p5, ncol = 2, nrow = 2, labels = LETTERS[2:5])
p7 <- ggarrange(p1, p6, ncol = 1, nrow = 2, labels = c(LETTERS[1], ""), heights = c(0.37, 0.63))

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p8 <- ggplot(metadata_filtered) +
  geom_point(aes(x=nUMI,y=nGene,fill=MitoRatio > 0.05),shape=21,alpha=0.4) + 
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text())  +
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~orig.ident, ncol = 9) +
  geom_vline(xintercept = 500,color="red",linetype="dotted")+
  geom_hline(yintercept=250,color="red", linetype="dotted")

## violin plot of different metrics 
p9a <- VlnPlot_scCustom(filtered_seurat, features = "nUMI", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9b <- VlnPlot_scCustom(filtered_seurat, features = "nGene", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9c <- VlnPlot_scCustom(filtered_seurat, features = "log10GenesPerUMI", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9d <- VlnPlot_scCustom(filtered_seurat, features = "MitoRatio", alpha = 0.4, raster = F, plot_boxplot = TRUE) & NoLegend()
p9  <- ggarrange(p9a, p9b, p9c, p9d, ncol = 2, nrow = 2, labels = LETTERS[1:4])

## save plots
ggsave("Seurat/results/QC_RidgePlot_filtered.png", p7, width = 240, height =350, unit = "mm")
ggsave("Seurat/results/combined_qc_scatter_filtered.png", p8, width = 300, height =250, unit = "mm")
ggsave("Seurat/results/QC_ViolinPlot_filtered.png", p9, width = 300, height =250, unit = "mm")

## Create .rds object to load at any time
saveRDS(filtered_seurat, file = "../../scRNAseq_data/GSE171524/rds/03_GSE171524_seurat_filtered.rds")
