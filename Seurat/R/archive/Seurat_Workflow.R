#*****************************************
# 1. Follow tutorial here: https://github.com/quadbio/scRNAseq_analysis_vignette/blob/master/Tutorial.pdf

# Have a look at supplementary from the following paper to adjust workflow
# Excerpt taken from "Single-cell RNA sequencing reveals melanoma cell statedependent
# heterogeneity of response to MAPK inhibitors" 
# Lim et al. 2024

#*****************************************
# 2. Install and load the data set
# install.packages("Seurat")

# load packages
library(tidyverse)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggpubr)

## create seurat object
counts <- readMM("Seurat/data/archive/DS1/matrix.mtx.gz")
barcodes <- read.table("Seurat/data/archive/DS1/barcodes.tsv.gz", stringsAsFactors=F)[,1]
features <- read.csv("Seurat/data/archive/DS1/features.tsv.gz", stringsAsFactors=F, sep="\t",
                     header=F)
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes
seurat <- CreateSeuratObject(counts, project="DS1")

#*****************************************
# 3. Quality control

## calculate percentage of mitochondrial transcripts (which need to be filtered later)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

## plot QC metrics 
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## calculate correlation between number of counts and detected transcripts as well as mtRNA percentage
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Due to the correlation of gene number and transcript number, we only need to set a cutoff to either one of these
# metrics, combined with an upper threshold of mitochondrial transcript percentage, for the QC. For instance, for this
# data set, a detected gene number between 500 and 5000, and a mitochondrial transcript percentage lower than
# 5% would be quite reasonable, but it is fine to use different thresholds
seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &  percent.mt < 5)

# It is worth to mention that sometimes more QC may need to be applied. One potential issue is the presence of
# doublets. As the amount of captured RNA varies a lot from cell to cell, doublets don't always show a higher number
# of detected genes or transcripts. There are several tools available now, which are designed to predict whether a
# 'cell' is indeed a singlet or actually a doublet/multiplet. DoubletFinder, for instance, predicts doublets by first
# constructing artificial doublets by randomly averaging cells in the data, and then for each cell testing whether it is
# more similar to the artificially doublets or not. This helps with the decision whether a cell is likely a doublet or not.
# Similarly, mitochondrial transcript percentage may not be sufficient to filter out stressed or unhealthy cells.
# Sometimes one would needs to do extra filtering, e.g. based on the machine learning based prediction.

#*****************************************
# 4. Data normalization

# Similar to bulk RNA-seq, the amount of captured RNA is different from cell to cell, and one should therefore not
# directly compare the number of captured transcripts for each gene between cells. A normalization step, aiming to
# make gene expression levels between different cells comparable, is therefore necessary. The most commonly used
# normalization in scRNA-seq data analysis is very similar to the concept of TPM (Transcripts Per Million reads) - one
# normalizes the feature expression measurements for each cell to the total expression, and then multiplies this by a
# scale factor (10000 by default). At the end, the resulting expression levels are log-transformed so that the
# expression values better fit a normal distribution. It is worth to mention that before doing the log-transformation,
# one pseudocount is added to every value so that genes with zero transcripts detected in a cell still present values of
# zero after log-transform.
seurat <- NormalizeData(seurat)

#*****************************************
# 5. Feature selection for heterogeneity analysis

# By default, Seurat calculates the standardized variance of each gene across cells, and picks the top 2000 ones as
# the highly variable features. One can change the number of highly variable features easily by giving the nfeatures
# option (here the top 3000 genes are used).
# There is no good criteria to determine how many highly variable features to use. Sometimes one needs to go
# through some iterations to pick the number that gives the most clear and interpretable result. Most often, a value
# between 2000 to 5000 is OK and using a different value doesn't affect the results too much.
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)

top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)

#*****************************************
# 6. Scale data

# Since different genes have different base expression levels and distributions, the contribution of each gene to the
# analysis is different if no data transformation is performed. This is something we do not want as we don't want our
# analysis to only depend on genes that are highly expressed. Therefore a scaling is applied to the data using the
# selected features, just like one usually does in any data science field.

seurat <- ScaleData(seurat)

### IMPORTANT NOTE ###
# Variables which are commonly considered to be regressed out include the number of detected genes/transcripts
# (nFeature_RNA / nCount_RNA), mitochondrial transcript percentage (percent.mt), and cell cycle related variables
# (see below). What it tries to do is to first fit a linear regression model, using the normalized expression level of a
# gene as the dependent variable, and the variables to be regressed out as the independent variables. Residuals of
# the linear model are then taken as the signals with the linear effect of the considered variables removed. I should
# note that this process of regressing out variables dramatically slows down the whole process, and it is not clear that
# the result will be satisfactory as the unwanted variation may be far from linear. Therefore, a common suggestion is
# not to perform any regress-out in the first iteration of data exploration, but first check the result, and if any
# unwanted source of variation dominates the cellular heterogeneity, try to regress out the respective variable and
# see whether the result improves.
#######################

### OPTIONAL ###
# Use SCTransform to better approximate normalization and remove zero-inflation
################

#*****************************************
# 7. Dimensionality reduction using PCA

## run PCA
seurat <- RunPCA(seurat, npcs = 50)
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))
png("Seurat/results/PCA_top20_Heatmap.png", width = 10, height = 8, unit = "in", res = 300)
PCHeatmap(seurat, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
dev.off()

# In this example, we would use the top-20 PCs for the following analysis. Again, it is absolutely fine to use more or
# fewer PCs, and in practice this sometimes needs some iteration to make the final decision. Meanwhile, for most of
# the data, a PC number ranging from 10 to 50 would be reasonable and in many cases it won't affect the conclusion
# very much (but sometimes it will, so still be cautious).

#*****************************************
# 8. Non-linear dimension reduction for visualization

seurat <- RunTSNE(seurat, dims = 1:20)
seurat <- RunUMAP(seurat, dims = 1:20)

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2

plot1 <- FeaturePlot(seurat,
                    c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                    ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat,
                     c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "umap")
plot1 / plot2

#*****************************************
# 9. Cluster the cells

# Doing feature plot of markers is usually a good way to start with when exploring scRNA-seq data. However, to more
# comprehensively understand the underlying heterogeneity in the data, it is necessary to identify cell groups with an
# unbiased manner. This is what clustering does. In principle, one can apply any clustering methods, including those
# widely used in bulk RNA-seq data analysis such as hierarchical clustering and k-means, to the scRNA-seq data.
# However, in practice, this is very difficult, as the sample size in scRNA-seq data is too much larger (one 10x
#                                                                                                        experiment usually gives several thousands of cells). It would be extremely slow to use these methods. In addition,
# due to the intrinsic sparseness of scRNA-seq data, even if data is denoised by dimension reduction like PCA,
# differences between different cells are not as well quantitative as those of bulk RNA-seq data. Therefore, the more
# commonly used clustering methods in scRNA-seq data analysis is graph-based community identification algorithm.
# Here, graph is the mathematical concept, where there is a set of objects, and some pairs of these objects are
# related with each other; or in a simplified way, a network of something, and here, a network of cells.
# First of all, a k-nearest neighbor network of cells is generated. Every cells is firstly connected to cells with the
# shortest distances, based on their corresponding PC values. Only cell pairs which are neighbors of each other are
# considered as connected. Proportion of shared neighbors between every cell pairs is then calculated and used to
# describe the strength of the connection between two cells. Weak connections are trimmed. This gives the resulted
# Shared Nearest Neighbor (SNN) network. In practice, this is very simple in Seurat.

## cluster data
seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 1)

## plot clusters
plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2

#*****************************************
# 10. Annotate cell clusters

## test which clusters express known markers
ct_markers <- c("MKI67","NES","DCX","FOXG1", # G2M, NPC, neuron, telencephalon
                "DLX2","DLX5","ISL1","SIX3","NKX2.1","SOX6","NR2F2", # ventral telencephalon related
                "EMX1","PAX6","GLI3","EOMES","NEUROD6", # dorsal telencephalon  related
                "RSPO3","OTX2","LHX9","TFAP2A","RELN","HOXB2","HOXB5") # non-telencephalon related
DoHeatmap(seurat, features = ct_markers) + NoLegend()

# Next, in order to do annotation in a more unbiased way, we should firstly identify cluster markers for each of the cell
# cluster identified. In Seurat, this can be done using the FindAllMarkers function. What it does is for cell cluster, to
# do differential expression analysis (with Wilcoxon's rank sum test) between cells in the cluster and cells in other
# clusters

## calculate marker genes for each cluster
cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% print(n="all")

## visualize markers
top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


png("Seurat/results/Marker_Heatmap.png", width = 15, height = 12, unit = "in", res = 300)
DoHeatmap(seurat, features = top10_cl_markers$gene) + NoLegend()
dev.off()

plot3 <- FeaturePlot(seurat, c("NEUROD2","NEUROD6"), ncol = 1)
plot4 <- VlnPlot(seurat, features = c("NEUROD2","NEUROD6"), pt.size = 0)
(plot4 + plot3) + plot_layout(ncol = 1, widths = c(1,2))

## arrange plots
p1 <- ggarrange(plot3, plot4)
p2 <- ggarrange(plot1, plot2, nrow = 1)
ggarrange(p1, p2, nrow = 2)

## replace cluster names by cell types
new_ident <- setNames(c("Dorsal telen. NPC", #0 0 
                        "Dien. and midbrain excitatory neuron", #3 1
                        "Dorsal telen. neuron", #2 2
                        "Midbrain-hindbrain boundary neuron", #1 3
                        "Dien. and midbrain IP and excitatory early neuron",# 8 4
                        "MGE-like neuron",#4 5
                        "G2M dorsal telen. NPC",  #5 6
                        "Dien. and midbrain NPC",  #7 7 
                        "Dorsal telen. IP", #6 8 
                        "Dien. and midbrain IP and early inhibitory neuron", #12 9
                        "G2M Dien. and midbrain NPC", #9 10
                        "G2M dorsal telen. NPC", #10 11
                        "Dien. and midbrain inhibitory neuron", #11 12
                        "Ventral telen. neuron", #13 13
                        "Unknown 1"), 
                      levels(seurat))
seurat <- RenameIdents(seurat, new_ident)
DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()

#*****************************************
# 11. Pseudotemporal cell ordering

# As you may have noticed, there are two clusters of dorsal telencephalic NPCs separated from the third cluster
# because they are at different phase of cell cycle. Since we are interested in the general molecular changes during
# differentiation and maturation, the cell cycle changes may strongly confound the analysis. We can try to reduce cell
# cycle effect by excluding cell cycle related genes from the identified highly variable gene list.
seurat_dorsal <- subset(seurat, subset = RNA_snn_res.1 %in% c(0,2,5,6,10))
seurat_dorsal <- FindVariableFeatures(seurat_dorsal, nfeatures = 2000)
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), unlist(cc.genes))

seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)

## remove cell cycle effect
seurat_dorsal <- CellCycleScoring(seurat_dorsal,
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes,
                                  set.ident = TRUE)
seurat_dorsal <- ScaleData(seurat_dorsal, vars.to.regress = c("S.Score",  "G2M.Score"))

seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)

