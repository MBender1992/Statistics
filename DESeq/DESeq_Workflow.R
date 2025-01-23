#*****************************************
# 1. Find a data set on this page: http://www.bioconductor.org/packages/release/data/experiment/

#*****************************************
# 2. Install and load the data set
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install(c("DEGreport", "tximport"))
# install.packages("pheatmap")

# load packages
library(tidyverse)
library(DESeq2)
library(airway)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(apeglm)
library(EnhancedVolcano)
library(pheatmap)
library(ComplexHeatmap)
library(ReportingTools)
library(PoiClaClu)
library(msigdbr)
library(clusterProfiler)

#*****************************************
# 3. Convert the data to a DESeq dataset object with an appropriate design formula
# installing/loading the package:

data(airway)
dds <- DESeqDataSet(airway, design = as.formula("~cell+dex"))
assay(dds)[rownames(dds) == "ENSG00000002586",]


# plot heatmap of Poisson distances between samples
# IMPORTANT: use Poisson distance for raw (non-normalized) count data
# IMPORTATN: use Euclidean distance for data normalized by regularized-logarithm transformation (rlog) or variance stablization transfromation (vst)
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd)
rownames(samplePoisDistMatrix) <- paste(dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

#*****************************************
# 3a. Run EDA

smallestGroupSize <- 4 # A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 4 samples in each group.
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize # Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. The count of 10 is a reasonable choice for bulk RNA-seq.
dds <- dds[keep,]
nrow(dds)

# Which transformation to choose? The VST is much faster to compute and is less sensitive to high count outliers
# than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST 
# when there is a wide range of sequencing depth across samples (an order of magnitude difference). 
# We therefore recommend the VST for medium-to-large datasets (n > 30). You can perform both transformations and 
# compare the meanSdPlot or PCA plots generated, as described below.

# variance stabilized transformation
vsd <- vst(dds, blind = FALSE)

# rlog transformation (regularized log transformation)
rld <- rlog(dds, blind = FALSE)

# To show the effect of the transformation, in the figure below we plot the first sample against the second,
# first simply using the log2 function (after adding 1, to avoid taking the log of zero), and then using the 
# VST and rlog-transformed values. For the log2 approach, we need to first estimate size factors to account for
# sequencing depth, and then specify normalized=TRUE. Sequencing depth correction is done automatically for the 
# vst and rlog.

# estimate size factors
dds <- estimateSizeFactors(dds)

# build data set
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

# assign column names
colnames(df)[1:2] <- c("x", "y")  

# transform to factor
lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

# plot values
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

###
# plot PCA
pcaData <- plotPCA(vsd, intgroup = c(factor1, factor2), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# plot data with ggplot
ggplot(pcaData, aes_string(x = "PC1", y = "PC2", color = factor1, shape = factor2)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data") +
  theme_bw() + 
  ggsci::scale_color_npg()

#*****************************************
# 4. Run DESeq
dds <- DESeq(dds)

#*****************************************
# 5. Get the results and plot Volcano plot without filter within results function
res_volcano <- results(dds)

### annotate genes
# for all genes
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(res_volcano), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")

# adding ENSEMBL gene ID as a column in significant differentially expression gene table.
dat_volcano <- cbind(ENSEMBL = rownames(res_volcano), res_volcano)
dat_volcano <- left_join(as.data.frame(dat_volcano), anno)

# define thresholds
lfcThres <- log2(1.5) ## always specify this filter criterion to take the desired difference into account in significance testing
pThres   <- 0.05

# plot Volcano
# this only serves as diagnostic and representation tool. the actual fold change thresholds and pvalues change later due to the way the results function works
EnhancedVolcano(as.data.frame(dat_volcano), lab = dat_volcano$SYMBOL, 
                x = 'log2FoldChange', y = 'padj',
                xlim = c(-8, 8), title = 'Treated vs untreated',
                pCutoff = pThres, FCcutoff = lfcThres, pointSize = 1.5, labSize = 3)


#*****************************************
# 5a. Get the results and apply specified filters

### run the desired results with the specified p and lfc threshold 
# apply lfcThreshold and alpha of 0.05
res <- results(dds, lfcThreshold = lfcThres, alpha = pThres)
summary(res)
# This procedure is very strict. If too few statistically significant genes are present it might be an option to test against lfc !=0 and later on filter the lfc.
# However this approach is preferred.

# res <- results(dds, contrast=c("dex","trt","untrt"))
# summary(res)
# 
# # adjust pvalue cutoff 
# res.05 <- results(dds, alpha = 0.05)
# table(res.05$padj < 0.05)
# 
# # adjust FC cutoff, this changes the Null hypothesis as now it is specified whether there are differences >0.58 --> different pvals than with a subsequent filter
# resLFC1 <- results(dds, lfcThreshold=0.58)
# summary(resLFC1)

#*****************************************
# 5b. Explore results and run some diagnostics

### plot top gene
# extract top gene
topGene <- rownames(res)[which.min(res$padj)]
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
# plot data 
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

# An MA-plot (Dudoit et al. 2002) provides a useful overview for the distribution of the estimated
# coefficients in the model, e.g. the comparisons of interest, across all genes. On the y-axis, the “M”
# stands for “minus” – subtraction of log values is equivalent to the log of the ratio – and on the x-axis,
# the “A” stands for “average”. You may hear this plot also referred to as a mean-difference plot, or a
# Bland-Altman plot.
# 
# Before making the MA-plot, we use the lfcShrink function to shrink the log2 fold changes for the comparison
# of dex treated vs untreated samples. There are three types of shrinkage estimators in DESeq2, which are covered
# in the DESeq2 vignette. Here we specify the apeglm method for shrinking coefficients, which is good for shrinking
# the noisy LFC estimates while giving low bias LFC estimates for true large differences (Zhu, Ibrahim, and Love 2018).
# To use apeglm we specify a coefficient from the model to shrink, either by name or number as the coefficient appears
# in resultsNames(dds).

### plot MA plot
resultsNames(dds)

# plot with original lfc
plotMA(res, ylim = c(-5, 5))

# shrink estimates and exlpore differences in MA plots
resShrunk <- lfcShrink(dds, coef="dex_untrt_vs_trt", type="apeglm")
plotMA(resShrunk, ylim = c(-5, 5))

###
# WARNING: DO NOT ALWAS SHRINK EVERYTHING BY DEFAULT. Look at MA plot and Volcano plot to see if it makes sense
### 

#*****************************************
# 6. Convert gene symbols
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(res), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")
#remove duplicated ENSEMBL ID entries
anno <- anno %>% filter(!duplicated(ENSEMBL))

# adding ENSEMBL gene ID as a column in significant differentially expression gene table.
allRes <- cbind(ENSEMBL = rownames(res), res)
allRes <- if(dim(anno)[1] == dim(allRes)[1]){
  left_join(as.data.frame(allRes), anno)
} else {
  stop("Dimensions of annotation object and result object are different.")
}

# subset significant results
sigRes <- subset(allRes, padj < pThres)

#*****************************************
# 7. Plot heatmap

# If using either of Euclidean distance or Pearson correlation, your data should follow a Gaussian /
#   normal (parametric) distribution. So, if coming from a microarray, anything from RMA normalisation is fine,
# whereas, if coming from RNA-seq, any data deriving from a transformed normalised count metric should be fine,
# such as variance-stabilised, regularised log, or log CPM expression levels.
# 
# If you are performing clustering on non-normal data, like 'normalised' [non-transformed] RNA-seq counts,
# FPKM expression units, etc., then use Spearman correlation (non-parametric).

# Heatmap of the top significantly enriched genes
nGenes <- 100
topSignif <- head(order(allRes$padj, decreasing = FALSE),nGenes)
labels <- allRes[topSignif,]$SYMBOL
mat  <- assay(vsd)[topSignif, ]
annoHeatmap <- as.data.frame(colData(vsd)[, c("cell","dex")])
pheatmap(mat, annotation_col = annoHeatmap, labels_row = labels, scale = "row", clustering_distance_rows = "correlation", 
         main= paste("Top ",nGenes," differentially expressed genes", sep =""))

# Heatmap of all DE genes
mat <- assay(vsd)[sigRes$ENSEMBL,]
annoHeatmap <- as.data.frame(colData(rld)[, c("cell","dex")])
pheatmap(mat, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation",
         annotation_col = annoHeatmap, main="Differentially Expressed genes")


## add example with ComplexHeatmap
mat  <- assay(vsd)[topSignif, ]
matScaled <- mat %>% t() %>%  scale() %>% t()

# define colors for groups and colorbar
colFun <- colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100)
colorbar <- HeatmapAnnotation(
  df =anno,
  col = list(cell = c("N61311" = "#0073C2FF","N052611" = "#EFC000FF","N080611"="#868686FF" ,"N061011"="#CD534CFF"),
             dex  = c("untrt" = colFun[1], trt = colFun[100])),
  annotation_legend_param = list(cell = list(nrow=1), 
                                 dex  = list(nrow=1))
  )

# draw Heatmap
Ht <- Heatmap(matScaled,
              col= colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
              top_annotation = colorbar,
              row_labels = labels,
              clustering_method_row = "complete",
              clustering_method_columns = "complete",
              clustering_distance_row = "pearson",
              clustering_distance_column = "euclidean",
              rect_gp = gpar(col = "grey60",lty = 1, lwd = 1),
              row_names_gp = gpar(fontsize = 10),
              show_column_names = FALSE,
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(
                title = "row Z-score",
                at = seq(-2,2,by=1),
                color_bar="continuous",
                title_position ="topcenter",
                legend_direction = "horizontal",
                legend_width = unit(4, "cm")
              ))

# draw Heatmap
draw(Ht, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")


#*****************************************
# 8. Overrepresentation analysis

### group overexpressed genes
over_expressed_genes <- allRes %>%
  filter(padj < pThres & log2FoldChange > 0) %>%
  pull(SYMBOL)

### group underexpressed genes
under_expressed_genes <- allRes %>%
  filter(padj < pThres & log2FoldChange < 0) %>%
  pull(SYMBOL)

# Get the gene sets and wrangle
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
# gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
gene_sets <- gene_sets %>%
  select(gs_name, gene_symbol)

# Run over-representation analysis
upReg <- enricher(gene = over_expressed_genes,
                 TERM2GENE = gene_sets)
# Plot results with clusterProfiler
dotplot(upReg)

# Run under-representation analysis
downReg <- enricher(gene = under_expressed_genes,
                 TERM2GENE = gene_sets)
# Plot results with clusterProfiler
dotplot(downReg)

#*****************************************
# 8a. GSEA

# Deal with inf
allResGSEA <- allRes %>%
  mutate(padj = case_when(padj == 0 ~ .Machine$double.xmin,
                          TRUE ~ padj)) %>%
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange))

# Remove NAs and order by GSEA
allResGSEA <- allResGSEA %>%
  filter(! is.na(gsea_metric)) %>%
  arrange(desc(gsea_metric))

# Get the ranked GSEA vector
ranks <- allResGSEA %>%
  select(SYMBOL, gsea_metric) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  deframe()

# Run GSEA
gseares <- GSEA(geneList = ranks, 
                TERM2GENE = gene_sets)
gsearesdf <- as.data.frame(gseares)

# Make summary plots
dotplot(gseares)

# Example of GSEA plot
gseaplot(gseares, geneSetID = "GOBP_RNA_SPLICING",
         title = "GOBP_RNA_SPLICING")

#*****************************************
# 9. accounting for batch effects

#*****************************************
# 10. Export data
resOrdered <- allRes[order(allRes$padj),]
resOrderedDF <- as.data.frame(resOrdered)[1:n_genes,]
write.csv(resOrderedDF, file = paste("DESeq_results_top_",nGenes ,".csv", sep = ""))

# write html report
htmlRep <- HTMLReport(shortName="report", title="DESeq2_results",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)





