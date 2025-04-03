#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ & https://www.youtube.com/watch?v=uvyG9yLuNSE

# Have a look at supplementary from the following paper to adjust workflow
# Excerpt taken from "Single-cell RNA sequencing reveals melanoma cell state dependent
# heterogeneity of response to MAPK inhibitors" 
# Lim et al. 2024

#*****************************************
# 2. Install and load packages 

# load packages
library(tidyverse)
library(AnnotationHub)
library(ensembldb)

#*****************************************
# 3. Extract annotations

### NOT RUN (only re run if upated annotation should be extracted )
# # Connect to AnnotationHub
# ah <- AnnotationHub()
# 
# # Access the Ensembl database for organism
# ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE)
# 
# # Check versions of databases available
# ahDb %>%  mcols()
# 
# # Acquire the latest annotation files
# id <- ahDb %>%
#   mcols() %>%
#   rownames() %>%
#   tail(n = 1)
# 
# # Download the appropriate Ensembldb database
# edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, return.type = "data.frame")

# # Explore biotypes to analyze mitochondrial gene expression
# annotations$gene_biotype %>%   factor() %>%  levels()
# 
# # Extract IDs for mitochondrial genes
# mt <- annotations %>% dplyr::filter(seq_name == "MT") %>% dplyr::pull(gene_id)

saveRDS(annotations, file = "Seurat/data/gene_annotations.rds")


