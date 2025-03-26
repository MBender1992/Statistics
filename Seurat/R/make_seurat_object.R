#*****************************************
# 1. Follow tutorial here: https://github.com/hbctraining/ and https://www.youtube.com/watch?v=uvyG9yLuNSE

# Have a look at supplementary from the following paper to adjust workflow
# Excerpt taken from "Single-cell RNA sequencing reveals melanoma cell state dependent
# heterogeneity of response to MAPK inhibitors" 
# Lim et al. 2024

#*****************************************
# 2. Install and load packages 

# load packages
library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)

## extract files based on project folder within scRNAseq folder (located within home folder)
project <- "GSE171524"

## define path where files are located
path_to_files <- list.files(path = paste0("~/scRNAseq/", project, "/csv"), full.names = TRUE)

# ## define path where metadata is located
# path_to_metadata <- list.files(path = paste0("~/scRNAseq/", project, "/metadata"), full.names = TRUE)
# 
# ## load metadata
# metadata_orig <- readr::read_delim(file = path_to_metadata, delim = "\t")
# metadata <- metadata_orig %>% select(NAME, biosample_id, donor_id, disease, disease__ontology_label)

## function to generate seurat object for 1 file
load_scData <- function(file){
  sample <- unlist(strsplit(file, "_"))[2]
  countsData <- read.csv(file = file, header = TRUE, row.names = 1)
  seurat_obj <- CreateSeuratObject(counts = countsData, project = sample)
}

## apply load_scData function over list of files
list_seurat <- lapply(path_to_files, load_scData)
names(list_seurat) <- as.data.frame(do.call(rbind, strsplit(path_to_files, "_")))[,2]

## save seurat object
saveRDS(list_seurat, file = paste0("~/scRNAseq/", project, "/rds/", paste0(project, "_seurat_list.rds")))


