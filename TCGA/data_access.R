## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

# https://www.biostars.org/p/181470/ # Ideas to calculate survival | days_to_sample_procrurement or days to collection? 
# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=sample # GDC dictionary

## load packages
library(TCGAbiolinks)
library(tidyverse)

## specify directory where to store data downloaded from GDC
dir <- "./TCGA/data"

##*********************************************************************************************************
## Clinical data

## access clinical data from melanoma samples
skcm_clin <- GDCquery_clinic(project = "TCGA-SKCM", type = "Clinical")

## access clinical data from melanoma samples in xml format (only for 2 patients to speed up the processs)
query_clin <- GDCquery(project = "TCGA-SKCM", data.category = "Clinical", data.format = "BCR XML")
# GDCdownload(query_clin, directory = dir, files.per.chunk = 50)
clinical.patient <- GDCprepare_clinic(query_clin, "patient", directory = dir)
clinical.drug <-  GDCprepare_clinic(query_clin, "drug", directory = dir)
clinical.followup <-  GDCprepare_clinic(query_clin, "follow_up", directory = dir)
clinical.admin <-  GDCprepare_clinic(query_clin, "admin", directory = dir)
clinical.nte <-  GDCprepare_clinic(query_clin, "new_tumor_event", directory = dir)
clinical.stage <-  GDCprepare_clinic(query_clin, "stage_event", directory = dir)

## restructure data to be able to calculate survival for different timepoints (e.g. corresponding to time of sample collection)


##*********************************************************************************************************
## sample data

## download data for sample preparation
query_specimen <- GDCquery(project = "TCGA-SKCM", data.category = "Biospecimen", data.format = "BCR XML")
# GDCdownload(query_specimen, directory = dir, files.per.chunk = 1)

##*********************************************************************************************************
## mRNA expression data

## access gene expression data
query_exp_skcm <- GDCquery(project = "TCGA-SKCM", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
#                            workflow.type = "STAR - Counts")
# GDCdownload(query_exp_skcm, directory = dir, files.per.chunk = 1)
skcm_rna <- GDCprepare(query_exp_skcm, save = FALSE, summarizedExperiment = TRUE, directory = dir)
meta_data <- colData(skcm_rna)

##*********************************************************************************************************
## miRNA expression data





data.frame(OS = ifelse(test$vital_status == "Alive", test$days_to_last_follow_up, test$days_to_death), OSCENS = ifelse(test$vital_status == "Alive", 0, 1))


View(getResults(query_clin))


dat_surv <- meta_data %>% as.data.frame() %>% select(patient, sample, days_to_collection, days_to_last_follow_up, days_to_death, vital_status)
dat_drug <- clinical.drug %>% select(bcr_patient_barcode, days_to_drug_therapy_start, days_to_drug_therapy_end) %>% left_join(dat_surv, by = c("bcr_patient_barcode" = "patient"))

clinical.patient %>% .$days_to_submitted_specimen_dx
meta_data %>% as.data.frame() %>% filter(patient == "TCGA-D3-A1Q1")

