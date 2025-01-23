## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## load packages
library(TCGAbiolinks)

## access gene expression data
query.exp.skcm <- GDCquery(project = "TCGA-SKCM", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
                           workflow.type = "STAR - Counts")
GDCdownload(query.exp.skcm)

View(getResults(query.exp.skcm))



## access clinical data from melanoma samples
skcm_clin <- GDCquery_clinic(project = "TCGA-SKCM", type = "Clinical")

## access clinical data from melanoma samples in xml format (only for 2 patients to speed up the processs)
query_clin <- GDCquery(project = "TCGA-SKCM", data.category = "Clinical", barcode = c("TCGA-EE-A29R", "TCGA-D3-A2JK"), data.format = "BCR")
GDCdownload(query_clin)
clinical.patient <- GDCprepare_clinic(query_clin, "patient")
clinical.drug <-  GDCprepare_clinic(query_clin, "drug")
clinical.followup <-  GDCprepare_clinic(query_clin, "follow_up")
clinical.admin <-  GDCprepare_clinic(query_clin, "admin")
clinical.nte <-  GDCprepare_clinic(query_clin, "new_tumor_event")
clinical.stage <-  GDCprepare_clinic(query_clin, "stage_event")
