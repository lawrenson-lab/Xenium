library(Seurat)
library(tidyverse)
library(SpatialExperiment)
library(scuttle)
slide<- commandArgs(trailingOnly=TRUE)
#slide<-"/media/Xenium_On_NAS/20250428__202040__042825_Batch4Endometrioma_Rerun/output-XETG00426__0034237__EDV013__20250428__202100"

xenium.obj <- LoadXenium(slide,fov="fov",
                         molecule.coordinates = F,
                         assay = "counts",
                         cell.centroids = T)
slide<-str_remove(slide,".+output-")
slide<-paste(unlist(str_split(slide,"__"))[3:2],collapse = "_")

labels<-read_csv(paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",
                        slide,"_labels.csv"))
i<-rownames(xenium.obj)
i<-i[str_detect(i,"fung|bact",negate=T)]
xenium.obj<-subset(xenium.obj,cells=labels$barcode,features=i)

se <- SpatialExperiment(assay = list(counts = xenium.obj@assays$counts$counts), spatialCoords = xenium.obj[["fov"]]$centroids@coords)
rm(xenium.obj)
gc()
# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE)
set.seed(1000)
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = 0.2)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = 0.2)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = 0.2, resolution = 2)#lambda 0.8 for domains
temp<-reducedDim(se,"UMAP_M1_lam0.2")
temp<-temp%>%rownames_to_column("barcode")
#temp<-unsm%>%select(barcode,labels)%>%inner_join(temp)
temp<-se@colData%>%as.data.frame()%>%rownames_to_column("barcode")%>%inner_join(temp)
data.table::fwrite(temp,
                   file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnskyCl.gz"))