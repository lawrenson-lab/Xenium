#!/usr/bin/Rscript --vanilla

library(Seurat)
library(SingleR)
library(tidyverse)
library(BiocParallel)
#library(Matrix)

slide<- commandArgs(trailingOnly=TRUE)
#slide<-"/media/Xenium_On_NAS/20250428__202040__042825_Batch4Endometrioma_Rerun/output-XETG00426__0034237__EDV013__20250428__202100"

xenium.obj <- LoadXenium(slide,fov="fov",
                         molecule.coordinates = F,
                         assay = "counts",
                         segmentations = "nucleus",
                         cell.centroids = F)
slide<-str_remove(slide,".+output-")
slide<-paste(unlist(str_split(slide,"__"))[3:2],collapse = "_")

i<-rownames(xenium.obj)
i<-i[str_detect(i,"fung|bact",negate=T)]
xenium.obj<-subset(xenium.obj,features=i)
xenium.obj$nFeature_Xenium<-colSums(xenium.obj@assays$Xenium$counts>0)
xenium.obj$nCount_Xenium<-colSums(xenium.obj@assays$Xenium$counts)
xenium.obj$ratio<-(xenium.obj$nCount_Xenium+0.01)/(xenium.obj$nFeature_Xenium+0.01)
xenium.obj<-subset(xenium.obj,subset=nFeature_Xenium>4)
j<-xenium.obj@meta.data%>%slice_max(ratio,n=100)%>%rownames()
j<-which(!colnames(xenium.obj)%in%j)
xenium.obj<-subset(xenium.obj,cells=j)
gc()
#scDblFinder??????????????????

xenium.obj <- SCTransform(xenium.obj,conserve.memory=TRUE,assay = "Xenium")
gc()
####### comment when getting unsmoothed labels##################
xenium.obj<-RunPCA(xenium.obj)
xenium.obj <- FindNeighbors(xenium.obj, dims = 1:10,k.param=5)
#xenium.obj<-FindClusters(xenium.obj,resolution = 5)
g <- xenium.obj@graphs$SCT_snn
expr <- xenium.obj@assays$SCT$data
expr_smooth <- g %*% t(expr)
#rm(xenium.obj)
gc()
################################################################

#############ran once and saved as bulkref.RDS##################
#reference<-readRDS("/media/cedars/caipirinha/speed/share/fixed_annotated_aux.seurat.rds")
#reference<-subset(reference,cells=which(!reference$harmonized_major_cell_type%in%
#                                          c("Exclude","Erythrocytes")),
#                  features=i)
#gc()
#reference@meta.data<-reference@meta.data%>%
#  mutate(
#    label=ifelse(!is.na(published_epithelial_subtype)&str_detect(published_epithelial_subtype,"KRT",negate = T),
#                 as.character(as.factor(published_epithelial_subtype)),
#                 harmonized_major_plus_immune_cell_type),
#    label=ifelse(!is.na(published_mesenchymal_subtype),
#                 as.character(as.factor(published_mesenchymal_subtype)),
#                 label)%>%
#      str_replace("_Muscle"," muscle cells")%>%
#      str_replace("_Pro.+|_Sec.+"," cells")%>%
#      trimws()%>%str_remove(" \\(.+"))
#reference@meta.data<-reference@meta.data%>%
#  mutate(nosub=case_when(str_detect(label,"Fibr|Mesen|GAS")~"Mesenchymal",
#                         str_detect(label,"T-|NK")~"T/NK",
#                         str_detect(label,"Macro|Dend")~"Myeloid",
#                         str_detect(label,"EnEp|MUC5B|SOX9|IHH|Cili|Gland")~"EnEpi",
#                         str_detect(label,"B|Plas")~"B/Plasma",TRUE~label))

#g<-reference@graphs$SCT_snn
#temp<-reference@meta.data%>%rownames_to_column("cell")%>%select(cell,nosub)%>%table()
#temp<-Matrix::Matrix(temp)
#temp<-temp[rownames(g),]
#temp<-g%*%temp
#i<-rowSums(temp>0)
#reference<-subset(reference,cells=names(which(i<3)))
#rm(g,temp)
#gc()
#
#bulkref<-AggregateExpression(reference,group.by = c("SampleName","label"),return.seurat = T)
#rm(reference);gc();
##############################################################

bulkref<-readRDS("bulkred.RDS")

#remove genes with poor correlation????????????
bp <- MulticoreParam(4)
set.seed(112358)

#start_time <- Sys.time()
res <- SingleR(test=xenium.obj@assays$SCT$data,#unsmoothed labels
	       #test=t(expr_smooth),#smoothed labels
               ref=bulkref@assays$SCT$data, 
               labels=bulkref$label, method = "cosine",
               #clusters = xenium.obj$seurat_clusters,
               BPPARAM=bp,quantile = 0.999,tune.thresh = .1)
#end_time <- Sys.time()

#fix filename smoothed vs unsmoothed
res%>%as.data.frame()%>%rownames_to_column("barcode")%>%
  data.table::fwrite(file=paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"_labels.csv"))

#xenium.obj<-RunUMAP(xenium.obj,dims = 1:10,n.neighbors = 50,min.dist = .001,spread = .5)
