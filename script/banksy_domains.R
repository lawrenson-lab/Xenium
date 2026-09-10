#!/usr/bin/env Rscript

start_time <- Sys.time()
options(future.globals.maxSize= 5242880000)
library(Banksy)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
#library(scuttle)
j<-c("KIT","TPSAB1","SIGLEC8","PTPRC","CD79A","BANK1","MS4A1",
     "CD27","CD38","BCL6","CXCR5","MME","JCHAIN","PRDM1","SDC1",
     "IGHA1","IGHG1","IGHG2","IGHG3","IGHG4","IGHM","NCAM1","FCGR3A",
     "CD3D","TRDC","TRGC1","TRGC2","CD4","IL4","GATA3","CCR4","CXCR4",
     "RORC","IL17A","STAT3","CCR6","FOXP3","IL2RA","CD8A","SELL",
     "CCR7","CD44","IL7R","CD69","CTLA4","LAG3","PDCD1","TIGIT","HAVCR2",
     "TOX","CD1C","HLA-DRA","HLA-DPA1","ITGAM","CD68","CXCL10","CCL2",
     "NOS2","ARG1","CD163","IL10","EPCAM","KRT10","KRT18","KRT8","WT1",
     "MUC5B","SOX9","SPDEF","LGR5",
     "PAX8","ESR1","FOXJ1","PAEP","PECAM1","PGR","IGF1","MMP11","FOXO1",
     "DCN","PDGFRA","THY1","COL1A1","ACTA2","TAGLN","MYL9","FAP","C7","GAS5")


slide<- commandArgs(trailingOnly=TRUE)
#slide<-"/media/Xenium_On_NAS/xenium_reanalysis/013125_Endometrium_Batch1-EDV015-Left_EDV015_EDV003/outs/"

xenium.obj <- LoadXenium(slide,fov="fov",
                         molecule.coordinates = F,
                         assay = "counts",
                         segmentations = "cell",
                         cell.centroids = T)
slide<-str_remove_all(slide,".+sis.|.ou.+")
print(slide) 

i<-rownames(xenium.obj)
i<-i[str_detect(i,"fung|bact",negate=T)]
xenium.obj<-subset(xenium.obj,features=i)
xenium.obj$nFeature_Xenium<-colSums(xenium.obj@assays$counts$counts>0)#keep changing assay
xenium.obj<-subset(xenium.obj,subset=nFeature_Xenium>4)#10?
xenium.obj$nCount_Xenium<-colSums(xenium.obj@assays$counts$counts)
model <- lm(log10(nCount_Xenium)~log10(nFeature_Xenium), data = xenium.obj@meta.data)
xenium.obj$residuals<-residuals(model)
thr<-mean(xenium.obj$residuals) + (3 * sd(xenium.obj$residuals))
xenium.obj<-subset(xenium.obj,
                   cells=which(xenium.obj$residuals<thr&xenium.obj$residuals>(-1*thr)))
gc()

xenium.obj<-NormalizeData(xenium.obj,assay = "counts",normalization.method = "LogNormalize")
xenium.obj<-ScaleData(xenium.obj)

set.seed(1000)
xenium.obj@assays$counts$data<-as.matrix(xenium.obj@assays$counts$data)
xenium.obj <- RunBanksy(xenium.obj, lambda = .8, verbose=TRUE, 
                        assay = 'counts',features = "all",
                        k_geom = 30,
                        use_agf=T)
#xenium.obj <- RunPCA(xenium.obj, assay = 'BANKSY',
#                     features=rownames(xenium.obj),npcs = 30)
#xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30,k.param = 30)
#xenium.obj <- FindClusters(xenium.obj, resolution = 0.5,algorithm = 4)
#cls<-unique(xenium.obj$BANKSY_snn_res.0.5)
#xenium.obj$subdomain<-xenium.obj$seurat_clusters
#for (x in cls) {
#  xenium.obj<-FindSubCluster(xenium.obj,graph.name = "BANKSY_snn",algorithm = 4,cluster = x)
#  xenium.obj@meta.data<-xenium.obj@meta.data%>%
#    mutate(subdomain=ifelse(str_detect(sub.cluster,"_"),sub.cluster,subdomain))
#}

#xenium.obj@meta.data%>%rownames_to_column("cell")%>%
#  mutate(BANKSY_snn_lambda0.8_res.0.5=BANKSY_snn_res.0.5)%>%
#  select(cell,BANKSY_snn_lambda0.8_res.0.5,subdomain)%>%
#  data.table::fwrite(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/",
#                                   slide,"bnksy_domains.gz"))
#gc()
#png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/",slide,"bnksy_domains.png"),
#    width=1200,height=300)
#DotPlot(xenium.obj,group.by = "BANKSY_snn_res.0.5",features = j)+RotatedAxis()+
#  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
#  theme(axis.text.x = element_text(size=8))
#dev.off()
#png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/",slide,"bnksy_domainSpatial.png"),
#    width=300,height=300)
#ImageDimPlot(xenium.obj,group.by = "BANKSY_snn_res.0.5",cols="polychrome")
#dev.off()

temp<-xenium.obj@assays$BANKSY$scale.data
temp%>%t()%>%as.data.frame()%>%rownames_to_column("cell")%>%
  data.table::fwrite(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/",
                                   slide,"l0.8_bnksy_mtrx.gz"))

end_time <- Sys.time()
end_time-start_time


