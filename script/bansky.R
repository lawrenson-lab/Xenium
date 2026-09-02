#!/usr/bin/env Rscript

start_time <- Sys.time()
library(Banksy)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
#library(scuttle)
library(SingleR)
#library(BiocParallel)
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
#slide<-"/media/Xenium_On_NAS/20260424__165225__042426_Endometrioma_Batch15_Seg/output-XETG00426__0102750__EDV059__20260424__165251"#v4 & myeloid issue

xenium.obj <- LoadXenium(slide,fov="fov",
                         molecule.coordinates = F,
                         assay = "counts",
                         segmentations = "nucleus",
                         cell.centroids = T)
slide<-str_remove_all(slide,".+sis.|.ou.+")

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

#lambda <- c(0.2,0.8)##lambda 0.8 for domains
set.seed(1000)
xenium.obj@assays$counts$data<-as.matrix(xenium.obj@assays$counts$data)
xenium.obj <- RunBanksy(xenium.obj, lambda = .2, verbose=TRUE, 
                        assay = 'counts',features = "all",use_agf=T)
xenium.obj <- RunPCA(xenium.obj, assay = 'BANKSY', features = i, npcs = 30)
#xenium.obj<-RunUMAP(xenium.obj,dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30,k.param = 15)
#clusterBanksy(se, use_agf = TRUE, lambda = 0.8, resolution=c(.1))
gc()
#se <-  Banksy::connectClusters(se)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.5,algorithm = 4)
cls<-unique(xenium.obj$BANKSY_snn_res.0.5)
xenium.obj$clus<-xenium.obj$seurat_clusters
for (x in cls) {
  xenium.obj<-FindSubCluster(xenium.obj,graph.name = "BANKSY_snn",algorithm = 4,cluster = x)
  xenium.obj@meta.data<-xenium.obj@meta.data%>%
    mutate(clus=ifelse(str_detect(sub.cluster,"_"),sub.cluster,clus))
}

bulkref<-readRDS("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/bulkref.RDS")
bulkref<-subset(bulkref,cells=which(!bulkref$label%in%c("Erythrocytes","Mesothelial")))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(sublabel=str_remove_all(label," \\(.+|Early ")%>%
           str_replace("emmo","emo")%>%str_replace("Th.+","CD4 T-cells"))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(toplabel=case_when(str_detect(label,"Smoo|Fibr|Mesen|GAS")~"Mesenchymal",
                            str_detect(label,"Th")~"CD4 T-cells",
                            str_detect(label,"CD8")~"CD8 T-cells",
                            str_detect(label,"Gamm|NK")~"Gamma_NK",
                            str_detect(label,"B-|Plas")~"B_Plasma",
                            str_detect(label,"EnEp|MUC5B|SOX9|IHH|Cili|Gland|KRT")~"Epithelial",
                            str_detect(label,"Macro|Dend")~"Myeloid",
                            TRUE~label))
get_labels<-function(scobj,cluster_column,level){
  ref<-bulkref
  bulkquery<-AggregateExpression(scobj,group.by = cluster_column,return.seurat = T)
  bulkquery<-ScaleData(bulkquery,assay = "counts",features = i)
  if(level=="sublabel"){
    ref<-subset(bulkref,subset=toplabel==scobj$percl_toplabel[1])
  }
  res <- SingleR(test=bulkquery@assays$counts$scale.data,#this
                 ref=ref@assays$SCT$data, 
                 labels=ref@meta.data[,level], method = "cosine",
                 quantile = 0.999,tune.thresh = .1)
  res<-as.data.frame(res)%>%rownames_to_column(cluster_column)
  res<-res%>%mutate(pruned.labels=str_replace_na(pruned.labels,"unknown"))
  if(sum(str_detect(res[,cluster_column],"^g"))>0){
    res[,cluster_column]=str_remove(res[,cluster_column],"g")%>%str_replace("-","_")
  }
  res<-scobj@meta.data%>%rownames_to_column("cell")%>%right_join(res)
  res<-res[order(match(res$cell,colnames(scobj))),]
  return(res)}
#labels_percl<-get_labels("BANKSY_snn_res.0.5","toplabel")
#xenium.obj$res0.5Label<-labels_percl$pruned.labels
labels_percl<-get_labels(xenium.obj,"clus","toplabel")
xenium.obj$percl_toplabel<-labels_percl$pruned.labels

#subgroup
subxen<-lapply(unique(xenium.obj$percl_toplabel),function(x) subset(xenium.obj,subset=percl_toplabel==x))
names(subxen)<-unique(xenium.obj$percl_toplabel)
subxen<-subxen[sapply(subxen,function(x) length(unique(x$clus)))>1]
labels_percl<-lapply(subxen,function(x) get_labels(x,"clus","sublabel"))
labels_percl<-labels_percl%>%bind_rows()%>%
  mutate(percl_sublabel=pruned.labels)%>%distinct(cell,percl_sublabel)
labels_percl<-xenium.obj@meta.data%>%rownames_to_column("cell")%>%left_join(labels_percl)
labels_percl<-labels_percl[order(match(labels_percl$cell,colnames(xenium.obj))),]
labels_percl<-labels_percl%>%
  mutate(percl_sublabel=ifelse(is.na(percl_sublabel),percl_toplabel,percl_sublabel))
xenium.obj$percl_sublabel<-labels_percl$percl_sublabel

png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnksy_toplabels.png"),
    width=1200,height=400)
DotPlot(xenium.obj,group.by = "percl_toplabel",features = j)+RotatedAxis()+
  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(size=8))
dev.off()
#png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnksy_sublabels.png"),
#    width=1200,height=400)
#DotPlot(xenium.obj,group.by = "percl_sublabel",features = j)+RotatedAxis()+
#  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
#  theme(axis.text.x = element_text(size=8))
#dev.off()

labels_percl%>%select(cell,clus,percl_toplabel,percl_sublabel)%>%
  data.table::fwrite(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",
                                   slide,"bnksy_labels.gz"))


end_time <- Sys.time()
print(paste("Lasted",end_time-start_time))


