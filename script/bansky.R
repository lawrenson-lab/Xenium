library(Seurat)
library(tidyverse)
library(Banksy)
library(SeuratWrappers)
#library(scuttle)
library(SingleR)
library(BiocParallel)
j<-c("KIT","TPSAB1","SIGLEC8","PTPRC","CD79A","BANK1","MS4A1",
     "CD27","CD38","BCL6","CXCR5","MME","JCHAIN","PRDM1","SDC1",
     "IGHA1","IGHG1","IGHG2","IGHG3","IGHG4","IGHM","NCAM1","FCGR3A",
     "CD3D","TRDC","TRGC1","TRGC2","CD4","IL4","GATA3","CCR4","CXCR4",
     "RORC","IL17A","STAT3","CCR6","FOXP3","IL2RA","CD8A","SELL",
     "CCR7","CD44","IL7R","CD69","CTLA4","LAG3","PDCD1","TIGIT","HAVCR2",
     "TOX","CD1C","HLA-DRA","HLA-DPA1","ITGAM","CD68","CXCL10","CCL2",
     "NOS2","ARG1","CD163","IL10","EPCAM","KRT10","KRT18","KRT8","WT1",
     "PAX8","ESR1","FOXJ1","PAEP","PECAM1","PGR","IGF1","MMP11","FOXO1",
     "DCN","PDGFRA","THY1","COL1A1","ACTA2","TAGLN","MYL9","FAP","C7","GAS5")


slide<- commandArgs(trailingOnly=TRUE)
#slide<-"/media/Xenium_On_NAS/20260424__165225__042426_Endometrioma_Batch15_Seg/output-XETG00426__0102750__EDV059__20260424__165251"#v4 & myeloid issue

xenium.obj <- LoadXenium(slide,fov="fov",
                         molecule.coordinates = F,
                         assay = "counts",
                         segmentations = "nucleus",
                         cell.centroids = T)
slide<-str_remove(slide,".+output-")
slide<-paste(unlist(str_split(slide,"__"))[3:2],collapse = "_")

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
bulkref<-subset(bulkref,cells=which(bulkref$label%in%c("Erythrocytes","Mesothelial")))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(sublabel=str_remove_all(label," \\(.+|Early ")%>%
           str_replace("emmo","emo")%>%str_replace("Th.+","CD4 T-cells"))%>%
  mutate(sublabel=if_else(str_detect(sublabel,"Macro"),"Macrophages",
                          sublabel))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(toplabel=case_when(str_detect(label,"Smoo|Fibr|Mesen|GAS")~"Mesenchymal",
                            str_detect(label,"Th")~"CD4 T-cells",
                            str_detect(label,"CD8")~"CD8 T-cells",
                            str_detect(label,"Gamm|NK")~"Gamma_NK",
                            str_detect(label,"B-|Plas")~"B_Plasma",
                            str_detect(label,"EnEp|MUC5B|SOX9|IHH|Cili|Gland")~"Epithelial",
                            str_detect(label,"Macro|Dend")~"Myeloid",
                            TRUE~label))
get_labels<-function(scobj,cluster_column,level){
  bulkquery<-AggregateExpression(scobj,group.by = cluster_column,return.seurat = T)
  bulkquery<-ScaleData(bulkquery,assay = "counts",features = i)
  res <- SingleR(test=bulkquery@assays$counts$scale.data,#this
                 ref=bulkref@assays$SCT$data, 
                 labels=bulkref@meta.data[,level], method = "cosine",
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

i<-xenium.obj@meta.data%>%
  split(f = xenium.obj@meta.data$percl_toplabel)%>%
  lapply(function(x) rownames(x))
g<-lapply(i,function(x) xenium.obj@graphs$BANKSY_snn[x,x])
expr<-xenium.obj@assays$counts$scale.data
expr<-lapply(i,function(x) expr[,x])
expr <- lapply(names(i), function(x) g[[x]] %*% t(expr[[x]]))
names(expr)<-names(i)
rm(g);gc()
#bp <- MulticoreParam(6)
labels_percl <- lapply(expr,function(x)
  SingleR(test=t(x),#this
          ref=bulkref@assays$SCT$scale.data, 
          labels=bulkref$sublabel, method = "cosine",
          BPPARAM=bp,quantile = 0.999,tune.thresh = .1)%>%
    as.data.frame()%>%
    mutate(pruned.labels=ifelse(is.na(pruned.labels),
                                "unknown",pruned.labels)))

labels_percl<-labels_percl%>%bind_rows()%>%rownames_to_column("cell")
labels_percl<-xenium.obj@meta.data%>%rownames_to_column("cell")%>%
  select(cell,clus,percl_toplabel)%>%inner_join(labels_percl)
labels_percl<-labels_percl%>%#pivot_longer(scores.Activated.Fibro:scores.Th2.T.cells,
  #     names_to = "ct",values_to = "score")%>%
  mutate(nosub=case_when(str_detect(pruned.labels,"Smoo|Fibr|Mesen|GAS|EnS")~"Mesenchymal",
                         str_detect(pruned.labels,"Th")~"CD4 T-cells",
                         str_detect(pruned.labels,"CD8")~"CD8 T-cells",
                         str_detect(pruned.labels,"Gamm|NK")~"Gamma_NK",
                         str_detect(pruned.labels,"B-|Plas")~"B_Plasma",
                         str_detect(pruned.labels,"EnEp|MUC5B|SOX9|IHH|Cili|Gland|Meso|KRT")~"Epithelial",
                         str_detect(pruned.labels,"Macro|Dend")~"Myeloid",
                         TRUE~pruned.labels))





base_res<-0.5
while(flag>0|base_res>0){
  base_res<-base_res-.1
  xenium.obj <- FindClusters(xenium.obj, resolution = base_res,algorithm = 4)
  print(paste("Checking resolution",base_res))
  labels_expre<-get_labels(paste0("BANKSY_snn_res.",base_res))
  flag<-labels_expre%>%filter(str_detect(id,"Mese|^E",negate = T)&
                                str_detect(features.plot,"EPCAM|KRT")&
                                avg.exp.scaled>.1)%>%nrow()
  gc()
}


bulkquery<-AggregateExpression(xenium.obj,group.by = c("clust_M1_lam0.2_k15_res2"),return.seurat = T)
bulkquery<-ScaleData(bulkquery,assay = "counts",features = i)
res <- SingleR(test=bulkquery@assays$counts$scale.data,#this
               ref=bulkref@assays$SCT$data, 
               labels=bulkref$sublabel, method = "cosine",
               BPPARAM=bp,quantile = 0.999,tune.thresh = .1)
res<-as.data.frame(res)
temp$clust_M1_lam0.2_k15_res2<-as.numeric(temp$clust_M1_lam0.2_k15_res2)
temp<-res%>%rownames_to_column("clust_M1_lam0.2_k15_res2")%>%
  mutate(clust_M1_lam0.2_k15_res2=str_remove(clust_M1_lam0.2_k15_res2,"g")
         %>%as.numeric(),
         res2Labels=str_replace_na(pruned.labels,"missing"))%>%
  select(clust_M1_lam0.2_k15_res2,res2Labels)%>%inner_join(temp)

temp<-temp[order(match(temp$cell,colnames(xenium.obj))),]
xenium.obj<-AddMetaData(xenium.obj,metadata = temp)

png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnksy_labels.png"),
    width=1200,height=400)
DotPlot(xenium.obj,group.by = "res0.5Label",features = j)+RotatedAxis()+
  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(size=8))
dev.off()

temp%>%data.table::fwrite(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnksy_labels.gz"))