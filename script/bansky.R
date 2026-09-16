#!/usr/bin/env Rscript

start_time <- Sys.time()
options(future.globals.maxSize= 5242880000)
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
#slide<-"/media/Xenium_On_NAS/xenium_reanalysis/042825_Batch4Endometrioma_Rerun-EDV013-Left_EDV002_EDV010_EDV013/outs/"
transcripts <- arrow::read_parquet(paste0(slide,"/transcripts.parquet"))
coords<-data.table::fread(paste0(slide,"/cells.csv.gz"))
slide<-str_remove_all(slide,".+sis.|.ou.+")
print(paste("Loaded",slide)) 

transcripts <- transcripts %>% 
  filter(qv >= 20 & overlaps_nucleus == 1 & cell_id != "UNASSIGNED")
transcripts<-transcripts%>%
  filter(str_detect(feature_name,"Code|fung|bact|Contro",negate = T))
mat<-transcripts%>%
  dplyr::count(cell_id,feature_name)%>%
  pivot_wider(names_from = cell_id,values_from = n,values_fill = 0)%>%
  column_to_rownames("feature_name")
mat<-Matrix::Matrix(as.matrix(mat))
sc<-CreateSeuratObject(counts = mat)
rm(transcripts);gc()

qced<-data.table::fread(paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/Banksy_matrix/",
                               slide,"l0.8_bnksy_mtrx.gz"))#comment for 5k panel
#uncomment for 5k panel
#qced<-data.table::fread("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_domains5k.gz")#for 5k panel
#qced<-qced%>%filter(Section==slide)%>%mutate(cell=str_remove(cell_id,slide)%>%str_remove("_"))#for 5k panel
sc<-subset(sc,cells = qced$cell)
rm(mat,qced);gc()

coords<-coords%>%filter(cell_id%in%colnames(sc))
coords<-coords[order(match(coords$cell_id,colnames(sc))),]
sc<-AddMetaData(sc,
                metadata = coords%>%column_to_rownames("cell_id"))
#print("Seurat object ready")

#sc<-NormalizeData(sc,scale.factor = median(sc$nucleus_area))
f<-median(sc$nucleus_area)/sc$nucleus_area
sc@assays$RNA["data"]<-t(t(sc@assays$RNA$counts)*f)
#sc<-ScaleData(sc)
set.seed(1000)
sc@assays$RNA$data<-as.matrix(sc@assays$RNA$data)
sc <- RunBanksy(sc, lambda = .2, verbose=TRUE, 
                assay = 'RNA',features = "all",
                use_agf=T,dimx = "x_centroid",dimy = "y_centroid")
sc <- RunPCA(sc, assay = 'BANKSY',
             rownames(sc),npcs = 30)
#sc<-RunUMAP(sc,dims=1:20,spread = 3,min_dist = 0.1)
sc <- FindNeighbors(sc, dims = 1:30,k.param = 5)
#clusterBanksy(se, use_agf = TRUE, lambda = 0.8, resolution=c(.1))
gc()
#se <-  Banksy::connectClusters(se)
sc <- FindClusters(sc, resolution = 0.1,algorithm = 4)
cls<-unique(sc$BANKSY_snn_res.0.1)
sc$clus<-sc$seurat_clusters
for (x in cls) {
  sc<-FindSubCluster(sc,graph.name = "BANKSY_snn",algorithm = 4,cluster = x)
  sc@meta.data<-sc@meta.data%>%
    mutate(clus=ifelse(str_detect(sub.cluster,"_"),sub.cluster,clus))
}

bulkref<-readRDS("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/bulkref.RDS")
bulkref<-subset(bulkref,subset=label!="Erythrocytes")
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(sublabel=str_remove_all(label," \\(.+|Early ")%>%
           str_replace("emmo","emo"))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(toplabel=case_when(str_detect(label,"Smoo|Fibr|Mesen|GAS")~"Mesenchymal",
                            str_detect(label,"Th")~"CD4 T-cells",
                            str_detect(label,"CD8")~"CD8 T-cells",
                            str_detect(label,"Gamm|NK")~"Gamma_NK",
                            str_detect(label,"B-|Plas")~"B_Plasma",
                            str_detect(label,"EnEp|MUC5B|SOX9|IHH|Cili|Gland|KRT|Meso")~"Epithelial",
                            str_detect(label,"Macro|Dend")~"Myeloid",
                            TRUE~label))
get_labels<-function(scobj,cluster_column,level){
  ref<-bulkref
  bulkquery<-AggregateExpression(scobj,group.by = cluster_column,return.seurat = T)
  bulkquery<-ScaleData(bulkquery,assay = "RNA",features = i)
  if(level=="sublabel"){
    ct<-scobj$percl_toplabel[1]
    ref<-subset(bulkref,subset=toplabel==ct)
  }
  res <- SingleR(test=bulkquery@assays$RNA$scale.data,#this
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
#sc$res0.5Label<-labels_percl$pruned.labels
i<-rownames(sc@assays$RNA)
labels_percl<-get_labels(sc,"clus","toplabel")
sc$percl_toplabel<-labels_percl$pruned.labels

#subgroup
subxen<-lapply(unique(sc$percl_toplabel),function(x) subset(sc,subset=percl_toplabel==x))
names(subxen)<-unique(sc$percl_toplabel)
subxen<-subxen[sapply(subxen,function(x) length(unique(x$clus)))>1]
subxen<-subxen[names(subxen)!="unknown"]
labels_percl<-lapply(subxen,function(x) get_labels(x,"clus","sublabel"))
labels_percl<-labels_percl%>%bind_rows()%>%
  mutate(percl_sublabel=pruned.labels)%>%distinct(cell,percl_sublabel)
labels_percl<-sc@meta.data%>%rownames_to_column("cell")%>%left_join(labels_percl)
labels_percl<-labels_percl[order(match(labels_percl$cell,colnames(sc))),]
labels_percl<-labels_percl%>%
  mutate(percl_sublabel=ifelse(is.na(percl_sublabel),percl_toplabel,percl_sublabel))
sc$percl_sublabel<-labels_percl$percl_sublabel

png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/",slide,"bnksy_toplabels.png"),
    width=1200,height=300)
DotPlot(sc,group.by = "percl_toplabel",features = j)+RotatedAxis()+
  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(size=8))
dev.off()
png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/",slide,"bnksy_sublabels.png"),
    width=1200,height=400)
DotPlot(sc,group.by = "percl_sublabel",features = j)+RotatedAxis()+
  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(size=8))
dev.off()

labels_percl%>%select(cell,clus,percl_toplabel,percl_sublabel)%>%
  data.table::fwrite(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/",
                                   slide,"bnksy_labels.gz"))

end_time <- Sys.time()
end_time-start_time


