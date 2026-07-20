#!/usr/bin/Rscript --vanilla

library(Seurat)
library(tidyverse)
#library(Matrix)
cellcols <- c("gray13", "gold2", "plum4", "darkorange1", "lightskyblue2", "firebrick",
              "burlywood3", "gray51", "springgreen4", "lightpink2", "deepskyblue4", 
              "lightsalmon2", "mediumpurple4", "orange", "maroon", "yellow3", "brown4", 
              "yellow4", "sienna4", "chocolate","#9e0142", "#d53e48",
              "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd",
              "cornflowerblue","forestgreen","hotpink","green3","blue4","brown","orange4",
              "yellow1","violet","purple2")
nosubcols<-c("B/Plasma"="gold2", "Endothelial cells"="plum4", "EnEpi"="darkorange1","T/NK"= "deepskyblue4",
             "EnS cells"="lightskyblue2", "Mast cells"="firebrick",Mesenchymal="burlywood3",
             Mesothelial="gray51", Myeloid="springgreen4", "Smooth muscle cells"="lightpink2",
             "low counts"="yellow1",non_human="hotpink","large cell"="purple2")

slide<- commandArgs(trailingOnly=TRUE)
#slide<-"/media/Xenium_On_NAS/20250201__003824__013125_Endometrium_Batch1/output-XETG00426__0033952__EDV003__20250201__003842/"
xenium.obj <- LoadXenium(slide,fov="fov")

slide<-str_remove(slide,".+output-")
slide<-paste(unlist(str_split(slide,"__"))[3:2],collapse = "_")
labels<-read_csv(paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"_labels.csv"))

#i<-rownames(xenium.obj)
#i<-i[str_detect(i,"fung|bact",negate=T)]
#xenium.obj<-subset(xenium.obj,features=i)
#xenium.obj$nFeature_Xenium<-colSums(xenium.obj@assays$Xenium$counts>0)
#xenium.obj$nCount_Xenium<-colSums(xenium.obj@assays$Xenium$counts)
#xenium.obj$ratio<-(xenium.obj$nCount_Xenium+0.01)/(xenium.obj$nFeature_Xenium+0.01)
#j<-xenium.obj@meta.data%>%slice_max(ratio,n=100)%>%rownames()

#xenium.obj@meta.data<-xenium.obj@meta.data%>%
#  mutate(non_human_only=nFeature_Xenium==0,
#         nCount_outl=colnames(xenium.obj)%in%j,
#         low_human_feats=nFeature_Xenium>0&nFeature_Xenium<5)

labels<-labels%>%mutate(smoo_labels=ifelse(is.na(smoo_labels),labels,smoo_labels))
labels<-labels%>%mutate(smoo_labels=str_remove_all(smoo_labels,"scores.| [1-8]$| \\(.+|\\+")%>%
                          str_replace_all('-|\\.+',' ')%>%trimws())
#labels<-labels%>%mutate(barcode=str_remove(barcode,".+_"))

#temp<-xenium.obj@meta.data%>%rownames_to_column("barcode")
#temp<-temp%>%left_join(labels)
#temp<-temp%>%mutate(nosub=case_when(non_human_only==T~"non_human",
#                                    nCount_outl==T~"high nCount",
#                                    low_human_feats==T~"low counts",
#                                    TRUE~nosub),
#                    integrated_label=case_when(non_human_only==T~"non_human",
#                                    nCount_outl==T~"high nCount",
#                                    low_human_feats==T~"low counts",
#                                    TRUE~integrated_label))
#temp<-temp[order(match(temp$barcode,colnames(xenium.obj))),]
xenium.obj<-subset(xenium.obj,cells=labels$barcode)
labels<-labels[order(match(labels$barcode,colnames(xenium.obj))),]
xenium.obj<-AddMetaData(xenium.obj,metadata=labels,col.name = "smoo_labels")
#xenium.obj<-subset(xenium.obj,subset = !is.na(nosub))

#xenium.obj<-SCTransform(xenium.obj,assay = "Xenium")
#xenium.obj<-RunPCA(xenium.obj)
#xenium.obj<-FindNeighbors(xenium.obj)
#xenium.obj<-FindClusters(xenium.obj,resolution = 1)

#temp<-xenium.obj@meta.data%>%count(seurat_clusters,smoo_labels)%>%
#  group_by(seurat_clusters)%>%mutate(n=100*n/sum(n))%>%
#  slice_max(n)
#temp<-temp%>%filter(n<50)
#xenium.obj@meta.data$clus<-xenium.obj@meta.data$seurat_clusters
#for(i in temp$seurat_clusters){
#  xenium.obj<-FindSubCluster(xenium.obj,cluster = 0,graph.name = "SCT_snn")
#  xenium.obj@meta.data<-xenium.obj@meta.data%>%
#    mutate(clus=ifelse(str_detect(sub.cluster,"_"),sub.cluster,clus))}

#temp<-xenium.obj@meta.data%>%count(clus,smoo_labels)%>%
#  group_by(clus)%>%mutate(n=100*n/sum(n))%>%
#  slice_max(n)
#temp<-temp%>%
#  mutate(nosub=case_when(str_detect(smoo_labels,"Fibr|Mesen|GAS")~"Mesenchymal",
#                         str_detect(smoo_labels,"T |NK")~"T/NK",
#                         str_detect(smoo_labels,"Macro|Dend")~"Myeloid",
#                         str_detect(smoo_labels,"EnEp|MUC5B|SOX9|IHH|Cili|Gland")~"EnEpi",
#                         str_detect(smoo_labels,"B |Plas")~"B/Plasma",TRUE~smoo_labels))
#temp<-temp%>%group_by(clus)%>%mutate(total=sum(n))%>%group_by(clus,nosub)%>%
#  mutate(nosub_prop=100*sum(n)/total)%>%filter(nosub_prop>50)%>%
#  group_by(clus)%>%slice_max(n)%>%mutate(cl_ct=smoo_labels)%>%
#  select(smoo_labels,nosub)
#temp<-temp%>%add_count(clus)%>%mutate(cl_ct=ifelse(n==1,smoo_labels,nosub))%>%distinct(cl_ct)%>%ungroup()
#temp<-xenium.obj@meta.data%>%rownames_to_column("cell")%>%
#  select(cell,clus)%>%left_join(temp)
#temp<-temp[order(match(temp$cell,colnames(xenium.obj))),]
#xenium.obj<-AddMetaData(xenium.obj,metadata = temp,col.name = "cl_ct")
#xenium.obj@meta.data<-xenium.obj@meta.data%>%mutate(cl_ct=ifelse(is.na(cl_ct),smoo_labels,cl_ct))
xenium.obj@meta.data<-xenium.obj@meta.data%>%
  mutate(nosub=case_when(str_detect(smoo_labels,"Fibr|Mesen|GAS")~"Mesenchymal",
                        str_detect(smoo_labels,"T |NK")~"T/NK",
                        str_detect(smoo_labels,"Macro|Dend")~"Myeloid",
                        str_detect(smoo_labels,"EnEp|MUC5B|SOX9|IHH|Cili|Gland")~"EnEpi",
                        str_detect(smoo_labels,"B |Plas")~"B/Plasma",TRUE~smoo_labels))



psize=ifelse(ncol(xenium.obj)>=100000,.1,1)
png(paste0("docs/",slide,"_nosub.png"),width=1200,height = 1000)
ImageDimPlot(xenium.obj,group.by = "nosub",axes = T,dark.background = F,size = psize,cols=nosubcols,coord.fixed = F)+
         ggtitle(slide)
dev.off()       
png(paste0("docs/",slide,"_sub.png"),width=1200,height = 1000)
ImageDimPlot(xenium.obj,group.by = "smoo_labels",axes = T,dark.background = F,size = psize,cols=cellcols,coord.fixed = F)+
  ggtitle(slide)
dev.off()       

bulk<-AggregateExpression(
  xenium.obj,
  assays = "Xenium",
  return.seurat = T,
  #  group.by = "integrated_label")
  group.by = "smoo_labels")

mat<-GetAssayData(bulk,assay="Xenium",layer = "counts")
mat%>%as.data.frame()%>%rownames_to_column("gene")%>%
  write_csv(paste0("docs/",slide,"_pseudobulk.csv"))
