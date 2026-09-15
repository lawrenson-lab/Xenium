#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(dbscan)


files<-list.files("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_matrix",full.names = T)
files<-files[str_detect(files,"Pro")]
files<-files[str_detect(files,"UTP",negate = T)]

mats<-lapply(files,data.table::fread)
names(mats)<-str_remove_all(files,".+rix.|l0.8.+")
mats<-lapply(mats,function(x) x%>%column_to_rownames("cell")%>%t())
cs<-lapply(names(mats),function(x) paste0(x,'_',colnames(mats[[x]])))
concate<-do.call(cbind,mats)
colnames(concate)<-unlist(cs)
rm(mats,cs)
###############################################################################################
#                                         ADDING COORDINATES
###############################################################################################
files<-str_remove_all(files,".+rix.|l0.8.+")
files<-paste0("/media/Xenium_On_NAS/xenium_reanalysis/",files,"/outs/cells.csv.gz")

coords<-lapply(files,data.table::fread)
names(coords)<-str_remove_all(files,".+sis.|.outs.+")

temp<-coords$`022426_TMA-Endo_GE-Pro_Ba1-TMA-1-Right_TMA`
temp$cl<-dbscan(coords$`022426_TMA-Endo_GE-Pro_Ba1-TMA-1-Right_TMA`[,2:3],eps = 100)$cluster
temp<-temp%>%group_by(cl)%>%summarise(max=max(y_centroid),min=min(y_centroid),size=max-min)%>%filter(cl!=0)
sepd<-mean(temp$size)

coords<-lapply(names(coords),function(x) coords[[x]]%>%mutate(Section=x))%>%bind_rows()
#paste all together?  circle will be rotated
coords<-coords%>%
  mutate(y_centroid=case_when(str_detect(Section,"TMA-3")~y_centroid+14*sepd,
                              str_detect(Section,"TMA-4")~y_centroid+6*sepd,
                              TRUE~y_centroid),
         x_centroid=case_when(str_detect(Section,"TMA-1")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-1")])),
                              str_detect(Section,"TMA-2")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-2")]))+5*sepd,
                              str_detect(Section,"TMA-3")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-3")]))+10*sepd,
                              str_detect(Section,"TMA-4")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-4")]))+8*sepd,
                              TRUE~x_centroid))
coords<-coords%>%mutate(cell_id=paste0(Section,"_",cell_id))%>%filter(cell_id%in%colnames(concate))
coords<-coords[order(match(coords$cell_id,colnames(concate))),]
###############################################################################################
#                                   REDUCTIONS
###############################################################################################
pcres<-irlba::irlba(A = t(concate),nv=30)
hares<-harmony::RunHarmony(pcres$u,meta_data=coords$Section)
#umap_results <- uwot::umap(hares)
rownames(hares)<-colnames(concate)
colnames(hares)<-paste0("pc",1:30)
neighs<-FindNeighbors(hares,k.param = 30)
clus<-FindClusters(neighs$snn, resolution = 0.5,algorithm = 4)

coords<-clus%>%rownames_to_column("cell_id")%>%inner_join(coords)
coords%>%data.table::fwrite("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_domains.gz")
