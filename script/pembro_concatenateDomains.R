#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)


files<-list.files("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_matrix",full.names = T)
files<-files[str_detect(files,"l0.8_bnksy_mtrx")]

mats<-lapply(files,data.table::fread)
names(mats)<-str_remove_all(files,".+rix.|l0.8.+")
i<-colnames(mats$`5K`)
mats<-lapply(mats,function(x) x%>%select(i))
temp<-lapply(1:4,function(x) mats[[x]]%>%mutate(cell=paste0(names(mats)[x],"_",cell),
                                                Section=names(mats)[x]))%>%bind_rows()
mats<-mats$`5K`%>%mutate(Section=str_remove(cell,"_[a-z]+\\-1$"))%>%bind_rows(temp)

########################################################
# CONCATENATED DOMAINS
########################################################
temp<-mats%>%select(cell,Section)
mats<-mats%>%column_to_rownames("cell")%>%select(-Section)%>%as.matrix()
pcres<-irlba::irlba(A = mats,nv=30)#crashed with 50
hares<-harmony::RunHarmony(pcres$u,meta_data=temp$Section)
rownames(hares)<-temp$cell
colnames(hares)<-paste0("comp",1:30)
# umap_results <- uwot::umap(hares)
neighs<-FindNeighbors(hares,k.param = 50)
clus<-FindClusters(neighs$snn, resolution = 0.1,algorithm = 4)#increase resolution?
temp<-clus%>%rownames_to_column("cell")%>%inner_join(temp)
temp%>%data.table::fwrite("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_domains_concate.gz")