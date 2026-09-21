#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(Banksy)
library(SeuratWrappers)
library(dbscan)
files<-list.files("/media/Xenium_On_NAS/xenium_reanalysis",full.names = T)
files<-files[str_detect(files,"TMA")]
files<-files[1:8]#5k panel, all come from the same slide


###############################################################################################
#                           LOADING AND FILTERING
###############################################################################################
xenium.list<-lapply(files, function(x) Read10X_h5(paste0(x,"/outs/cell_feature_matrix.h5")))
xenium.list<-lapply(xenium.list,function(x) x$`Gene Expression`)
names(xenium.list)<-str_remove(files,".+sis.")

xenium.list<-lapply(xenium.list,function(x) CreateSeuratObject(counts = x))
xenium.list<-lapply(xenium.list,function(x) subset(x,subset=nFeature_RNA>4))#10?
model <- lapply(xenium.list,function(x)
  lm(log10(nCount_RNA)~log10(nFeature_RNA), data = x@meta.data))
thr<-lapply(model,function(x) mean(x$residuals) + (3 * sd(x$residuals)))
xenium.list<-lapply(names(xenium.list),function(x) 
  subset(xenium.list[[x]],
         cells=which(model[[x]]$residuals<thr[[x]]&model[[x]]$residuals>(-1*thr[[x]]))))
names(xenium.list)<-names(model)
rm(model,thr);gc()

xenium.list<-lapply(xenium.list,NormalizeData)
xenium.list<-lapply(xenium.list,FindVariableFeatures,nfeatures=500)#low so banksy is small and easy to handle
i<-lapply(xenium.list,VariableFeatures)
i<-unique(unlist(i))

###############################################################################################
#                                         ADDING COORDINATES
###############################################################################################
coords<-lapply(files,function(x) data.table::fread(paste0(x,"/outs/cells.csv.gz")))
names(coords)<-str_remove(files,".+sis.")

temp<-coords$`020626_Prime_TMA_Endo_ba1_rerun-TMA-1-Right_TMA`
temp$cl<-dbscan(coords$`020626_Prime_TMA_Endo_ba1_rerun-TMA-1-Right_TMA`[,2:3],eps = 100)$cluster
temp<-temp%>%group_by(cl)%>%summarise(max=max(y_centroid),min=min(y_centroid),size=max-min)%>%filter(cl!=0)
sepd<-mean(temp$size)

coords<-lapply(names(coords),function(x) coords[[x]]%>%mutate(Section=x))%>%bind_rows()
#paste all together?  circle will be rotated
coords<-coords%>%
  mutate(y_centroid=case_when(str_detect(Section,"TMA-2")~y_centroid+max(y_centroid[str_detect(Section,"TMA-1")])+.5*sepd,
                              str_detect(Section,"TMA-3")~y_centroid+max(y_centroid[str_detect(Section,"TMA-2")])+5*sepd,
                              str_detect(Section,"TMA-4")~y_centroid+sepd,str_detect(Section,"TMA-5")~y_centroid-sepd,
                              str_detect(Section,"TMA-6")~y_centroid+max(y_centroid[str_detect(Section,"TMA-5")])-sepd,
                              str_detect(Section,"TMA-7")~y_centroid+max(y_centroid[str_detect(Section,"TMA-1")])+2*sepd,
                              str_detect(Section,"TMA-8")~y_centroid+max(y_centroid[str_detect(Section,"TMA-1")])+3*sepd,
                              TRUE~y_centroid),
         x_centroid=case_when(str_detect(Section,"TMA-1")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-1")])),
                              str_detect(Section,"TMA-2")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-2")])),
                              str_detect(Section,"TMA-3")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-3")])),
                              str_detect(Section,"TMA-4")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-4")]))+2*sepd,
                              str_detect(Section,"TMA-5")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-5")]))+3*sepd,
                              str_detect(Section,"TMA-6")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-6")]))+4.5*sepd,
                              str_detect(Section,"TMA-7")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-7")]))+8*sepd,
                              str_detect(Section,"TMA-8")~abs(x_centroid-max(x_centroid[str_detect(Section,"TMA-8")]))+10*sepd,
                              TRUE~x_centroid))
###############################################################################################
#                                   BANKSY
###############################################################################################
#targeted_mat<-data.table::fread("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_matrix/022426_TMA-Endo_GE-Pro_Ba1-TMA-1-Right_TMAl0.8_bnksy_mtrx.gz")
#i<-colnames(targeted_mat)
#j<-rownames(xenium.list$`020626_Prime_TMA_Endo_ba1_rerun-TMA-1-Right_TMA`)
#i<-intersect(i,j)

xenium.list<-lapply(xenium.list,subset,features=i)
inte<-merge(xenium.list[[1]],y = xenium.list[names(xenium.list)[-1]],
            add.cell.ids =names(xenium.list))
coords<-coords%>%mutate(cell_id=paste0(Section,"_",cell_id))
coords<-coords%>%filter(cell_id%in%colnames(inte))
coords<-coords[order(match(coords$cell_id,colnames(inte))),]
inte<-AddMetaData(inte,metadata = coords)

inte[["RNA"]]<-JoinLayers(inte[["RNA"]])
inte@assays$RNA$data<-as.matrix(inte@assays$RNA$data)#or banksy will fail

set.seed(1000)
inte <- RunBanksy(inte, lambda = .8, verbose=TRUE, assay = 'RNA',features = "all",
                  k_geom = 30,use_agf=T,dimx = "x_centroid",dimy="y_centroid")
inte@assays$BANKSY$data%>%t()%>%as.data.frame()%>%rownames_to_column("cell")%>%
  data.table::fwrite("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_matrix/5Kl0.8_bnksy_mtrx.gz")

mat<-inte@assays$BANKSY$data
rm(inte);gc()
pcres<-irlba::irlba(A = mat,nv=30)
hares<-harmony::RunHarmony(pcres$u,meta_data=coords$Section)
rownames(hares)<-rownames(mat)
colnames(hares)<-paste0("comp",1:30)
 umap_results <- uwot::umap(hares)
neighs<-FindNeighbors(hares,k.param = 30)
clus<-FindClusters(neighs$snn, resolution = 0.5,algorithm = 4)
coords<-clus%>%rownames_to_column("cell_id")%>%inner_join(coords)
coords%>%data.table::fwrite("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_domains5k.gz")
