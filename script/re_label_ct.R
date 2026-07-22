#!/usr/bin/Rscript --vanilla

library(tidyverse)
library(Seurat)

ct<- commandArgs(trailingOnly=TRUE)
print(ct)

files<-list.files("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",full.names = T)
files<-files[str_detect(files,"EDV")]

annot<-lapply(files,read_csv)
names(annot)<-str_remove_all(files,".+\\/|_labe.+")
annot<-lapply(names(annot),function(x) annot[[x]]%>%mutate(Sample=x))%>%bind_rows()
annot<-annot%>%mutate(patient=str_remove(Sample,"_.+"))
annot<-annot%>%mutate(nosub=case_when(str_detect(labels,"Fibr|Mesen|GAS")~"Mesenchymal",
                                      str_detect(labels,"EnEp|MUC5B|SOX9|IHH|Cili|Gland")~"EnEpi",
                                      str_detect(labels,"T-|NK")~"T/NK",
                                      str_detect(labels,"Macro|Dend")~"Myeloid",
                                      str_detect(labels,"Meso")~"Mesothelial",
                                      str_detect(labels,"B-|Plas")~"B/Plasma",TRUE~labels))
i<-annot%>%filter(nosub==ct&str_detect(Sample,"unsm",negate = T))%>%distinct(barcode)%>%unlist()
annot<-annot%>%filter(barcode%in%i&str_detect(Sample,"unsm"))
annot<-annot%>%add_count(Sample)%>%filter(n>100)%>%select(-n)              
gc()

annot<-annot%>%mutate(Sample=str_remove(Sample,"unsm"))
#files<-list.files("/media/Xenium_On_NAS",full.names = T)
#files<-files[str_detect(files,"ndometrioma")]
#files<-list.files(files,full.names = T)
files<-readLines("temp")
i<-str_split_i(files,pattern = "/",i = 5)%>%str_split(pattern = "__")%>%lapply(function(x) paste(x[3:2],collapse = "_"))
i<-unlist(i)
j<-which(i%in%annot$Sample)
files<-files[j]
i<-i[j]

mats<-lapply(files,function(x) Read10X_h5(paste0(x,"/cell_feature_matrix.h5")))
names(mats)<-unlist(i)
mats<-lapply(mats,function(x) x[[1]])#filter control matrixes
i<-annot%>%split(f = annot$Sample)%>%lapply(function(x) str_remove(x$barcode,".+_"))
i<-i[names(mats)]
mats<-lapply(names(i),function(x) mats[[x]][,colnames(mats[[x]])%in%i[[x]]])#filter myeloid cells
names(mats)<-names(i)
gc()

sc<-lapply(mats,CreateSeuratObject)
rm(mats);gc()
inte<-merge(sc[[1]],y = sc[names(i)[-1]],add.cell.ids =names(i))
rm(sc);gc()
annot<-annot%>%mutate(barcode=paste0(Sample,'_',barcode))
annot<-annot%>%filter(barcode%in%colnames(inte))
annot<-annot[order(match(annot$barcode,colnames(inte))),]
inte<-AddMetaData(inte,metadata = annot)
i<-rownames(inte)
i<-i[str_detect(i,pattern = "fung|bacte",negate = T)]
inte<-subset(inte,features = i)
gc()


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
#bulkref<-AggregateExpression(reference,group.by = c("SampleName","label"),return.seurat = T)
#rm(reference);gc()
#degREF<-FindAllMarkers(bulkref,group.by = "label")
#write_csv(degREF,"bulkref.csv")
ref<-read_csv("bulkref.csv")
i<-annot%>%filter(nosub==ct)%>%distinct(labels)%>%unlist()
ref<-ref%>%filter(cluster%in%i)
i<-ref%>%filter(cluster%in%unlist(i)&avg_log2FC>0)%>%add_count(gene)%>%filter(n<6)%>%distinct(gene)%>%unlist()

inte<-NormalizeData(inte)
inte<-ScaleData(inte,features = i)
gc()

inte<-RunPCA(inte,features = i)#is not using all features
inte <- harmony::RunHarmony(inte, group.by.vars="Sample")
#inte[["RNA"]]<-JoinLayers(inte[["RNA"]])
#inte<-RunUMAP(inte,min.dist = 0.001,repulsion.strength = 3,reduction = "harmony",features = i)
gc()

inte<-FindNeighbors(inte,reduction = "harmony",features = i)
#scores<-scale(t(annot[,2:48]))
#scores[scores<0]<-0
#g<-inte@graphs$RNA_snn
#scores_smooth <- g %*% t(scores)
#inte$smoo_labels<-colnames(scores_smooth)[apply(scores_smooth,1,which.max)]
inte<-FindClusters(inte,algorithm = 4,resolution = .5)
temp<-inte@meta.data%>%count(seurat_clusters,labels)%>%group_by(seurat_clusters)%>%mutate(n=100*n/sum(n))%>%slice_max(n)
temp<-temp%>%filter(n<75)
inte@meta.data$clus<-inte$seurat_clusters
for(i in temp$seurat_clusters){
  inte<-FindSubCluster(inte,cluster = i,graph.name = "RNA_snn",algorithm = 4)
  inte@meta.data<-inte@meta.data%>%mutate(clus=ifelse(str_detect(sub.cluster,"_"),sub.cluster,clus))
}

percl<-inte@meta.data%>%count(clus,nosub)%>%group_by(clus)%>%mutate(n=100*n/sum(n))%>%slice_max(n)
i<-percl%>%filter(nosub==ct)%>%distinct(clus)%>%unlist()
temp<-percl%>%filter(clus%in%i)%>%select(-n)%>%inner_join(inte@meta.data)%>%
  mutate(labels=str_remove(labels," \\(.+"))%>%
  count(clus,labels)%>%group_by(clus)%>%slice_max(n)%>%
  mutate(cl_label=labels)%>%select(cl_label)%>%ungroup()
temp<-percl%>%filter(!clus%in%i)%>%select(-n)%>%inner_join(inte@meta.data)%>%
  mutate(labels=str_remove(labels," \\(.+"))%>%count(clus,labels)%>%group_by(clus)%>%
  slice_max(n)%>%mutate(cl_label=labels)%>%select(cl_label)%>%ungroup()%>%bind_rows(temp)
#solve duplicates
i<-temp%>%add_count(clus)%>%filter(n>1)%>%distinct(clus)%>%unlist()
temp<-temp%>%filter(!clus%in%i)
for(x in i){
  temp<-temp%>%filter(str_detect(clus,str_replace(x,"_[0-9]+$","_")))%>%
    count(cl_label)%>%slice_max(n)%>%select(cl_label)%>%mutate(clus=x)%>%
    bind_rows(temp)}
temp<-temp%>%add_count(clus)%>%mutate(cl_label=ifelse(n>1,"unknown",cl_label))%>%
  select(-n)%>%distinct()
#add
temp<-inte@meta.data%>%select(barcode,clus)%>%inner_join(temp)
temp<-temp[order(match(temp$barcode,colnames(inte))),]
inte<-AddMetaData(inte,metadata = temp,col.name = "cl_label")

png(file=paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",
                str_replace_all(ct,"/| ","_"),"_markers.png"),height = 600,width = 1100)
DotPlot(inte,group.by = "cl_label",features =c("CD79A","CR2","MS4A1","CD1C","HLA-DRA","CLEC10A","VSIG4","S100A9","FCN1","IL1A","IL1B","IL1RN","TPSAB1","SIGLEC8","KIT","GZMB","GNLY","KLRC2","IGHG2","IGHG1","IGHA1","CD3D","CD8A","CD69","IL2RA","IL17A","TRDC","TRGC1","TRGC2","FAP","FN1","DCN","COL1A1","C7","FOXJ1","CLU","EPCAM","A2M","ACTA2"))+RotatedAxis()+scale_color_gradient2(low="blue",mid="white",high = "red")
dev.off()

temp<-inte@meta.data%>%select(barcode,labels,cl_label,clus)
#temp<-inte@reductions$umap@cell.embeddings%>%bind_cols(temp)
data.table::fwrite(temp,file=paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",
                                    str_replace_all(ct,"/| ","_"),"_labels.csv"))