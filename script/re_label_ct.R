#!/usr/bin/Rscript --vanilla

####################################LOADING####################################
library(tidyverse)
library(Seurat)
library(SingleR)
library(cowplot)
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

ct<- commandArgs(trailingOnly=TRUE)
print(ct)

files<-list.files("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels",full.names = T)
files<-files[str_detect(files,"bnksy_labels.gz")]

annot<-lapply(files,read_csv)
annot<-lapply(1:length(files),function(x) 
  annot[[x]]%>%mutate(Sample=str_remove_all(files[x],".+m_labels.|bnksy.+")))%>%bind_rows()
#annot<-annot%>%mutate(patient=str_remove(Sample,"_.+"))
annot<-annot%>%filter(percl_toplabel==ct)
i<-annot%>%dplyr::count(Sample)%>%filter(n>1)
annot<-annot%>%filter(Sample%in%i$Sample)
gc()

files<-list.files("/media/Xenium_On_NAS/xenium_reanalysis",full.names = T)
i<-str_remove(files,"^.+sis.")
files<-files[i%in%annot$Sample]
files<-paste0(files,"/outs")
i<-i[i%in%annot$Sample]
print(head(paste0(files,"/cell_feature_matrix.h5")))
mats<-lapply(files,function(x) Read10X_h5(paste0(x,"/cell_feature_matrix.h5")))
names(mats)<-unlist(i)
mats<-lapply(mats,function(x) x[[1]])#filter control matrixes
i<-annot%>%split(f = annot$Sample)%>%lapply(function(x) str_remove(x$cell,".+_"))
i<-i[names(mats)]
mats<-lapply(names(i),function(x) mats[[x]][,colnames(mats[[x]])%in%i[[x]]])#filter myeloid cells
names(mats)<-names(i)
gc()

####################################INTEGRATING SAMPLES####################################

sc<-lapply(mats,CreateSeuratObject)
rm(mats);gc()
inte<-merge(sc[[1]],y = sc[names(i)[-1]],add.cell.ids =names(i))
rm(sc);gc()
annot<-annot%>%mutate(cell=paste0(Sample,'_',cell))
annot<-annot%>%filter(cell%in%colnames(inte))
annot<-annot[order(match(annot$cell,colnames(inte))),]
i<-colSums(is.na(annot))
annot<-annot[,i<nrow(annot)]
inte<-AddMetaData(inte,metadata = annot)
i<-rownames(inte)
i<-i[str_detect(i,pattern = "fung|bacte",negate = T)]
inte<-subset(inte,features = i)
rm(annot)
gc()


inte<-NormalizeData(inte,normalization.method = "LogNormalize")
inte<-ScaleData(inte,features = i)
gc()

inte<-RunPCA(inte,features = rownames(inte))#is not using all features
inte <- IntegrateLayers(inte, method=HarmonyIntegration,
                        new.reduction = "harmony",features = i)
inte[["RNA"]]<-JoinLayers(inte[["RNA"]])
#inte<-RunUMAP(inte,dims = 1:15)
inte<-RunUMAP(inte,reduction = "harmony",dims = 1:15,min.dist = 0.01,repulsion.strength = 5)
gc()

####################################CLUSTERING####################################

inte<-FindNeighbors(inte,reduction = "harmony",k.param = 50)
inte<-FindClusters(inte,algorithm = 4,resolution = .1)
cls<-as.numeric(unique(inte$RNA_snn_res.0.1))
inte$clus<-inte$RNA_snn_res.0.1
for (x in cls) {#overcluster?
  inte<-FindSubCluster(inte,graph.name = "RNA_snn",algorithm = 4,cluster = x,resolution = 1)
  inte@meta.data<-inte@meta.data%>%
    mutate(clus=ifelse(str_detect(sub.cluster,"_"),sub.cluster,clus))
}

####################################CELL TYPES####################################

bulkref<-readRDS("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/bulkref.RDS")
bulkref<-subset(bulkref,subset=label!="Erythrocytes")
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(sublabel=str_remove_all(label," \\(.+|Early ")%>%
           str_replace("emmo","emo")%>%str_replace("Th.+","CD4 T-cells"))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(toplabel=case_when(str_detect(label,"Smoo|Fibr|Mesen|GAS")~"Mesenchymal",
                            str_detect(label,"Th")~"CD4 T-cells",
                            str_detect(label,"CD8")~"CD8 T-cells",
                            str_detect(label,"Gamm|NK")~"Gamma_NK",
                            str_detect(label,"B-|Plas")~"B_Plasma",
                            str_detect(label,"EnEp|MUC5B|SOX9|IHH|Cili|Gland|KRT|Meso")~"Epithelial",
                            str_detect(label,"Macro|Dend")~"Myeloid",
                            TRUE~label))

bulkquery<-AggregateExpression(inte,group.by = "clus",return.seurat = T)
bulkquery<-ScaleData(bulkquery,assay = "RNA",features = i)
cormat<-cor(bulkquery@assays$RNA$scale.data,method = "spearman")
cormat<-1-cormat
bq_clus<-hclust(as.dist(cormat),method = "ward.D2")
res<-as.data.frame(cbind(clus=bq_clus$labels,
                         alt_clus=dynamicTreeCut::cutreeDynamic(dendro = bq_clus,
                                                                minClusterSize = 2,distM = cormat)))
res<-res%>%mutate(clus=str_remove(clus,"g")%>%str_replace("-","_"))
res<-inte@meta.data%>%right_join(res,by="clus")
res<-res[order(match(res$cell,colnames(inte))),]
inte<-AddMetaData(inte,metadata = res,col.name = "alt_clus")
get_labels<-function(cluster_column){
  bulkquery<-AggregateExpression(inte,group.by = cluster_column,return.seurat = T)
  bulkquery<-ScaleData(bulkquery,assay = "RNA",features = i)
  res <- SingleR(test=bulkquery@assays$RNA$scale.data,#this
                 ref=bulkref@assays$SCT$data, 
                 labels=bulkref@meta.data$sublabel, method = "cosine",
                 quantile = 0.999,tune.thresh = .1)
  res<-as.data.frame(res)%>%rownames_to_column(cluster_column)
  res<-res%>%mutate(per_intecl_sublabel=str_replace_na(pruned.labels,"unknown"))
  if(sum(str_detect(res[,cluster_column],"^g"))>0){
    res[,cluster_column]=str_remove(res[,cluster_column],"g")%>%str_replace("-","_")
  }
  res<-inte@meta.data%>%right_join(res,by = cluster_column)
  res<-res[order(match(res$cell,colnames(inte))),]
  return(res)}
labels_percl<-get_labels("alt_clus")
inte<-AddMetaData(inte,metadata = labels_percl,col.name = "per_intecl_sublabel")
labels_percl<-get_labels("clus")
inte@meta.data$per_intesubcl_sublabel<-labels_percl$per_intecl_sublabel.y
####################################OUTPUT####################################
inte@meta.data<-inte@meta.data%>%mutate(temp=paste(alt_clus,per_intecl_sublabel,sep = '-'))
p1<-DotPlot(inte,features = j,group.by ="temp" )+RotatedAxis()+
  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(size=8))
inte@meta.data<-inte@meta.data%>%mutate(temp=paste(alt_clus,clus,per_intesubcl_sublabel,sep = '-'))
p2<-DotPlot(inte,features = j,group.by ="temp" )+RotatedAxis()+
  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(size=8))
p3<-ggdendro::ggdendrogram(bq_clus)    
p4<-inte@meta.data%>%dplyr::count(clus,alt_clus)%>%
  arrange(alt_clus)%>%mutate(clus=factor(clus,levels=clus))%>%
  ggplot(aes(x=n,y=clus,fill=alt_clus))+geom_col()

png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",
                  str_replace_all(ct,"/| ","_"),"_labels.png"),
    width=1300,height=800)
plot_grid(p1,p2,plot_grid(p3,p4),nrow=3,rel_heights = c(.25,.5,.25))
dev.off()
#temp<-inte@reductions$umap@cell.embeddings%>%bind_cols(temp)
inte@meta.data%>%select(cell,alt_clus,per_intecl_sublabel,clus,per_intesubcl_sublabel)%>%
  data.table::fwrite(file=paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",
                                 str_replace_all(ct,"/| ","_"),"_labels.gz"))


