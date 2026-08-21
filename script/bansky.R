library(Seurat)
library(tidyverse)
library(SpatialExperiment)
library(scuttle)
library(SingleR)
library(BiocParallel)
slide<- commandArgs(trailingOnly=TRUE)
#slide<-"/media/Xenium_On_NAS/20260310__213926__031026_Endometrioma_Batch13/output-XETG00426__0081712__EDV086__20260310__213950"#v4 & myeloid issue

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
xenium.obj<-subset(xenium.obj,subset=nFeature_Xenium>4)
xenium.obj$nCount_Xenium<-colSums(xenium.obj@assays$counts$counts)
model <- lm(log10(nCount_Xenium)~log10(nFeature_Xenium), data = xenium.obj@meta.data)
xenium.obj$residuals<-residuals(model)
thr<-mean(xenium.obj$residuals) + (3 * sd(xenium.obj$residuals))
xenium.obj<-subset(xenium.obj,
                   cells=which(xenium.obj$residuals<thr&xenium.obj$residuals>(-1*thr)))
gc()

se <- SpatialExperiment(assay = list(counts = xenium.obj@assays$counts$counts), spatialCoords = xenium.obj[["fov"]]$centroids@coords)
#rm(xenium.obj)
gc()
# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se)

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE)
lambda <- c(0.2)##lambda 0.8 for domains
set.seed(1000)
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, 
                            lambda = lambda, k_neighbors = 15,
                            resolution=c(.5,1,2))
gc()
se <-  Banksy::connectClusters(se)

temp<-se@colData%>%as.data.frame()%>%rownames_to_column("cell")
xenium.obj<-AddMetaData(xenium.obj,metadata = temp)
xenium.obj<-NormalizeData(xenium.obj,assay = "counts",scale.factor=1000,normalization.method = "LogNormalize")
bulkquery<-AggregateExpression(xenium.obj,group.by = c("clust_M1_lam0.2_k15_res0.5"),return.seurat = T)
bulkquery<-ScaleData(bulkquery,assay = "counts",features = i)

bulkref<-readRDS("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/bulkref.RDS")
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(sublabel=str_remove_all(label," \\(.+|Early |Resident ")%>%
           str_replace("emmo","emo")%>%str_replace("Effector M","M"))%>%
  mutate(sublabel=case_when(str_detect(sublabel,"Fibr|Mesen|GAS")~"Mesenchymal",
                            str_detect(sublabel,"EnEp|MUC5B|SOX9|IHH|Cili|Gland")~"EnEpi",
                            str_detect(sublabel,"Th")~"CD4 T-cells",
                            #str_detect(sublabel,"CD8")~"CD8 T-cells",
                            #str_detect(sublabel,"Gamm|NK")~"Gamma/NK",
                            str_detect(sublabel,"Macro")~"Macrophages",
                            TRUE~sublabel))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(toplabel=case_when(str_detect(sublabel,"Smoo|Mese")~"Mesenchymal",
                            str_detect(sublabel,"CD8")~"CD8 T-cells",
                            str_detect(sublabel,"Gamm|NK")~"Gamma_NK",
                            str_detect(sublabel,"B-|Plas")~"B_Plasma",
                            str_detect(sublabel,"EnEp|Meso")~"Epithelial",
                            str_detect(sublabel,"Macro|Dend")~"Myeloid",
                            TRUE~sublabel))

bp <- MulticoreParam(6)
set.seed(112358)
#start_time <- Sys.time()
res <- SingleR(test=bulkquery@assays$counts$scale.data,#this
               ref=bulkref@assays$SCT$data, 
               labels=bulkref$toplabel, method = "cosine",
               BPPARAM=bp,quantile = 0.999,tune.thresh = .1)
res<-as.data.frame(res)
temp<-res%>%rownames_to_column("clust_M1_lam0.2_k15_res0.5")%>%
  mutate(clust_M1_lam0.2_k15_res0.5=str_remove(clust_M1_lam0.2_k15_res0.5,"g"))%>%
  select(clust_M1_lam0.2_k15_res0.5,pruned.labels)%>%inner_join(temp)
temp<-temp%>%mutate(res0.5Label=str_replace_na(pruned.labels,"missing"))
temp<-temp[order(match(temp$cell,colnames(xenium.obj))),]
xenium.obj<-AddMetaData(xenium.obj,metadata = temp,col.name = "res0.5Label")

bulkquery<-AggregateExpression(xenium.obj,group.by = c("clust_M1_lam0.2_k15_res1"),return.seurat = T)
bulkquery<-ScaleData(bulkquery,assay = "counts",features = i)
res <- SingleR(test=bulkquery@assays$counts$scale.data,#this
               ref=bulkref@assays$SCT$data, 
               labels=bulkref$sublabel, method = "cosine",
               BPPARAM=bp,quantile = 0.999,tune.thresh = .1)
res<-as.data.frame(res)
temp<-res%>%rownames_to_column("clust_M1_lam0.2_k15_res1")%>%
  mutate(clust_M1_lam0.2_k15_res1=str_remove(clust_M1_lam0.2_k15_res1,"g"))%>%
  select(clust_M1_lam0.2_k15_res1,pruned.labels)%>%inner_join(temp)
temp<-temp%>%mutate(res1Label=str_replace_na(pruned.labels,"missing"))

bansky_pca<-reducedDim(se,"PCA_M1_lam0.2")
subg<-temp%>%split(f = temp$res0.5Label)%>%
  lapply(function(x) bansky_pca[x$cell,])
subg<-lapply(subg,FindNeighbors,k.param=5)
subg<-lapply(subg,function(x) x$snn)

labels<-temp%>%split(f = temp$res0.5Label)%>%
  lapply(function(x) x%>%select(cell,res1Label)%>%table()%>%Matrix::Matrix())
labels<-lapply(names(subg),function(x)
  subg[[x]]%*%labels[[x]])
labels<-lapply(labels,function(x) 
  as.data.frame(cbind(nei_label=colnames(x)[apply(x,1,which.max)],
                      cell=rownames(x))))%>%bind_rows()  
temp<-temp%>%left_join(labels)
temp<-temp%>%dplyr::count(clust_M1_lam0.2_k15_res1,nei_label)%>%
  group_by(clust_M1_lam0.2_k15_res1)%>%mutate(total=sum(n),prop=n/total)%>%
  slice_max(prop)%>%mutate(res1clLabel=ifelse(prop>=.7,nei_label,"mix"))%>%
  select(clust_M1_lam0.2_k15_res1,res1clLabel)%>%inner_join(temp)
temp<-temp%>%ungroup()%>%
  mutate(lab=ifelse(res1clLabel=="mix",nei_label,res1clLabel))%>%
  dplyr::count(clust_M1_lam0.2_k15_res0.5,lab)%>%
  group_by(clust_M1_lam0.2_k15_res0.5)%>%mutate(prop=n/sum(n))%>%
  slice_max(prop,n=2)%>%
  mutate(res0.5clLabel=ifelse(prop>.9&min(n)<100,lab,"mix"))%>%
  slice_max(prop)%>%distinct(clust_M1_lam0.2_k15_res0.5,res0.5clLabel)%>%
  right_join(temp)
temp<-temp%>%mutate(finaLabel=ifelse(res0.5clLabel=="mix",res1clLabel,res0.5clLabel))%>%
  mutate(finaLabel=ifelse(finaLabel=="mix",nei_label,finaLabel))

temp<-temp[order(match(temp$cell,colnames(xenium.obj))),]
xenium.obj<-AddMetaData(xenium.obj,metadata = temp,col.name = "finaLabel")
#finaLabel B cells, CD4, exhausted,... make no sense :(










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
png(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnsky_labels.png"),
    width=1200,height=400)
DotPlot(xenium.obj,group.by = "sublabel",features = j)+RotatedAxis()+
  scale_color_gradient2(low="#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(size=8))
dev.off()

temp%>%select(-cl)%>%
  data.table::fwrite(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnsky_labels.gz"))