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
xenium.obj$nCount_Xenium<-colSums(xenium.obj@assays$counts$counts)
xenium.obj$ratio<-(xenium.obj$nCount_Xenium+0.01)/(xenium.obj$nFeature_Xenium+0.01)
xenium.obj<-subset(xenium.obj,subset=nFeature_Xenium>4)
j<-xenium.obj@meta.data%>%slice_max(ratio,n=100)%>%rownames()
j<-which(!colnames(xenium.obj)%in%j)
xenium.obj<-subset(xenium.obj,cells=j)
gc()

se <- SpatialExperiment(assay = list(counts = xenium.obj@assays$counts$counts), spatialCoords = xenium.obj[["fov"]]$centroids@coords)
#rm(xenium.obj)
gc()
# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE)
lambda <- c(0.2,0.8)
set.seed(1000)
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = 0.2, resolution = 2)#lambda 0.8 for domains

temp<-se@colData%>%as.data.frame()%>%rownames_to_column("cell")
xenium.obj<-AddMetaData(xenium.obj,metadata = temp,col.name = "clust_M1_lam0.2_k50_res2")
xenium.obj<-NormalizeData(xenium.obj,assay = "counts",scale.factor=1000,normalization.method = "LogNormalize")
xenium.obj<-ScaleData(xenium.obj,assay = "counts",features = i)
bulkquery<-AggregateExpression(xenium.obj,group.by = c("clust_M1_lam0.2_k50_res2"),return.seurat = T)

bulkref<-readRDS("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/bulkref.RDS")
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(sublabel=str_remove_all(label," \\(.+|Resident |Early ")%>%
           str_replace("Effector M","M")%>%str_replace("emmo","emo"))%>%
  mutate(sublabel=case_when(str_detect(sublabel,"Fibr|Mesen|GAS")~"Mesenchymal",
                            str_detect(sublabel,"EnEp|MUC5B|SOX9|IHH|Cili|Gland")~"EnEpi",
                            str_detect(sublabel,"Th")~"CD4 T-cells",
                            TRUE~sublabel))
bulkref@meta.data<-bulkref@meta.data%>%
  mutate(toplabel=case_when(str_detect(sublabel,"Smoo|Mese")~"Mesenchymal",
                            str_detect(sublabel,"T-|NK")~"T_NK",
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
temp<-res%>%rownames_to_column("cl")%>%mutate(cl=str_remove(cl,"g"))%>%
  select(cl,pruned.labels)%>%inner_join(temp,
                                        by = c("cl"="clust_M1_lam0.2_k50_res2"))
temp<-temp%>%mutate(pruned.labels=str_replace_na(pruned.labels,"missing"))
#xenium.obj<-AddMetaData(xenium.obj,metadata = temp,col.name = "pruned.labels")


bansky_pca<-reducedDim(se,"PCA_M1_lam0.2")
subg<-temp%>%split(f = temp$pruned.labels)%>%
  lapply(function(x) bansky_pca[x$cell,])
subg<-lapply(subg,FindNeighbors,k.param=5)
subg<-lapply(subg,function(x) x$snn)
subc<-lapply(subg,function(x) 
  FindClusters(x,#initial.membership = as.numeric(as.factor(labels[[x]]$bynei)),
               algorithm = 4,resolution = 1))
subc<-lapply(names(subc),function(x) subc[[x]]%>%mutate(toplabel=x))%>%bind_rows()
subc<-subc%>%rownames_to_column("cell")%>%mutate(cl=paste0(toplabel,res.1))
subc<-subc[order(match(subc$cell,colnames(xenium.obj))),]
xenium.obj<-AddMetaData(xenium.obj,metadata = subc,col.name = "cl")
bulkquery<-AggregateExpression(xenium.obj,group.by = "cl",return.seurat = T)
res <- SingleR(test=bulkquery@assays$counts$scale.data,#this
               ref=bulkref@assays$SCT$data, 
               labels=bulkref$sublabel, method = "cosine",
               BPPARAM=bp,quantile = 0.999,tune.thresh = .1)
res<-as.data.frame(res)
temp<-res%>%rownames_to_column("cl")%>%mutate(sublabel=pruned.labels,
                                              cl=str_replace(cl,"-","_"))%>%
  select(cl,sublabel)%>%inner_join(subc)#
temp<-temp%>%mutate(sublabel=str_replace_na(sublabel,"missing"))
temp<-temp[order(match(temp$cell,colnames(xenium.obj))),]
xenium.obj<-AddMetaData(xenium.obj,metadata = temp,col.name = "sublabel")

temp%>%select(-cl)%>%
  data.table::fwrite(file = paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/Xenium_labels/",slide,"bnsky_labels.gz"))