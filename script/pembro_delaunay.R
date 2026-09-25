#!/usr/bin/env Rscript
library(tidyverse)
library(geometry)

sampl<-commandArgs(trailingOnly=TRUE)

coords<-data.table::fread("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_domains.gz")
coords<-coords%>%add_count(circle)%>%filter(n>500)
temp<-coords%>%split(f = coords$circle)

temp<-lapply(temp,function(x) 
  x%>%select(cell_id,x_centroid,y_centroid))
cell_index<-lapply(temp,function(x) 
  x%>%mutate(i=1:nrow(x))%>%select(cell_id,i)%>%deframe)
temp<-lapply(temp,function(x) x%>%column_to_rownames("cell_id"))
delaunay_triangles  <-lapply(temp,function(x) 
  geometry::delaunayn(p = x,options = "Pp"))
delaunay_edges <- lapply(delaunay_triangles,function(x)
  rbind(x[ ,c(1,2)],x[ ,c(1,3)],x[ ,c(2,3)])%>%
    unique()%>%
    as.data.frame()%>%
    setNames(c("from","to")))

delaunay_DT =lapply(names(delaunay_edges),function(x) 
  data_frame(from = names(cell_index[[x]][delaunay_edges[[x]]$from]),
             to = names(cell_index[[x]][delaunay_edges[[x]]$to]),
             xstart = temp[[x]][delaunay_edges[[x]]$from, "x_centroid"],
             ystart = temp[[x]][delaunay_edges[[x]]$from, "y_centroid"],
             xend = temp[[x]][delaunay_edges[[x]]$to, "x_centroid"],
             yend = temp[[x]][delaunay_edges[[x]]$to, "y_centroid"]))
delaunay_DT<-delaunay_DT%>%bind_rows()%>%
  mutate(xlen=xstart-xend,
         ylen=yend-ystart,
         distance=sqrt(xlen**2+ylen**2))
delaunay_DT%>%data.table::fwrite("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/edges.gz")

delaunay_DT<-coords%>%mutate(from=cell_id,domain_from=res.0.5)%>%select(from,domain_from)%>%inner_join(delaunay_DT)
delaunay_DT<-coords%>%mutate(to=cell_id,domain_to=res.0.5)%>%select(to,domain_to)%>%inner_join(delaunay_DT)

labels<-data.table::fread("/media/Lawrenson_Lab_NAS/uthscsa/group_data/Xenium_labels/pembro/Banksy_matrix/prote_l0.2_labels.gz")
labels<-labels%>%mutate(m1_cl=case_when(m1_cl%in%c(1:4,"Dendritic Cells")~"Macrophages M2",
                                        m1_cl==5~"Dendritic cells",TRUE~m1_cl))
labels<-labels%>%mutate(c7_cl=case_when(c7_cl==1~"Mesenchymal cells",
                                        c7_cl%in%c(2,4)~"Activated Fibro/CAF",
                                        c7_cl==3~"C7-Activated Fibro",
                                        c7_cl==5~"C7-Fibro",
                                        c7_cl==6~"Activated Fibro",TRUE~m1_cl))
labels<-labels%>%filter(str_detect(c7_cl,"disc",negate = T))%>%
  mutate(finalabel=str_remove(c7_cl,"Effector "))
delaunay_DT<-labels%>%mutate(from=cell,type_from=finalabel)%>%select(from,type_from)%>%inner_join(delaunay_DT)
delaunay_DT<-labels%>%mutate(to=cell,type_to=finalabel)%>%select(to,type_to)%>%inner_join(delaunay_DT)

temp<-delaunay_DT%>%filter(domain_to==domain_from)%>%count(type_to,type_from,domain_to)
temp<-temp%>%pivot_wider(names_from = type_to,values_from = n)
temp<-temp%>%mutate(across(`Activated Fibro`:`Smooth muscle cells`,scale))

library(ComplexHeatmap)
temp%>%filter(domain_to%in%c(7,10,14))%>%select(`Activated Fibro`:`Smooth muscle cells`)%>%
  Heatmap(name="z",row_labels = temp$type_from[temp$domain_to%in%c(7,10,14)],
          row_split = temp$domain_to[temp$domain_to%in%c(7,10,14)],
          column_title = "to",cluster_rows = F,cluster_columns = F,
          row_names_gp = gpar(fontsize=8),rect_gp = gpar(col = "white", lwd = 1))









