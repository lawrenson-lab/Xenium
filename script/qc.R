#!/usr/bin/env Rscript
library(tidyverse)
library(ggrepel)

#folder=commandArgs(trailingOnly=TRUE)

files<-list.files("/media/Xenium_On_NAS",full.names = T)
files<-files[str_detect(files,"atch")]
files<-unlist(lapply(files,list.files,full.names=T))
files<-sapply(files,list.files,full.names=T)

metrics_files<-files[str_detect(files,"metrics_summary.csv")]
metrics<-lapply(metrics_files,read_csv)%>%bind_rows()
i<-which(!is.na(apply(metrics,2,function(x) min(as.numeric(x)))))
plots<-lapply(names(i),function(v)
  metrics%>%ggplot(aes(x=(!!sym(v))))+geom_histogram()+
    ggtitle(paste("total runs=",sum(metrics[[v]]!=0),
                  ", avg=",round(mean(as.numeric(metrics[[v]])),4))))
ggsave(filename = "qc_metrics.png",width = 20,height = 20,
       plot = gridExtra::grid.arrange(grobs = plots, ncol = 5))

avgs<-colMeans(metrics[,i])
sds<-apply(metrics[,i],2,sd)       
temp<-lapply(names(i), function(x) 
  as.data.frame(cbind(var=paste0("Posible ",x," outlier"),
                      file=metrics_files[metrics[[x]]>avgs[x]+3*sds[x]|metrics[[x]]<avgs[x]-3*sds[x]])))%>%
  bind_rows()
temp<-temp%>%mutate(file=str_remove(file,"/metri.+"))%>%drop_na()
write_excel_csv(temp,"samples_to_verify.csv")

transcr_metrics_files<-files[str_detect(files,"transcripts.parquet")]
summarize_transcripts<-function(tmf){
  parquet_file <- arrow::read_parquet(tmf, as_data_frame = T)
  parquet_file<-parquet_file%>%filter(is_gene==T)
  gc()
  parquet_file<-parquet_file%>%mutate(extracellular=cell_id=="UNASSIGNED")
  parquet_file<-parquet_file%>%group_by(feature_name)%>%
    summarise(extracellular=sum(cell_id=="UNASSIGNED"),
              intracellular=sum(cell_id!="UNASSIGNED"),
              nuclear=sum(overlaps_nucleus==1),
              avg_qv=mean(qv),
              avg_nucleus_distance=mean(nucleus_distance))
  gc()
  parquet_file$id<-str_replace(tmf,".+Ba","Ba")%>%
    str_replace("\\/.+EDV","\\/EDV")%>%
    str_remove("__20.+")
  return(parquet_file)
}
transcr_metrics<-lapply(transcr_metrics_files,summarize_transcripts)%>%bind_rows()
transcr_metrics<-transcr_metrics%>%mutate(nuclear_pct=100*nuclear/intracellular)
transcr_metrics<-transcr_metrics%>%pivot_longer(-c(feature_name,id),names_to = "metric")
pdf(file = "transcript_qc_metrics.pdf",width = 30,height = 10)
transcr_metrics%>%filter(metric=="avg_qv")%>%ggplot(aes(x=id,y=value))+geom_boxplot()+
  geom_text_repel(data=transcr_metrics%>%group_by(id)%>%filter(metric=="avg_qv")%>%slice_min(value),
                  aes(x=id,y=value,label = feature_name),max.overlaps = 100)+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1))+ggtitle("avg_qv")
i<-unique(transcr_metrics$metric)
i<-i[str_detect(i,"qv",negate = T)]
lapply(i,function(x)
  transcr_metrics%>%filter(metric==x)%>%
    ggplot(aes(x=id,y=value))+geom_boxplot()+
    geom_text_repel(data=transcr_metrics%>%filter(metric==x)%>%group_by(id)%>%slice_max(value),
                    aes(x=id,y=value,label = feature_name),max.overlaps = 100)+
    geom_text_repel(data=transcr_metrics%>%filter(metric==x)%>%group_by(id)%>%slice_min(value),
                    aes(x=id,y=value,label = feature_name),max.overlaps = 100)+
    scale_y_log10()+theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1))+ggtitle(x))
dev.off()


cell_files<-files[str_detect(files,"cells.csv")]
cell_metrics<-lapply(cell_files,read_csv)
names(cell_metrics)<-cell_files%>%
  str_remove_all(".+output")%>%
  str_split(pattern = "__")%>%
  lapply(function(x) paste(x[3:2],collapse = "_"))
cell_metrics<-lapply(names(cell_metrics),function(x) 
  cell_metrics[[x]]%>%mutate(id=x))%>%bind_rows
i<-colnames(cell_metrics)[4:13]
cell_metrics<-cell_metrics%>%pivot_longer(i,names_to = "metric")
pdf(file = "cell_qc_metrics.pdf",width = 30,height = 10)
lapply(i,function(x)
  cell_metrics%>%filter(metric==x)%>%
    ggplot(aes(x=id,y=value))+geom_boxplot()+
    theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1))+ggtitle(x))
dev.off()
