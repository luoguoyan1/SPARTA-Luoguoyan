#火山图加大最富集字体
#新火山图(中位值)
library(Seurat)
library(limma)
SUM159<-readRDS("D:/单细胞测序/第三批/SUM159_filter_mt_aptamer_100_nc473_scree_0.7_rename_apt_motif.rds")

SUM159$group<-ifelse(grepl("NC",SUM159$group),'Control',str_split(SUM159$group,"\\-",simplify = T)[,1])
protein_colors <- c("CDCP1" = "#53ABD8",
                    "ITGA3" = "#FF8C00",
                    "ITGB1" = "#AB82FF",
                    "NRP1" = "#7CCD7C",
                    "NRP2" = "#DEB887",
                    "PTK7" = "#F08080",
                    'PTPRF'='#B291B5',
                    'PTPRD'='#BB6B84')
SUM159_fi_sg<-subset(SUM159,target!="ITGB1-2")
SUM159_fi_sg<-subset(SUM159_fi_sg,target!="CD151-3")
SUM159_fi_sg<-subset(SUM159_fi_sg,target!="PTPRF-3")
SUM159_fi_sg<-subset(SUM159_fi_sg,target!='ITGA3-2')
SUM159_fi_sg<-subset(SUM159_fi_sg,target!='CDCP1-1')
SUM159_fi_sg<-subset(SUM159_fi_sg,target!='NRP2-2')
#SUM159_fi_sg<-subset(SUM159_fi_sg,target!='NRP1-3')

apt_exp <- as.data.frame(t(as.matrix(SUM159_fi_sg@assays$aptamer@counts)))


need_motif<-paste0("Clust-",c(1,2,3,4,5,6,7,11,13,14,15,17))
DefaultAssay(SUM159_fi_sg)<-'motif1w'
rownames(SUM159_fi_sg@assays$motif1w@counts)<-gsub('\\.','\\-',rownames(SUM159_fi_sg@assays$motif1w@counts))
data_motif<-FetchData(SUM159_fi_sg,need_motif,slot = 'counts')

data_motif$group<-SUM159_fi_sg$group
table(data_motif$group)
new_apt<-read.csv("C:/Users/fuybo/Desktop/new_apt_motif_name.csv")
clust_sum<-aggregate(.~motif,new_apt[,c(5,9)],sum)
new_apt$clus<-clust_sum$count[match(new_apt$motif,clust_sum$motif)]
new_apt$apt_ratio<-new_apt$count/new_apt$clus

#都去掉对应target，除ptk7-600细胞，其余150细胞，其中ptk7采用clust1，nrp2使用clust14，其余都用
for (i in unique(names(protein_clust_map))){
  message("Processing ", i)
  print(i)
  i<-'PTPRD'
  ko_name<-i
  
  data_motif_f<-data_motif[data_motif$group %in% c(i,'Control'),c('Clust-5','Clust-13','Clust-17','group')]
  data_motif_f<-data_motif[data_motif$group %in% c(i,'Control'),c('Clust-2','Clust-6','Clust-7','group')]
  data_motif_f<-data_motif[data_motif$group %in% c(i,'Control'),c('Clust-11','Clust-14','group')]
  data_motif_f<-data_motif[data_motif$group %in% c(i,'Control'),c('Clust-3','group')]
  data_motif_f<-data_motif[data_motif$group %in% c(i,'Control'),c('Clust-4','group')]
  data_motif_f<-data_motif[data_motif$group %in% c(i,'Control'),c('Clust-1','Clust-15','group')]
  
  data_motif_f$cell<-rownames(data_motif_f)
  # 定义函数来选择中位值附近的100个细胞
  select_near_median <- function(sub_df) {
    median_value <- median(sub_df$value)
    sub_df <- sub_df %>%
      mutate(dist_to_median = abs(value - median_value)) %>%
      arrange(dist_to_median) %>%
      head(150)
    sub_df$dist_to_median <- NULL  # 移除临时列
    return(sub_df)
  }
  
  # 按组进行过滤
  filtered_df<-data.frame()
  head(data_motif_f)
  for (t in colnames(data_motif_f)[1:(ncol(data_motif_f)-3)]){
    data_ff<-data_motif_f[,c(t,'group','cell')]
    colnames(data_ff)[1]<-'value'
    c <- data_ff %>%
      group_by(group) %>%
      group_modify(~ select_near_median(.x)) %>%
      ungroup()
    filtered_df<-rbind(filtered_df,c)
  }
  filtered_df<-filtered_df[!duplicated(filtered_df$cell),]
  dim(filtered_df)
  
  exp_unique<-apt_exp[filtered_df$cell,]
  exp_unique<-as.data.frame(t(exp_unique))
  exp_unique<-log2(exp_unique+1)
  
  df<-data.frame(df=filtered_df$group,cell=filtered_df$cell)
  
  group<-df
  
  
  design<-model.matrix(~0+factor(group$df))
  head(design)
  
  colnames(design)<-ifelse(grepl(i,colnames(design)),'ss','Control')
  
  rownames(design)<-rownames(group)
  head(design)
  contrast.matrix<-makeContrasts(ss- Control,levels=design)
  fit<-lmFit(exp_unique,design)
  fit2<-contrasts.fit(fit,contrast.matrix)
  fit2<-eBayes(fit2)
  options(digits=5)
  DEG<-topTable(fit2,coef=1,n=Inf)
  DEG$gene<-rownames(DEG)
  rownames(DEG)<-gsub("\u00B2",'',rownames(DEG))
  
  DEG['Apt-15-27',]
  DEG['Apt-11-15',]
  DEG['Apt-14-26',]
  DEG['Apt-3-3',]
  DEG['Apt-3-42',]
  DEG['Apt-4-5',]
  DEG<-DEG[grepl("-1-",rownames(DEG))|grepl("-15-",rownames(DEG)),]
  DEG<-DEG[grepl("-14-",rownames(DEG)),]
  DEG<-DEG[grepl("-5-",rownames(DEG)) | grepl("-13-",rownames(DEG)) | grepl("-17-",rownames(DEG)) ,]
  DEG<-DEG[grepl("-2-",rownames(DEG)) |grepl("-6-",rownames(DEG)) | grepl("-7-",rownames(DEG)),]
  DEG<-DEG[grepl("-11-",rownames(DEG)) | grepl("-14-",rownames(DEG)),]
  DEG<-DEG[grepl("-3-",rownames(DEG)) | grepl("-3-",rownames(DEG)),]
  DEG<-DEG[grepl("-4-",rownames(DEG)) | grepl("-4-",rownames(DEG)),]
  
  DEG<-DEG[order(DEG$logFC,decreasing = F),]
  write.csv(DEG,paste0("fig4火山图/",i,".csv"))
  DEG$gg<-ifelse(DEG$logFC>=0,"no_change",'change')
  DEG$gene<-rownames(DEG)
  
  
  DEG$label<-''
  DEG$group<-''
  DEG$label[seq(1,10,2)]<-rownames(DEG)[seq(1,10,2)]
  DEG$label[c(1,2,4,6,8,10)]<-rownames(DEG)[c(1,2,4,6,8,10)]
  DEG$label[c(1,3,7)]<-rownames(DEG)[c(1,3,7)]#nrp1
  DEG$label[c(1,5,7,9)]<-rownames(DEG)[c(1,5,7,9)]#itga3
  DEG$label[c(1,3,5,9)]<-rownames(DEG)[c(1,3,5,9)]#ptprd
  
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-9"))]<-paste0("Apt","-1-9")
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-21"))]<-paste0("Apt","-1-21")
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-22"))]<-paste0("Apt","-1-22")
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-24"))]<-paste0("Apt","-1-24")
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-1"))]<-paste0("Apt","-1-1 (29.78%)")
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-19"))]<-paste0("Apt","-1-19")
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-14"))]<-paste0("Apt","-1-14")
  DEG$label[which(rownames(DEG)==paste0("Apt","-1-25"))]<-paste0("Apt","-1-25")
  DEG$label[which(rownames(DEG)==paste0("Apt","-15-27"))]<-paste0("Apt","-15-27 (91.67%)")
  
  DEG$label[which(rownames(DEG)==paste0("Apt","-13-18"))]<-paste0("Apt","-13-18 (92.68%)")
  DEG$label[which(rownames(DEG)==paste0("Apt","-5-4"))]<-paste0("Apt","-5-4 (93.71%)")
  DEG$label[which(rownames(DEG)==paste0("Apt","-17-39"))]<-paste0("Apt","-17-39 (77.07%)")
  
  DEG$label[which(rownames(DEG)==paste0("Apt","-2-2"))]<-paste0("Apt","-2-2 (80.96%)")
  DEG$label[which(rownames(DEG)==paste0("Apt","-6-6"))]<-paste0("Apt","-6-6 (93.39)")
  DEG$label[which(rownames(DEG)==paste0("Apt","-7-10"))]<-paste0("Apt","-7-10 (93.11%)")
  
  DEG$label[which(rownames(DEG)==paste0("Apt","-11-15"))]<-paste0("Apt","-11-15 (90.78%)")
  DEG$label[which(rownames(DEG)==paste0("Apt","-14-26"))]<-paste0("Apt","-14-26 (93.31%)")
  
  DEG$label[which(rownames(DEG)==paste0("Apt","-3-3"))]<-paste0("Apt","-3-3 (90.04%)")
  DEG$label[which(rownames(DEG)==paste0("Apt","-3-42"))]<-paste0("Apt","-3-42")
  
  DEG$label[which(rownames(DEG)==paste0("Apt","-4-5"))]<-paste0("Apt","-4-5 (86.25%)")
  
  p3<-ggplot(data=DEG,aes(x=logFC,y=-log10(DEG$P.Value),color=gg))+geom_point()+
    scale_color_manual(values=c(as.character(protein_colors[ko_name]),"#696969"))
  data<-DEG#[grepl('A',DEG$label),]
  p3
  colo<-rep(as.character(protein_colors[ko_name]),nrow(data))
  colo[which(rownames(data)==c('Apt-1-1','Apt-15-27'))]<-"#ec3333"#PTK7
  
  colo<-rep(as.character(protein_colors[ko_name]),nrow(data))
  colo[rownames(data) %in% c("Apt-13-18",'Apt-5-4','Apt-17-39')]<-"#28769c"#cdcp1 "#2ca0d9"
  
  colo<-rep(as.character(protein_colors[ko_name]),nrow(data))
  colo[rownames(data) %in% c("Apt-2-2",'Apt-6-6','Apt-7-10')]<-"#5d975d"#nrp1
  
  colo<-rep(as.character(protein_colors[ko_name]),nrow(data))
  colo[rownames(data) %in% c("Apt-11-15",'Apt-14-26')]<-"#bb9c74"#nrp2
  
  colo<-rep(as.character(protein_colors[ko_name]),nrow(data))
  colo[rownames(data) %in% c("Apt-3-3",'Apt-3-42')]<-"#dc7a02"#itga3
  
  colo<-rep(as.character(protein_colors[ko_name]),nrow(data))
  colo[rownames(data) %in% c("Apt-4-5")]<-"#b169b6"#ptprf
  
  p3<-p3+geom_text_repel(data = data, aes(x = data$logFC,
                                          y = -log10(data$P.Value),
                                          label = label),
                         size = 5,box.padding = unit(0.5, "lines"),#6
                         point.padding = unit(0.8, "lines"),
                         colour = colo,
                         max.overlaps = 4000,
                         show.legend = FALSE)+ theme_bw () + 
    theme (panel.grid=element_blank ())
  p3
  
  p3<-p3+labs(title =paste0(ko_name," sgRNA vs Control"), y = "log10(P.Value)",x="logFC")+
    theme(panel.background = element_blank(),
          text = element_text(size = 12,family = 'sans'),
          axis.title = element_text(size = 12,family = 'sans'),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(size=12,angle = 0,vjust = 0.85,hjust = 0.75))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size = 12,family = 'sans'),
          axis.title.x = element_text(size = 12,family = 'sans'),
          axis.title.y = element_text(size = 12,family = 'sans'))+
    NoLegend()
  p3
  p3<-p3+theme(legend.position = c(0.95, 0.95),
               legend.justification = c(1,1))+
    theme(legend.title = element_blank())+
    NoLegend()
  p3
  
  ggsave(p3,filename = file.path("fig4火山图/", paste0(ko_name, "_150_.pdf")),width = 5,height = 5)
}