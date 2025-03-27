##同一家族不同aptamer流式图
#需要的八条
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(forcats)

SUM159<-readRDS("./fyb/cripr/mnra_grna_aptamer.rds")
DefaultAssay(object = SUM159) <- 'aptamer'
SUM159 <- NormalizeData(object = SUM159) %>% FindVariableFeatures() %>% ScaleData()
need_aptamer<-paste0(paste0("Apt",'-1-'),c(21,1,25,14,19,9,22,24))
apt_n<-FetchData(SUM159,need_aptamer,slot = "data")
apt_n$group<-SUM159$target
apt_n$group<-str_split(apt_n$group,"\\-",simplify = T)[,1]
apt_n<-apt_n[apt_n$group=="PTK7"|apt_n$group=="Control",]

#转为长数据
apt_data <- melt(apt_n)

need_aptamer<-colnames(apt_n)[1:(ncol(apt_n)-1)]

apt_sort<-rev(need_aptamer)
apt_data$variable<-factor(apt_data$variable,
                          levels = apt_sort)
apt_data$group<-ifelse(grepl("PTK7",apt_data$group),"PTK7 sg cell",'Control cell')


#计算均值
apt_data %>% 
  group_by(group,variable) %>% 
  summarise(mean_value=mean(value)) %>% 
  ungroup() %>% 
  mutate(variable = fct_relevel(variable,
                                apt_sort)) %>% 
  mutate(new_col03=as.numeric(variable)) -> new.apt_data


#对PTK7组归一化
new_f<-new.apt_data[grepl("PTK7",new.apt_data$group),c(1,2,3)]
diff_val<-data.frame(variable=new_f$variable,value=new_f$mean_value-min(new_f$mean_value))

merged_df <- merge(diff_val,apt_data , by="variable")
merged_df$diff <- merged_df$value.y - merged_df$value.x

apt_data<-merged_df[,c(1,3,4,5)]
colnames(apt_data)[ncol(apt_data)]<-'value'
#新配色

rs<-rep('#9E9E9EFF',8)
rs_ptk7<-rep("#F08080",8)

# 创建组合水平的向量，确保顺序与颜色向量匹配  
combined_levels <- c(  
  paste("Control cell", levels(factor(apt_data$variable))),  
  paste("PTK7 sg cell", levels(factor(apt_data$variable)))  
)  

# 将 apt_data 中的 group 和 variable 组合成一个新列  
apt_data$combined_var <- with(apt_data, paste(group, variable))  

# 将 combined_var 转换为因子，并确保其水平顺序与上面创建的 combined_levels 一致  
apt_data$combined_var <- factor(apt_data$combined_var, levels = combined_levels)  

# 合并颜色向量  
all_colors <- c(rs, rs_ptk7)  

# 确保颜色向量的长度与 combined_var 的水平数相匹配  
stopifnot(length(all_colors) == length(levels(apt_data$combined_var)))  

p4.1<-ggplot(data=apt_data,aes(x=value,y=variable,fill=combined_var))+
  geom_density_ridges(
    bandwidth=0.04,
    quantile_lines=TRUE, 
    quantile_fun=function(x,...)mean(x),
    #linetype="dashed",
    scale=1,
    vline_linetype="dashed")+
  scale_fill_manual(values = c(rs,rs_ptk7), labels = c("Conrtrol cell", "PTK7 sg cell"))+#
  theme_classic()+
  scale_color_manual(values = c(rs,rs_ptk7), labels = c("Conrtrol cell", "PTK7 sg cell"))+
  theme_classic() +
  guides(fill="none",color="none")#rel_min_height=0.005


#不加图列
p4<-p4.1+
  coord_cartesian(clip="on")+
  labs(title = "", y = "",x="") +
  theme(text = element_text(face = "bold",size = 25,colour = 'black'),
        axis.text.y = element_text(size = 25,face = "bold",colour = 'black'),
        plot.title = element_text(size = 27, face = "bold",colour = 'black'))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0.1, 6), breaks = seq(0.1, 6, 1))+ 
  theme(axis.text.x = element_blank(),#去x轴
        axis.ticks = element_blank(),#去x轴刻度
        axis.text.y = element_text(size = 25)#y轴字体
  )+
  coord_cartesian(ylim = c(1.57,8.5))+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        plot.margin = margin(t = 0, r = 0, b = -10, l = 0, 
                             unit = "pt")
  )


ggsave(p4,filename ="fig5/家族1_8个aptamer_nc.pdf",width = 5,height = 6)
#write.csv(apt_n,"fig5/八条序列原始值.csv")
#write.csv(apt_data,"fig5/照PTK7组均值归一化（图数据）.csv")
