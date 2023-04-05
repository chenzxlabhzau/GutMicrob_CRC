setwd("F:/A_MetaboAnalysis/FMT/")
library(limma)
library(reshape2)
library(dplyr)
FMT_species <- read.delim("F:/A_MetaboAnalysis/FMT/merged_abundance_table_species.txt",row.names = 1 ,header=T)
head(FMT_species)

#View(FMT_species)
FMT_species <- FMT_species[,colnames(FMT_species)[order(colnames(FMT_species))]]
colnames(FMT_species) <- gsub("F_F_","F_FMT_AF_",colnames(FMT_species))%>%
  gsub("F_M_","F_FMT_AM_",.)%>%gsub("M_F_","M_FMT_AF_",.) %>%gsub("M_M_","M_FMT_AM_",.)

group <- substr(colnames(FMT_species),1,8)
FMT_species$sum <- rowSums(FMT_species)
FMT_species_order <- FMT_species[order(FMT_species$sum, decreasing=TRUE), ]
colnames(FMT_species)
FMT_species_order$sum <- NULL
#colnames(test) <- substr(colnames(test),3,3)
FMT_species$sum <- NULL
FMTInfo <- data.frame(group=substr(colnames(FMT_species),1,8))
rownames(FMTInfo) <- colnames(FMT_species)
FMTInfo$source <- substr(colnames(FMT_species),8,8)
FMTInfo$target <- substr(colnames(FMT_species),1,1)
dim(FMTInfo)
colnames(FMT_species_order)
colnames(FMT_species)
FMTInfo
mat1 <- FMT_species_order[1:20,]
mat2 <- FMT_species_order[1:30,]
dim(mat1)
colnames(mat1)
dim(FMTInfo)
tmp1 = apply(mat1,1,function(x){tapply(x,FMTInfo$group,median)}) %>% t() %>% as.data.frame()
tmp2 = apply(mat2,1,function(x){tapply(x,FMTInfo$group,median)}) %>% t() %>% as.data.frame()
tmp = apply(FMT_species_order,1,function(x){tapply(x,FMTInfo$group,median)}) %>% t() %>% as.data.frame()
n <- t(scale(t(mat2)))
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#4D7587","#85AEC9","#976276","#E1AB8D")),
                       labels = c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(FMTInfo$group,levels =c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"))
?Heatmap
colnames(n)
Heatmap(n,name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = T,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))

#colnames(FMT_species)
#tmp = apply(FMT_species[,1:36],1,function(x){tapply(x,FMTInfo$group, mean)}) %>% t() %>% as.data.frame()
#head(tmp)
#ggdata = lapply(tmp, function(x){rownames(tmp)[x != 0] })
#ggven <- data.frame(F_F=c(ggdata$F_F,rep("",22)),F_M=c(ggdata$F_M,rep("",10)),M_F=c(ggdata$M_F,rep("",13)),M_M=ggdata$M_M)
#write.table(ggven,"F:/A_MetaboAnalysis/ggven.txt", sep = '\t', col.names = NA, quote = FALSE)
#BiocManager::install("ggvenn")
#library(ggvenn)
#?ggvenn
#p=ggvenn(ggdata,       
#         show_percentage = F,show_elements=F,text_size = 8,
#         stroke_color = "white",
#         fill_color = c("#DA4F43","#FCC541","#4E8AF1","#0096A6"),
#         set_name_color = c("#DA4F43","#FCC541","#4E8AF1","#0096A6"))  #最后一个颜色要改
library(ggplot2)
#堆积图
##top20
#devtools::install_github("microbiota/amplicon")
#library(amplicon)
#library("tidyr")
library(tidyverse)
others <-apply(tmp[!rownames(tmp)%in%rownames(tmp1),],2,sum)%>%as.data.frame(.)
colnames(others) <- "Others"
others <- as.data.frame(t(others))
mat <- rbind(tmp1,others[,-5])
mat
ggdata = mat %>% rownames_to_column(var = "species") %>% gather(key = "group",value = "value", -"species")
unique(ggdata$species)
ggdata$species <- factor(ggdata$species,levels = unique(ggdata$species))
ggplot(ggdata) +				
  geom_bar(aes(x = group,y=value,fill = species),stat = "identity",position = "fill",colour= 'black') +						
  scale_fill_manual(values  = c("#6686c6","#a0c1db","gray","#E65F92","#c77364","#ce8f5c",						
                                "#7bac80","#75a5c0","#b5181a","#b72d9e",						
                                "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",						
                                "#3e6926","#0a0883","#49ddd0","#e0f8f7",						
                                "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+	
  xlab("Group")+ylab("% of species")+guides(fill = guide_legend( ncol = 1, byrow = TRUE))	
##################### Diff Microbiome Source M vs F ######################
colnames(FMT_species_order)
otu1<- FMT_species_order[,c(26:33,8:17,18:25,1:7)]
colnames(otu1)
p_value=c() 
for (i in 1:nrow(otu1)){
  data_new <- melt(otu1[i,])
  data_new$variable<-substr(colnames(otu1),8,8)
  P=wilcox.test(value~variable,data_new)
  p_value[i]=P$p.value
  p_value[i]<- round(p_value[i],3) #设置小数点位数
}
FC <- c()
for(i in 1:nrow(otu1)){
  FC[i]<-"NA"
  FC[i]<-mean(as.numeric(otu1[i,1:18]))/mean(as.numeric(otu1[i,19:33])) 
}
otu1$pvalue <- p_value
otu1$FC <- FC
otu1$log2FC <- log2(as.numeric(otu1$FC))
#View(otu1)
sourceMF<- subset(otu1,pvalue<0.05)%>%as.data.frame()
sourceMF <- sourceMF[order(sourceMF$log2FC),]
#View(sourceMF)
sourceMF$group <- 0
sourceMF$group[which(sourceMF$log2FC>0)]="Source_M"
sourceMF$group[which(sourceMF$log2FC<0)]="Source_F"
data <- otu1[rownames(otu1)%in%rownames(sourceMF),]
n <- t(scale(t(data[,1:33])))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#4D7587","#85AEC9","#976276","#E1AB8D")),
                       labels = c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(substr(colnames(otu1[,1:33]),1,8),levels =c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"))
?Heatmap
dim(n)
Heatmap(n[rownames(sourceMF),],name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = F,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))


###################### Species ##########################
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
# library(reshape)
colnames(FMT_species_order)
data <- FMT_species_order[rownames(FMT_species_order)%in%c("Akkermansia_muciniphila",
                                                           "Parabacteroides_goldsteinii","Dubosiella_newyorkensis","Bifidobacterium_pseudolongum"),]
data$species <- rownames(data)
ggdata <- melt(data)
ggdata$group <-substr(ggdata$variable,1,8)%>%gsub("_","-",.)%>%
  factor(.,levels = c("M-FMT-AM","F-FMT-AM","M-FMT-AF","F-FMT-AF"))
ggdata$source <- paste0("Source","-",substr(ggdata$variable,8,8))%>%
  factor(.,levels = c("Source-M","Source-F"))
data1 <- ggdata[,c(1,3,4)]
colnames(data1) <- c("species","value","group") 
data2 <- ggdata[,c(1,3,5)]
colnames(data2) <- c("species","value","group")
data <- rbind(data1,data2)
data$group <- factor(data$group,levels = c("Source-M","Source-F","M-FMT-AM","F-FMT-AM","M-FMT-AF","F-FMT-AF"))

ggplot(data[data$species=="Akkermansia_muciniphila",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+
  theme_classic()+
  scale_fill_manual(values=c("#85AEC9","#E1AB8D","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  geom_signif(comparisons = list(c("Source-M","Source-F"),c("M-FMT-AM","M-FMT-AF"),c("F-FMT-AM","F-FMT-AF"),
                                 c("M-FMT-AM","F-FMT-AM"),c("M-FMT-AF", "F-FMT-AF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative abundance")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+ggtitle("A.muciniphila")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(data[data$species=="Parabacteroides_goldsteinii",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+
  theme_classic()+
  scale_fill_manual(values=c("#85AEC9","#E1AB8D","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  geom_signif(comparisons = list(c("Source-M","Source-F"),c("M-FMT-AM","M-FMT-AF"),c("F-FMT-AM","F-FMT-AF"),
                                 c("M-FMT-AM","F-FMT-AM"),c("M-FMT-AF", "F-FMT-AF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative abundance")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+ggtitle("P.goldsteinii")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')



ggplot(ggdata[ggdata$species=="Parabacteroides_goldsteinii",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+
  theme_classic()+
  scale_fill_manual(values=c("#4D7587","#85AEC9","#976276","#E1AB8D"))+
  geom_signif(comparisons = list(c("M-FMT-AM","M-FMT-AF"),c("F-FMT-AM","F-FMT-AF"),
                                 c("M-FMT-AM","F-FMT-AM"),c("M-FMT-AF", "F-FMT-AF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative abundance")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+ggtitle("P.goldsteinii")+ 
  theme(plot.title = element_text(size = 24, face = "italic"))+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')



ggplot(ggdata[ggdata$species=="Dubosiella_newyorkensis",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+
  theme_classic()+
  scale_fill_manual(values=c("#4D7587","#85AEC9","#976276","#E1AB8D"))+
  geom_signif(comparisons = list(c("M-FMT-AM","M-FMT-AF"),c("F-FMT-AM","F-FMT-AF"),
                                 c("M-FMT-AM","F-FMT-AM"),c("M-FMT-AF", "F-FMT-AF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative abundance")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+ggtitle("D.newyorkensis")+ 
  theme(plot.title = element_text(size = 24, face = "italic"))+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')




##################################### M→M vs F→M #######################################
colnames(FMT_species_order)
otu1<- FMT_species_order[,c(26:33,18:25)]
colnames(otu1)
p_value=c() 
for (i in 1:nrow(otu1)){
  data_new <- melt(otu1[i,])
  data_new$variable<-substr(colnames(otu1),8,8)
  P=wilcox.test(value~variable,data_new)
  p_value[i]=P$p.value
  p_value[i]<- round(p_value[i],3) #设置小数点位数
}
FC <- c()
for(i in 1:nrow(otu1)){
  FC[i]<-"NA"
  FC[i]<-mean(as.numeric(otu1[i,1:8]))/mean(as.numeric(otu1[i,9:16])) 
}
otu1$pvalue <- p_value
otu1$FC <- FC
otu1$log2FC <- log2(as.numeric(otu1$FC))
#View(otu1)
MM_FM<- subset(otu1,pvalue<0.05)%>%as.data.frame()
MM_FM <- MM_FM[order(MM_FM$log2FC),]
#View(sourceMF)
MM_FM$group <- 0
MM_FM$group[which(MM_FM$log2FC>0)]="Source_M"
MM_FM$group[which(MM_FM$log2FC<0)]="Source_F"
#View(sourceMF)
intersect(rownames(sourceMF[sourceMF$group=="Source_M",]),rownames(MM_FM[MM_FM$group=="Source_M",]))
intersect(rownames(sourceMF[sourceMF$group=="Source_F",]),rownames(MM_FM[MM_FM$group=="Souece_F",]))

data <- FMT_species_order[rownames(FMT_species_order)%in%rownames(MM_FM),]
n <- t(scale(t(data[,1:33])))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#4D7587","#85AEC9","#976276","#E1AB8D")),
                       labels = c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(substr(colnames(FMT_species_order[,1:33]),1,8),levels =c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"))
?Heatmap
dim(n)
Heatmap(n[rownames(MM_FM),],name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = F,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))















##################################### M→F vs F→F #######################################
colnames(FMT_species_order)
otu1<- FMT_species_order[,c(1:17)]
colnames(otu1)
p_value=c() 
for (i in 1:nrow(otu1)){
  data_new <- melt(otu1[i,])
  data_new$variable<-substr(colnames(otu1),8,8)
  P=wilcox.test(value~variable,data_new)
  p_value[i]=P$p.value
  p_value[i]<- round(p_value[i],3) #设置小数点位数
}
FC <- c()
for(i in 1:nrow(otu1)){
  FC[i]<-"NA"
  FC[i]<-mean(as.numeric(otu1[i,8:17]))/mean(as.numeric(otu1[i,1:7])) 
}
otu1$pvalue <- p_value
otu1$FC <- FC
otu1$log2FC <- log2(as.numeric(otu1$FC))
#View(otu1)
MF_FF<- subset(otu1,pvalue<0.05)%>%as.data.frame()
MF_FF <- MF_FF[order(MF_FF$log2FC),]
#View(sourceMF)
MF_FF$group <- 0
MF_FF$group[which(MF_FF$log2FC>0)]="Source_M"
MF_FF$group[which(MF_FF$log2FC<0)]="Source_F"


data <- FMT_species_order[rownames(FMT_species_order)%in%rownames(MF_FF),]
n <- t(scale(t(data[,1:33])))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#4D7587","#85AEC9","#976276","#E1AB8D")),
                       labels = c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(substr(colnames(FMT_species_order[,1:33]),1,8),levels =c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"))
?Heatmap
dim(n)
Heatmap(n[rownames(MF_FF),],name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = F,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))











ggdata = list(rownames(MM_FM),rownames(MF_FF),rownames(sourceMF))
names(ggdata) <- c("M-FMT(MvsF)","F-FMT(MvsF)","FMT-MvsFMT-F")
ggdata
#BiocManager::install("ggvenn")
library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 10,
         stroke_color = "black",
         fill_color = c("#85AEC9","#A0B3A1","#B2662A"),
         set_name_color = c("black","black","black")) 
p


a=intersect(rownames(sourceMF[sourceMF$group=="Source_M",]),rownames(MF_FF[MF_FF$group=="Source_M",]))
b=intersect(rownames(sourceMF[sourceMF$group=="Source_F",]),rownames(MF_FF[MF_FF$group=="Source_F",]))

c=intersect(rownames(sourceMF[sourceMF$group=="Source_M",]),rownames(MM_FM[MM_FM$group=="Source_M",]))
d=intersect(rownames(sourceMF[sourceMF$group=="Source_F",]),rownames(MM_FM[MM_FM$group=="Source_F",]))
e=intersect(rownames(sourceMF),c(rownames(MM_FM),rownames(MF_FF)))
e
dim(sourceMF)

data <- FMT_species_order[rownames(FMT_species_order)%in%e,]
n <- t(scale(t(data[,1:33])))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#4D7587","#85AEC9","#976276","#E1AB8D")),
                       labels = c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(substr(colnames(FMT_species_order[,1:33]),1,8),levels =c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"))
?Heatmap
n
Heatmap(n[e,],name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = T,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))

FMT_select_diff <- sourceMF[rownames(sourceMF)%in%e,]
View(FMT_select_diff)










data <- FMT_species_order[rownames(FMT_species_order)%in%rownames(MF_FF),]
n <- t(scale(t(data[,1:33])))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#4D7587","#85AEC9","#976276","#E1AB8D")),
                       labels = c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(substr(colnames(FMT_species_order[,1:33]),1,7),levels =c("M_FMT_AM","F_FMT_AM","M_FMT_AF","F_FMT_AF"))
?Heatmap
n
Heatmap(n[rownames(MF_FF),],name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = F,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))

write.table(MF_FF,"F:/A_MetaboAnalysis/FMT/MF_FF.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(MM_FM,"F:/A_MetaboAnalysis/FMT/MM_FM.txt", sep = '\t', col.names = NA, quote = FALSE)


##################### α diversity ##########################
library(vegan)
FMT_counts <- read.delim("F:/A_MetaboAnalysis/FMT/merged_abundance_counts_species.txt", row.names=1)
colnames(FMT_counts)
FMT_counts$F_F4_modify2_counts <- NULL
FMT_counts$M_M4_modify2_counts <- NULL
FMT_counts$F_F7_modify2_counts <- NULL
otu_count <- FMT_counts
colnames(otu_count)<- gsub("F_F","F-FMT-AF-",colnames(otu_count))%>%
  gsub("F_M","F-FMT-AM-",.)%>%gsub("M_F","M-FMT-AF-",.) %>%gsub("M_M","M-FMT-AM-",.)

alpha_diversity <- function(x, tree = NULL) {
  otus <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- vegan::diversity(x, index = 'shannon',base = 2)
  Simpson <- vegan::diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  
  
  result <- data.frame(otus, ACE,Chao1, Shannon, Simpson, goods_Coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(otus, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_Coverage)
  }
  
  
  result
}
library(reshape2)
otu_t <- t(otu_count)
alpha <- alpha_diversity(otu_t)
alpha$ACE<-NULL
alpha$goods_Coverage <- NULL
alpha$Simpson <- NULL
alpha$Chao1 <- NULL
alpha$group=substr(rownames(alpha),1,8)
alpha$source=paste0("Source","-",substr(rownames(alpha),8,8))
head(alpha)
otu_diversity_long <- melt(alpha,id.vars = c("source","group")) 
otu_diversity_long$value <- as.numeric(otu_diversity_long$value)

data1 <- otu_diversity_long[,c(1,3,4)]
colnames(data1) <- c("group","index","value") 
data2 <- otu_diversity_long[,c(2,3,4)]
colnames(data2) <- c("group","index","value")
data <- rbind(data1,data2)
data$group <- factor(data$group,levels = c("Source-M","Source-F","M-FMT-AM","F-FMT-AM","M-FMT-AF","F-FMT-AF"))

## violin
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

Data_summary <- summarySE(otu_diversity_long, measurevar="value", groupvars=c("group","variable"))

otu_diversity_long$group <- factor(otu_diversity_long$group,
                                   levels = c("Source-M","Source-F","M-FMT-AM","F-FMT-AM","M-FMT-AF","F-FMT-AF"))

################ otus #########################
ggplot(data[data$index=="otus",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+theme_classic()+
  scale_fill_manual(values=c("#85AEC9","#E1AB8D","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  geom_signif(comparisons = list(c("Source-M","Source-F"),c("M-FMT-AM","M-FMT-AF"),c("F-FMT-AM","F-FMT-AF"),
                                 c("M-FMT-AM","F-FMT-AM"),c("M-FMT-AF", "F-FMT-AF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')







################## shannon #################
ggplot(otu_diversity_long[otu_diversity_long$variable=="Shannon",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+theme_classic()+
  scale_fill_manual(values=c("#85AEC9","#E1AB8D","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  geom_signif(comparisons = list(c("M-FMT-AM","M-FMT-AF"),c("F-FMT-AM","F-FMT-AF"),
                                 c("M-FMT-AM","F-FMT-AM"),c("M-FMT-AF", "F-FMT-AF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(data[data$index=="Shannon",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+theme_classic()+
  scale_fill_manual(values=c("#85AEC9","#E1AB8D","#4D7587","#85AEC9","#976276","#E1AB8D"))+
  geom_signif(comparisons = list(c("Source-M","Source-F"),c("M-FMT-AM","M-FMT-AF"),c("F-FMT-AM","F-FMT-AF"),
                                 c("M-FMT-AM","F-FMT-AM"),c("M-FMT-AF", "F-FMT-AF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')




################### β diversity #########################
library(ape)
metaphlan_anno = read.table("F:/A_MetaboAnalysis/metaphlan_anno.txt",header = T) 
#View(metaphlan_anno)
metaphlan_tree = ape::read.tree("F:/A_MetaboAnalysis/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")
metaphlan_tree$tip.label <- metaphlan_tree$tip.label %>% gsub(".+\\|s__","",.)
my_tree = ape::keep.tip(metaphlan_tree,intersect(rownames(otu),metaphlan_tree$tip.label))
library(ggplot2)
library(ggdendro)
library(phyloseq)
library(tidyverse)
library(ggfun)
library(ggtree)
library(ggstance)
library(amplicon)
dim(otu_count)
colnames(otu_count) <- gsub("_modify2_counts","",colnames(otu_count))
colnames(FMT_species_order)
meta <- data.frame(group=substr(colnames(otu_count),1,8))
rownames(meta) <- colnames(otu_count)
meta$source <- substr(colnames(otu_count),8,8)
meta$target <- substr(colnames(otu_count),1,1)
dim(meta)
my_phy = phyloseq(otu_table(as.matrix(otu_count),taxa_are_rows = TRUE),
                   sample_data(meta),phy_tree(my_tree),tax_table(metaphlan_anno %>% as.matrix()))  #这里是否要as.matrix很严格
 
weight_unifrac= distance(my_phy,method = "wunifrac") %>% as.matrix()
unweight_unifrac = distance(my_phy,method = "unifrac") %>% as.matrix()
library(vegan)
library(ape)
res = pcoa(unweight_unifrac,correction = "none",rn = NULL)
biplot.pcoa(res)
res
pcoachoose = c("Axis.1","Axis.2")
ggdata = res$vectors[,pcoachoose] %>% cbind(.,meta)

ggdata$group <- factor(gsub("_","-",ggdata$group),levels = c("M-FMT-AM","F-FMT-AM","M-FMT-AF","F-FMT-AF"))
ggdata$source <- factor(ggdata$source,levels = c("M","F"))
ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),shape=group,fill=group))+						
  geom_point(size=8) +stat_ellipse(type = "t", linetype = 1, show.legend = FALSE,aes(color = group))+						
  ylab("Bray-Curtis diversity") +						
  theme(axis.title.x=element_blank())+						
  theme_classic()+theme(axis.title =element_text(size = 20,face = "bold"),						
                        axis.text =element_text(size = 20,face = "bold", color = 'black'))+						
  theme(axis.text.x = element_text(size = 20,face = "bold"))+						
  theme(axis.text.y = element_text(size = 20,face = "bold"))+						
  theme(legend.text = element_text(size=14,face = "bold"),legend.title  = element_blank(),legend.position = "top")+						
  scale_shape_manual(values = c(21,21,22,22))+						
  scale_color_manual(values=c("#4D7587","#85AEC9","#976276","#E1AB8D"))+						
  scale_fill_manual(values=c("#4D7587","#85AEC9","#976276","#E1AB8D"))+						
  xlab("PCoA1 (61.2%)")+ylab("PCoA2 (13.9%)")						

ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),shape=source,fill=source))+						
  geom_point(size=8) +stat_ellipse(type = "t", linetype = 1, show.legend = FALSE,aes(color = source))+						
  ylab("Bray-Curtis diversity")  +						
  theme(axis.title.x=element_blank())+						
  theme_classic()+theme(axis.title =element_text(size = 20,face = "bold"),						
                        axis.text =element_text(size = 20,face = "bold", color = 'black'))+						
  theme(axis.text.x = element_text(size = 20,face = "bold"))+						
  theme(axis.text.y = element_text(size = 20,face = "bold"))+						
  theme(legend.text = element_text(size=14,face = "bold"),legend.title  = element_blank(),legend.position = "top")+				
  scale_shape_manual(values = c(21,22))+						
  scale_color_manual(values=c("#85AEC9","#E1AB8D"))+						
  scale_fill_manual(values=c("#85AEC9","#E1AB8D"))+						
  xlab("PCoA1 (61.2%)")+ylab("PCoA2 (13.9%)")	

#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",force = TRUE)
library(pairwiseAdonis)

# This is a wrapper function for multilevel pairwise comparison 
# using adonis() from package 'vegan'. 
# The function returns adjusted p-values using p.adjust().
unweight_unifrac = distance(my_phy,method = "unifrac") %>% as.matrix()
pcoa.au <- cmdscale(unweight_unifrac, k=3, eig=T)
points.au <- as.data.frame(pcoa.au$points)
colnames(points.au) <- c("x", "y", "z") 
eig.au <- pcoa.au$eig
points.au <-  cbind(points.au, meta[match(rownames(points.au), rownames(meta)), ])
## 计算显著性：
dune.pairwise.adonis <- pairwise.adonis(unweight_unifrac, factors=factor(meta$group),
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)
dune.pairwise.adonis <- pairwise.adonis(unweight_unifrac, factors=factor(meta$source),
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)
write.table(dune.pairwise.adonis,"F:/A_MetaboAnalysis/FMT/FMT_adonis.txt", sep = '\t', col.names = NA, quote = FALSE)







library(vegan)
library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)
pcoa<- pcoa(unweight_unifrac, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]

plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,meta$source)
colnames(plotdata) <-c("sample","PC1","PC2","Group")
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)

plotdata$Group <- factor(plotdata$Group)

length=length(unique(as.character(group)))
times1=length%%8
times1
res1=length%%8
times2=length%%5
res2=length%%5
col1=rep(1:8,times1)
col=c(col1,1:res1)
pich1=rep(c(21:24),times2)
pich=c(pich1,15:(15+res2))
cbbPalette <- c("#DA4F43","#FECA41","#4E8AF1","#0096A6","#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",                "#ADD1E5")

yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~Group,data = plotdata)
plotdata$Group
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)
fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)
test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
test$Group <- factor(test$Group)

p1 <- ggplot(plotdata,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 7,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=20,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")
p1
p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 7,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=20,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")

p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=8,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=28),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=28),
        axis.title.y=element_text(colour='black', size=28),
        axis.text=element_text(colour='black',size=28),
        legend.title=element_text(size = 24,face = "bold"),
        legend.text=element_text(size=20),
        legend.key=element_blank(),legend.position = c(0.9,0.18),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1,"cm")) +
  guides(fill = guide_legend(ncol = 1))

p2
?adonis


otu.adonis=adonis(unweight_unifrac~Group,data =FMTInfo,distance = "bray")

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],  "\nR2 = ",round(otu.adonis$aov.tab$R2[1],4),  "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),
            size = 7) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)

p5
#weight
res = pcoa(weight_unifrac,correction = "none",rn = NULL)
biplot.pcoa(res)

pcoachoose = c("Axis.1","Axis.2")
ggdata = res$vectors[,pcoachoose] %>% cbind(.,sampleInfo)

p = ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),colour= Group))+geom_point(size= 4)+
  stat_ellipse(type = "t", linetype = 1)
mi=c("#DA4F43","#FCC541","#4E8AF1","#0096A6")
p=p+theme_classic()+scale_colour_manual(values = mi)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text = element_text(size = 14),axis.text.x = element_text(angle = 0, hjust =0.5),
        legend.text = element_text(size=12),legend.position = "top")+xlab("PCoA1")+ylab("PCoA2")


#NMDS
library("amplicon", lib.loc="~/R/win-library/4.0")
library(phyloseq)
?BetaDiv
result=BetaDiv(otu=otu, map=FMTInfo, group="Group", 
               dist="bray", method="NMDS", Micromet="adonis")
result
#提取排序散点图(结果列表中的1)
p=result[[1]]
p
#ggsave(paste0("p3.NMDS.bray.jpg"), p, width=89, height=56, units="mm")
#ggsave(paste0("p3.NMDS.bray.pdf"), p, width=89, height=56, units="mm")

# 提取出图坐标
plotdata=result[[2]]
plotdata[1:3,1:3]

# 提取带标签排序散点图
p=result[[3]]
p
ggsave(paste0("p4.NMDS.bray.label.jpg"), p, width=89, height=56, units="mm")
ggsave(paste0("p4.NMDS.bray.label.pdf"), p, width=89, height=56, units="mm")

# 提取两两比较差异检测结果
pair=result[[4]]
pair
# 提取全部组整体差异检测结果
Mtest=result[[5]]

#otu_t <- otu
nmds1 = metaMDS(bray,k = 2,try = 100)
nmds1$stress
nmds_dis_species <- wascores(nmds1$points, otu_t)
stressplot(nmds1, main = "Shepard")
ggdata = cbind(nmds_dis_species%>% as.data.frame() , sampleInfo)

dim(nmds1$points %>% as.data.frame())
dim(sampleInfo)
nmdschoose = c("MDS1","MDS2")
ggplot(data = ggdata, aes(MDS1, MDS2)) +
  geom_point(size=2,aes(color =  Group)) +
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +   #添加置信椭圆，注意不是聚类
  scale_color_manual(values =c("#DA4F43","#FCC541","#4E8AF1","#0096A6")) +
  scale_fill_manual(values = c("#DA4F43","#FCC541","#4E8AF1","#0096A6")) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5),legend.title = element_blank()) +
  #, legend.position = 'none'
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5)

ggplot(data = ggdata,aes(x=get(nmdschoose[1]),y=get(nmdschoose[2])))+
  geom_point(size= 4 ,aes(shape= group,fill= group)) +
  stat_ellipse(aes(color = group),level = 0.95) + 
  scale_shape_manual(values = seq(from = 21,by = 1,length.out = length(levels(ggdata$group)))) +
  theme_half_open() +
  labs(x ="",y = "") +
  scale_x_continuous(labels = "",breaks = 0)+
  scale_y_continuous(breaks = 0,labels = "")+
  theme(legend.position ="top",          #legend.position =c(0.45,0.8),
        axis.ticks.y = element_line(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.direction = "horizontal" ) +
  guides(fill = guide_legend(override.aes = list(shape = seq(from = 21,by = 1,length.out = length(levels(ggdata$group))),size =3)),size =5)



nmds2 = metaMDS(bray,k =2 ,try = 100)
stressplot(nmds2, main = "Shepard图")
ggdata = cbind(nmds2$points %>% as.data.frame() , sampleInfo)
nmdschoose = c("MDS1","MDS2")
ggplot(data = ggdata,aes(x=get(nmdschoose[1]),y=get(nmdschoose[2])))+
  geom_point(size= 4 ,aes(shape= group,fill= group)) +
  stat_ellipse(aes(color = group),level = 0.95) + 
  scale_shape_manual(values = seq(from = 21,by = 1,length.out = length(levels(ggdata$group)))) +
  theme_half_open() +
  labs(x ="",y = "") +
  scale_x_continuous(labels = "",breaks = 0)+
  scale_y_continuous(breaks = 0,labels = "")+
  theme(legend.position ="top",          #legend.position =c(0.45,0.8),
        axis.ticks.y = element_line(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.direction = "horizontal" ) +
  guides(fill = guide_legend(override.aes = list(shape = seq(from = 21,by = 1,length.out = length(levels(ggdata$group))),size =3)),size =5)


## 稀释曲线
BiocManager::install("microbiome")
BiocManager::install("amplicon")
library(microbiome)
library(amplicon)
p = alpha_rare_curve(alpha_rare, metadata, groupID = "Group")


# ABT analysis(Aggregated boosted tree（ABT）评估变量的相对重要性)
install.packages(c("dismo","gbm"))
library(dismo)
library(gbm)

#处理之前要保证两个数据文件均是行为样本。
#先计算OTU丰度表的Bray-Curtis距离。
library(vegan)
library(reshape2)
bray.otu <- vegdist(otu,method = "bray")
bray.otu <- as.matrix(bray.otu)
bray.otu <- melt(bray.otu)
colnames(z) <- colnames(env)
ABT.data <- cbind(z,bray.otu$value)
ABT.data <- as.data.frame(ABT.data)
ABT.fit <- gbm.step(data = ABT.data,gbm.x = 1:10,gbm.y = 11,family = "laplace",
                    tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)
#gbm.x的数值为1至环境因子的数目，gbm.y的数值为环境因子数目加1

#将变量重要性的结果信息提取后输出
softcorals_var_influence <- summary(ABT.fit)
write.csv(softcorals_var_influence, 'softcorals.var_influence.csv', row.names = FALSE, quote = FALSE)

#然后打开一个已经安装 ggplot2 的较新的 R 版本
#加载 ggplot2，读取数据后重新绘制变量重要性的柱形图
library(ggplot2)

softcorals_var_influence <- read.csv('softcorals.var_influence.csv', stringsAsFactors = FALSE)
softcorals_var_influence <- softcorals_var_influence[order(softcorals_var_influence$rel.influence), ]
softcorals_var_influence$var <- factor(softcorals_var_influence$var, levels = softcorals_var_influence$var)

#颜色代表不同的变量
p <- ggplot(softcorals_var_influence, aes(var, rel.influence)) +
  coord_flip() +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '', y = 'Relative influence(%)', title = '')

p + geom_col(aes(fill = var), width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c('#8DD3C7', '#FFFFB3', '#BEBADA',
                               '#FB8072', '#80B1D3', '#FDB462', '#B3DE69'))

#或者颜色以渐变色表示变量的重要性
p + geom_col(aes(fill = rel.influence), width = 0.7, show.legend = FALSE) +
  scale_fill_gradientn(colours = colorRampPalette(c('#86CEF7', '#0000F6'))(10))

# 普氏分析结果可视化--物种与环境、物种与物种、物种与功能关联分析
# https://mp.weixin.qq.com/s/CzlW_I8fcHrHdpJ_dqVV-A
# 在微生物群落研究的过程中，我们经常需要评估微生物群落结构与环境因子整体之间是否具有显著的相关性，此时，通常使用的方式是Mantel test和普氏分析。
# 当然除了分析群落结构与环境因子的相关性之外，这两个分析还可以用于分析同一样品不同类型微生物群落之间的相关性，比如同一样品的稀有和丰富物种或者同一样品细菌和真菌群落结构的相关性。
# 不同微生物群落之间，该分析更多的还是用于分析微生物群落组成结构与其它功能基因组之间的关系，比如细菌组成结构与抗生素抗性组的相关性。
# 最后还有一种不太常用的用法，就是分析配对的两种不同类型样品微生物群落的相关性，比如河流或海洋同一位置水体和沉积物细菌群落组成结构的相关性，从而分析这两种关联样品类型中微生物群落的转移。
library(vegan)
data <- data[which(rowSums(data) > 0),]
data <- t(data)
s.dist <- vegdist(data,method = "bray")
r.dist <- vegdist(data)
mantel(s.dist,r.dist)
mantel(s.dist,r.dist,method = "spearman")
mds.s <- monoMDS(s.dist)
mds.r <- monoMDS(r.dist)
pro.s.r <- procrustes(mds.s,mds.r)
protest(mds.s,mds.r)
library(ggplot2)
Y <- cbind(data.frame(pro.s.r$Yrot), data.frame(pro.s.r$X))
X <- data.frame(pro.s.r$rotation)
Y$ID <- rownames(Y)
p <- ggplot(Y) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#B2182B", size = 1) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#56B4E9", size = 1) +
  geom_point(aes(X1, X2), fill = "#B2182B", size = 4, shape = 21) +
  geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 4, shape = 21) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between community and environment") + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\n    M2 = 0.8358, p-value = 0.035\nMantel test:\n    r = 0.1703, p-value = 0.04',
           x = -1.5, y = 1.2, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=16,colour = "black",hjust = 0,face = "bold"))

# donor receptor coreelation
############################## micro_metab cor ######################################

colnames(FMT_species_order)
data1 <- FMT_species_order[rownames(FMT_species_order)%in%e,]
colnames(data1)
diff_sourceMF_micrometab<- read.delim("F:/A_MetaboAnalysis/FMT/metabolism/diff_sourceMF_micrometab202201007.txt")
colnames(diff_sourceMF_micrometab)
data2 <- diff_sourceMF_micrometab
colnames(data2)
rownames(data2) <- data2$Metabolite
data2$Metabolite <- NULL
data2 <- data2[,3:40]
data2 <- data2[,colnames(data1)]
rownames(data1)
rownames(data2)

library(psych)
library(reshape2)
n1 <- t(scale(t(data1)))
n2 <- t(scale(t(data2)))
cor <-corr.test(t(n1),t(n2), method = "spearman",adjust= "none")						
cor <-corr.test(t(data1),t(data2), method = "spearman",adjust= "none")						
#cor <-corr.test(diff_microb, diff_metab, method = "pearson",adjust= "none")						
## 提取相关性、p值						
cor						
cmt <-cor$r						
#View(cmt)						
c=melt(cmt)						
pmt <- cor$p						
## 输出相关系数表格,第一行为代谢物信息，第一列为物种信息						
cmt.out<-cbind(rownames(cmt),cmt)						
#write.table(pmt.out,file= "pvalue.txt",sep= "t",row.names=F)						
## 第一列为物种名，第二列为代谢物名，第三、第四列对应显示相关系数与p值						
df <-melt(cmt,value.name= "cor")						
df$pvalue <- as.vector(pmt)						
#write.table(cmt.out,file= "cor.txt",sep= "t",row.names=F)						
## 输出p值表格，第一行为代谢物信息，第一列为物种信息						
pmt.out<-cbind(rownames(pmt),pmt)						
df$fdr <- p.adjust(df$pvalue,method = "fdr")						
df <- subset(df,abs(cor)>0.6& fdr<0.05)						
View(df)						

View(subset(df,abs(cor)>0.6))						
#write.table(df,file= "cor-p.txt",sep= "t")						
if(!is.null(pmt)){						
  ssmt <- pmt< 0.01						
  pmt[ssmt] <- '**'						
  smt <- pmt > 0.01& pmt < 0.05						
  pmt[smt] <- '*'						
  pmt[!ssmt&!smt]<- ''						
} else{						
  pmt <- F						
}						
pmt	


mycol <-  colorRampPalette(c("#336699", "white", "#CC3333"))(200)

p=pheatmap::pheatmap(t(na.omit(cmt)),fontsize = 16,show_colnames = T,legend_breaks = c(-1,-0.5,0,0.5,1),display_numbers = t(na.omit(pmt)),
                     color=mycol)	

annotation_col <- data.frame(group=sourceMF[e,]$group,row.names = e)
#colnames(cmt)
#table(annotation_col$group)
#View(annotation_col)
#head(diff_microb)
#dim(diff_microb)
#View(cmt)
annotation_row <- data.frame(group=diff_sourceMF_micrometab$group,row.names = diff_sourceMF_micrometab$Metabolite)
anno_colors=list(group = c(Source_F = "#CFB5C8", Source_M = "#B0C5C3"))
p=pheatmap::pheatmap(t(cmt),fontsize = 16,show_colnames = T,legend_breaks = c(-1,-0.5,0,0.5,1),display_numbers = t(pmt),
                     color=mycol,annotation_row = annotation_row,annotation_colors = anno_colors,
                     annotation_col = annotation_col,clustering_method = "ward.D2")	

colnames(data1)
ggdata <- t(rbind(data1,data2))%>%as.data.frame(.)						
#dim(data)						
#dim(logh[,-78])						
colnames(ggdata)
rownames(ggdata)
ggdata
ggplot(data=ggdata, aes(x=`Bifidobacterium_pseudolongum`, y=`D-Maltose`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					


ggplot(data=ggdata, aes(x=`Bifidobacterium_pseudolongum`, y=`Stachyose`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					
