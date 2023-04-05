setwd("F:/A_MetaboAnalysis/Microbiome/")
library(limma)
library(dplyr)
library(reshape2)
library(tidyverse)
# 种 species
merged_abundance_table_species <- read.delim("F:/A_MetaboAnalysis/Microbiome/merged_abundance_table_species.txt",row.names = 1 ,header=T)
#View(merged_abundance_table_species)
sampleInfo <- strsplit2(colnames(merged_abundance_table_species),"_") %>%as.data.frame()
colnames(sampleInfo) <- c("condition","sex","replicate")
sampleInfo$group <- paste0(sampleInfo$condition,"_",sampleInfo$sex)
rownames(sampleInfo) <- colnames(merged_abundance_table_species)
otu <- merged_abundance_table_species
otu$Apc_F_9 <- NULL
colnames(otu)

# top 20/30
otu$sum <- rowSums(otu)
otu_order <- otu[order(otu$sum, decreasing=TRUE), ]
colnames(otu)
otu_order$sum <- NULL
dim(otu_order[1:20,])
mat1 <- otu_order[1:20,]
mat2 <- otu_order[1:30,]
dim(sampleInfo)
tmp1 = apply(mat1,1,function(x){tapply(x,sampleInfo$group[-30],median)}) %>% t() %>% as.data.frame()
tmp2 = apply(mat2,1,function(x){tapply(x,sampleInfo$group[-30],median)}) %>% t() %>% as.data.frame()
tmp = apply(otu_order,1,function(x){tapply(x,sampleInfo$group[-30],median)}) %>% t() %>% as.data.frame()
n <- t(scale(t(tmp2)))
n[n>2]=2
n[n < -2]=-2
pheatmap::pheatmap(n[,c(4,2,1,3)],fontsize = 20,clustering_method = "ward.D",
                   border_color = "white",cluster_rows = T,cluster_cols = F,angle_col = 45,
                   color = colorRampPalette(c("#336699", "white", "#CC3333"))(50))
## complexheatmap
n <- t(scale(t(mat2)))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#51659C","#9A668C","#B0C5C3","#CFB5C8")),
                       labels = c("WT-M","WT-F","Apc-M","Apc-F"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

Group = factor(gsub("_","-",sampleInfo$group[-30]),levels =c("WT-M","WT-F","Apc-M","Apc-F"))
?Heatmap
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

others <-apply(tmp[!rownames(tmp)%in%rownames(tmp1),],2,sum)%>%as.data.frame(.)
colnames(others) <- "Others"
others <- as.data.frame(t(others))
mat <- rbind(tmp1,others)
mat
ggdata = mat %>% rownames_to_column(var = "species") %>% gather(key = "group",value = "value", -"species")
unique(ggdata$species)
ggdata$species <- factor(ggdata$species,levels = unique(ggdata$species))

tmp1 = apply(mat1,1,function(x){tapply(x,sampleInfo$group[-30],mean)}) %>% t() %>% as.data.frame()						
mat <- rbind(tmp1,others)						
mat						
ggdata = mat %>% rownames_to_column(var = "species") %>% gather(key = "group",value = "value", -"species")						
ggdata$group <- factor(gsub("_","-",ggdata$group),levels = c("WT-M","WT-F","Apc-M","Apc-F"))						
ggdata$species <- factor(ggdata$species,levels = unique(ggdata$species))						
ggplot(ggdata) +				
  geom_bar(aes(x = group,y=value,fill = species),stat = "identity",position = "fill",colour= 'black') +						
  scale_fill_manual(values  = c("#6686c6","#a0c1db","gray","#E65F92","#c77364","#ce8f5c",						
                                "#7bac80","#75a5c0","#b5181a","#b72d9e",						
                                "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",						
                                "#3e6926","#0a0883","#49ddd0","#e0f8f7",						
                                "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+	
  xlab("Group")+ylab("% of species")+guides(fill = guide_legend( ncol = 1, byrow = TRUE))						

## 参考https://mp.weixin.qq.com/s?__biz=MzI5MTcwNjA4NQ==&mid=2247523961&idx=1&sn=ddea2f5ad36a987b18c19fa5b29c2248&scene=21#wechat_redirect
#beta
library(ape)
metaphlan_anno = read.table("F:/A_MetaboAnalysis/metaphlan_anno.txt",header = T) 
View(metaphlan_anno)
metaphlan_tree = ape::read.tree("F:/A_MetaboAnalysis/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")
metaphlan_tree$tip.label <- metaphlan_tree$tip.label %>% gsub(".+\\|s__","",.)
my_tree = ape::keep.tip(metaphlan_tree,intersect(rownames(otu),metaphlan_tree$tip.label))
#BiocManager::install("ggdendro")
library(ggplot2)
library(ggdendro)
library(phyloseq)
library(tidyverse)
library(ggfun)
#BiocManager::install("amplicon")
library(ggtree)
library(ggstance)
#library(amplicon)
#taxonomy <- metaphlan_anno
#
#result <- tax_stack_clust(otu=otu,tax=taxonomy,map=FMTInfo,j = "Species",
#                          Top=10,tran=T,dist="bray",
#                          hcluter_method="ward.D",
#                          cuttree=3,Group="group")
#result[[4]]
#?tax_stack_clust

#unifrac距离

#mpa_table = otu[my_tree$tip.label,]/100.0
#rbiom_distmat <- rbiom::unifrac(as.matrix(mpa_table), weighted=F, tree=my_tree)
#pheatmap::pheatmap(rbiom_distmat)
#unifrac距离2
#BiocManager::install("phyloseq")
library(phyloseq)
library(vegan)
dim(sampleInfo)
my_phy = phyloseq(otu_table(as.matrix(otu_order),taxa_are_rows = TRUE),
                  sample_data(sampleInfo[-30,]),phy_tree(my_tree),tax_table(metaphlan_anno %>% as.matrix()))  #这里是否要as.matrix很严格

weight_unifrac= phyloseq::distance(my_phy,method = "wunifrac") %>% as.matrix()
unweight_unifrac = phyloseq::distance(my_phy,method = "unifrac") %>% as.matrix()
#合并展示
res = pcoa(unweight_unifrac,correction = "cailliez",rn = NULL)
biplot.pcoa(res)
res$vectors
pcoachoose = c("Axis.1","Axis.2")
ggdata = res$vectors[,pcoachoose] %>%as.data.frame()
dim(ggdata)
length(Group)
ggdata$group <- Group
ggdata$group <- gsub("_","-",ggdata$group)%>%factor(.,levels = c("WT-M","WT-F","Apc-M","Apc-F"))
ggdata$group2 <- strsplit2(ggdata$group,"-")[,1]%>%factor(.,levels = c("WT","Apc"))

ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),shape=group,fill=group))+						
  geom_point(size=8) +stat_ellipse(type = "t", linetype = 1, show.legend = FALSE,aes(color = group))+						
  ylab("Bray-Curtis diversity") +theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+						
  scale_shape_manual(values = c(21,21,22,22))+						
  scale_color_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+						
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+						
  xlab("PCoA1 (33.0%)")+ylab("PCoA2 (23.9%)")

ggplot(data = ggdata,aes(x=get(pcoachoose[1]),y=get(pcoachoose[2]),shape=group2,fill=group2))+						
  geom_point(size=8) +stat_ellipse(type = "t", linetype = 1, show.legend = FALSE,aes(color = group2))+						
  ylab("Bray-Curtis diversity") +theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+						
  scale_shape_manual(values = c(21,22))+						
  scale_color_manual(values=c("#CFD9D0","#D8B294"))+						
  scale_fill_manual(values=c("#CFD9D0","#D8B294"))+						
  xlab("PCoA1 (33.0%)")+ylab("PCoA2 (23.9%)")


# 计算加权bray-curtis距离
dune_dist <- vegdist(t(otu_order), method="bray", binary=F)
dune_pcoa <- cmdscale(dune_dist, k=3, eig=T)
dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)
set.seed(1)
dune.div <- adonis2(t(otu_order) ~ group, data = sampleInfo[-30,], permutations = 999, method="bray")
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",force = TRUE)
library(pairwiseAdonis)

# This is a wrapper function for multilevel pairwise comparison 
# using adonis() from package 'vegan'. 
# The function returns adjusted p-values using p.adjust().
dune.pairwise.adonis <- pairwise.adonis(x=t(otu_order), factors=Group, sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

dune.pairwise.adonis <- pairwise.adonis(x=t(otu_order), factors=strsplite2(Group,"-")[,1], sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

write.table(dune.pairwise.adonis,"F:/A_MetaboAnalysis/Apc_adonis.txt", sep = '\t', col.names = NA, quote = FALSE)


##############################Aoc_M vs Apc_F###################################
colnames(otu_order)
otu1<- otu_order[,c(22:37)]
colnames(otu1)
p_value=c() 
for (i in 1:nrow(otu1)){
  data_new <- melt(otu1[i,])
  data_new$variable<-Group[c(22:37)]
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
ApcMF<- subset(otu1,pvalue<0.05)%>%as.data.frame()
ApcMF <- ApcMF[order(ApcMF$log2FC),]
head(ApcMF)
ApcMF$group <- 0
ApcMF$group[which(ApcMF$log2FC>0)]="Apc_M"
ApcMF$group[which(ApcMF$log2FC<0)]="Apc_F"
View(ApcMF)
data <- otu_order[rownames(otu_order)%in%rownames(ApcMF),]
n <- t(scale(t(data)))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#51659C","#9A668C","#B0C5C3","#CFB5C8")),
                       labels = c("WT-M","WT-F","Apc-M","Apc-F"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(gsub("_","-",sampleInfo$group[-30]),levels =c("WT-M","WT-F","Apc-M","Apc-F"))
?Heatmap
Heatmap(n[rownames(ApcMF),],name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = F,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))

##################################### WT_M vs WT_F ####################################
otu2<- otu_order[,c(1:21)]
colnames(otu2)
p_value=c() 
for (i in 1:nrow(otu2)){
  data_new <- melt(otu2[i,])
  data_new$variable<-Group[c(1:21)]
  P=wilcox.test(value~variable,data_new)
  p_value[i]=P$p.value
  p_value[i]<- round(p_value[i],3) #设置小数点位数
}
FC <- c()
for(i in 1:nrow(otu2)){
  FC[i]<-"NA"
  FC[i]<-mean(as.numeric(otu2[i,1:9]))/mean(as.numeric(otu2[i,10:21])) 
}
otu2$pvalue <- p_value
otu2$FC <- FC
otu2$log2FC <- log2(as.numeric(otu2$FC))
#View(otu2)
WTMF<- subset(otu2,pvalue<0.05)
WTMF <- WTMF[order(WTMF$log2FC),]
#head(WTMF)
WTMF$group <- 0
WTMF$group[which(WTMF$log2FC>0)]="WT_M"
WTMF$group[which(WTMF$log2FC<0)]="WT_F"
write.table(WTMF,'WTMF_diff.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(ApcMF,'ApcMF_diff.txt', sep = '\t', col.names = NA, quote = FALSE)
dim(WTMF)
dim(ApcMF)
data <- otu_order[rownames(otu_order)%in%rownames(WTMF),]
n <- t(scale(t(data)))
n[n>2]=2
n[n < -2]=-2
library(ComplexHeatmap)
library(circlize)
?Heatmap
Heatmap(n[rownames(WTMF),],name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = F,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))



setdiff <- setdiff(rownames(ApcMF),rownames(WTMF))
setdiff
data <- ApcMF[rownames(ApcMF)%in%setdiff,]
dim(data)
n <- t(scale(t(otu_order[rownames(data),])))
Heatmap(n,name = " ",						
        col = col_fun,						
        top_annotation = top_annotation,						
        column_split = Group,cluster_column_slices = F,						
        show_heatmap_legend = T,						
        border = "black",cluster_columns = T,						
        show_column_names = F,column_names_gp = gpar(fontsize=14),						
        show_row_names = T,row_names_gp = gpar(fontsize=14),						
        column_title = NULL,cluster_rows = F,						
        row_names_max_width = max_text_width(rownames(n),gp = gpar(fontsize = 14)))

###################### Akkermansia muciniphila #########################
data <- otu_order[rownames(otu_order)%in%setdiff,]
data$species <- rownames(data)
ggdata <- melt(data)
ggdata$group <- paste0(strsplit2(ggdata$variable,"_")[,1],"-",strsplit2(ggdata$variable,"_")[,2])%>%
  factor(.,levels = c("WT-M","WT-F","Apc-M","Apc-F"))
setdiff

ggplot(ggdata[ggdata$species=="Akkermansia_muciniphila",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative abundance")+ggtitle("A.muciniphila")+ 
  theme(plot.title = element_text(size = 24, face = "italic.bold"))

########################### Parabacteroides_goldsteinii ###############################
ggplot(ggdata[ggdata$species=="Parabacteroides_goldsteinii",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative abundance")+ggtitle("P.goldsteinii")+ 
  theme(plot.title = element_text(size = 24, face = "italic"))

###################### Lactobacillus_fermentum ###########################
ggplot(ggdata[ggdata$species=="Lactobacillus_fermentum",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative abundance")+ggtitle("L.fermentum")+ 
  theme(plot.title = element_text(size = 24, face = "italic"))+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+guides(fill=FALSE)

###################### Alistipes_inops ###########################
ggplot(ggdata[ggdata$species=="Alistipes_inops",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative abundance")+ggtitle("A.inops")+ 
  theme(plot.title = element_text(size = 24, face = "italic"))+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+guides(fill=FALSE)

###################### Lactobacillus_fermentum ###########################
ggplot(ggdata[ggdata$species=="Lactobacillus_taiwanensis",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')+
  ylab("Relative abundance")+ggtitle("L.taiwanensis")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+guides(fill=FALSE)

#merged_model_metaphlan_species <- readRDS("F:/A_MetaboAnalysis/merged_model_metaphlan_species.rds")
#with_human <- intersect(setdiff,rownames(merged_model_metaphlan_species))
#ggdata = otu_order %>% rownames_to_column(var = "species") %>% gather(key = "group",value = "value", -"species")
#ggdata$Group <- paste0(strsplit2(ggdata$group,"_")[,1],"_",strsplit2(ggdata$group,"_")[,2])
#ggplot(ggdata[ggdata$species%in%intersect(setdiff,rownames(merged_model_metaphlan_species)),],aes(x = Group,y = value,fill=Group))+
#  geom_boxplot(outlier.color = NA)+
#  theme_classic()+theme(axis.title =element_text(size = 20),
#                        axis.text =element_text(size = 20, color = 'black'))+
#  theme(axis.text.x = element_text(size = 20,hjust = 1,angle = 45))+
#  theme(axis.text.y = element_text(size = 20))+
#  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.text = element_text(size=12),
#        strip.text = element_text(size=14))+
#  scale_fill_manual(values=c("#CFB5C8","#B0C5C3","#9A668C","#51659C"))+
#  geom_signif(comparisons = list(c("WT_M", "WT_F"),c("Apc_M", "Apc_F"),c("WT_M", "Apc_M"),c("WT_F", "Apc_F")),
#              map_signif_level = F,step_increase = 0.1)+ylab("")+
#  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
#  facet_wrap(~species,scales = "free",nrow =3)+scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))

############# alpha diversity #################
library(vegan)
#BiocManager::install("picante")
library(picante)      
otu_count <- read.delim("F:/A_MetaboAnalysis/Microbiome/merged_count_species.txt",row.names = 1 ,header=T)
alpha_diversity <- function(x, tree = NULL) {
  otus <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
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
library(reshape2)
otu_t <- t(otu_count)
alpha <- alpha_diversity(otu_t)
alpha$ACE<-NULL
alpha$goods_Coverage <- NULL
#alpha$observed_species <- NULL
#alpha$Simpson <- NULL
alpha$group <- sampleInfo$group

## 
otu_diversity_long <- alpha %>% gather(key = "index",value = "value", -"group")%>%
  as.data.frame(.)
class(otu_diversity_long$value)
otu_diversity_long$value <- as.numeric(otu_diversity_long$value)
Data_summary <- summarySE(otu_diversity_long, measurevar="value", groupvars=c("group","index"))
otu_diversity_long$group <- factor(gsub("_","-",otu_diversity_long$group),levels = c("WT-M","WT-F","Apc-M","Apc-F"))

otu_diversity_long$group2 <- factor(strsplit2(otu_diversity_long$group,"-")[,1],levels = c("WT","Apc"))
a=otu_diversity_long[,c(1,2,3)]
b=otu_diversity_long[,c(4,2,3)]
colnames(b)[1] <- "group"
c=rbind(a,b)
c$group <- factor(c$group,levels = c("WT","Apc","WT-M","WT-F","Apc-M","Apc-F"))

#### otus
ggplot(otu_diversity_long[otu_diversity_long$index=="otus",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(otu_diversity_long[otu_diversity_long$index=="otus",],aes(x = group2,y = value,fill=group2))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(width=0.5,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#CFD9D0","#D8B294"))+
  geom_signif(comparisons = list(c("WT","Apc")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT","Apc"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(c[c$index=="otus",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#CFD9D0","#D8B294","#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT","Apc"),c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT","Apc","WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')





#### shannon
ggplot(otu_diversity_long[otu_diversity_long$index=="Shannon",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(otu_diversity_long[otu_diversity_long$index=="Shannon",],aes(x = group2,y = value,fill=group2))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(width=0.5,lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#CFD9D0","#D8B294"))+
  geom_signif(comparisons = list(c("WT","Apc")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT","Apc"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

ggplot(c[c$index=="Shannon",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#CFD9D0","#D8B294","#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT","Apc"),c("WT-M", "WT-F"),c("Apc-M", "Apc-F"),
                                 c("WT-M", "Apc-M"),c("WT-F", "Apc-F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT","Apc","WT-M","WT-F","Apc-M","Apc-F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')

################### SparCC ########################
#for (i in c("tidyverse",'network',"tidyfst", 'psych',"ggnewscale","ggraph", 'sna', 'ggrepel', 'phyloseq', 'reshape2', 'packcircles',"WGCNA", 'ggpubr')) {
#  BiocManager::install(i)
#}
remotes::install_github("taowenmicro/ggClusterNet")
############ AB treatment##################
level.7 <- read.csv("F:/A_MetaboAnalysis/level-7.csv", dec=",",row.names = 1)
level.7 <- read.delim("C:/Users/tutu/Desktop/level7.txt",row.names = 1)
View(level.7)
otu_t <- t(level.7)
alpha <- alpha_diversity(otu_t)
alpha$ACE<-NULL
alpha$goods_Coverage <- NULL
alpha$group <- substr(gsub("FMT","AB",colnames(level.7)),1,4)
alpha$group <- c(rep("female",5),rep("male",5))

otu_diversity_long <- alpha %>% gather(key = "index",value = "value", -"group")%>%
  as.data.frame(.)
class(otu_diversity_long$value)
otu_diversity_long$value <- as.numeric(otu_diversity_long$value)
Data_summary <- summarySE(otu_diversity_long, measurevar="value", groupvars=c("group","index"))
otu_diversity_long$group <- factor(otu_diversity_long$group,levels = c("WT_M","WT_F","AB_M","AB_F"))
otu_diversity_long$group <- factor(otu_diversity_long$group,levels = c("male","female"))

#### otus
ggplot(otu_diversity_long[otu_diversity_long$index=="otus",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT_M", "WT_F"),c("AB_M", "AB_F"),
                                 c("WT_M", "AB_M"),c("WT_F", "AB_F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT_M","WT_F","AB_M","AB_F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')
ggplot(otu_diversity_long[otu_diversity_long$index=="Shannon",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("male","female")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("No. observed OTUs")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')


#### shannon
ggplot(otu_diversity_long[otu_diversity_long$index=="Shannon",],aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT_M", "WT_F"),c("AB_M", "AB_F"),c("WT_M", "AB_M"),c("WT_F", "AB_F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("WT_M","WT_F","AB_M","AB_F"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')


View(level.7)
dim(level.7)
tt <- as.data.frame(t(level.7))
level7=tt[,1:51]/tt$sum 
level7 <- as.data.frame(t(level7))
head(level7)
View(level7)
level7$FMT_M_3 <- NULL
## Reads
reads=data.frame(Count=tt$sum)
rownames(reads) <- rownames(tt)
reads$Sample=gsub("FMT","AB",rownames(reads))
reads$group=factor(substr(reads$Sample,1,4),levels = c("WT_M","WT_F","AB_M","AB_F"))

ggplot(reads,aes(group,Count))+
  stat_summary(mapping=aes(fill = group),fun=mean,geom = "bar",
             fun.args = list(mult=1),width=0.5)+
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2,lwd=1)+
  theme_classic()+scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT_M", "WT_F"),c("AB_M", "AB_F"),
                                 c("WT_M", "AB_M"),c("WT_F", "AB_F")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2)+
  scale_x_discrete(limits=c("WT_M","WT_F","AB_M","AB_F"))+ylab("Sequence Read Counts")+xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+theme(legend.position = 'none')


tmp = apply(level7,1,function(x){tapply(x, c(rep("AB_F",5),rep("AB_M",5),rep("WT_F",6),rep("WT_M",6)), mean)}) %>% t() %>% as.data.frame()
library("tidyr")
library(tibble)
library(dplyr)
library(reshape)
View(tmp)
ggdata = level7 %>% rownames_to_column(var = "species") %>% gather(key = "group",value = "value", -"species")
View(ggdata)
#FFFFFF
#DDDDDD	#AAAAAA	#888888	#666666	#444444	#000000
#FFB7DD	#FF88C2	#FF44AA	#FF0088	#C10066	#A20055	#8C0044
#FFCCCC	#FF8888	#FF3333	#FF0000	#CC0000	#AA0000	#880000
#FFC8B4	#FFA488	#FF7744	#FF5511	#E63F00	#C63300	#A42D00
#FFDDAA	#FFBB66	#FFAA33	#FF8800	#EE7700	#CC6600	#BB5500
#FFEE99	#FFDD55	#FFCC22	#FFBB00	#DDAA00	#AA7700	#886600
#FFFFBB	#FFFF77	#FFFF33	#FFFF00	#EEEE00	#BBBB00	#888800
#EEFFBB	#DDFF77	#CCFF33	#BBFF00	#99DD00	#88AA00	#668800
#CCFF99	#BBFF66	#99FF33	#77FF00	#66DD00	#55AA00	#227700
#99FF99	#66FF66	#33FF33	#00FF00	#00DD00	#00AA00	#008800
#BBFFEE	#77FFCC	#33FFAA	#00FF99	#00DD77	#00AA55	#008844
#AAFFEE	#77FFEE	#33FFDD	#00FFCC	#00DDAA	#00AA88	#008866
#99FFFF	#66FFFF	#33FFFF	#00FFFF	#00DDDD	#00AAAA	#008888
#CCEEFF	#77DDFF	#33CCFF	#00BBFF	#009FCC	#0088A8	#007799
#CCDDFF	#99BBFF	#5599FF	#0066FF	#0044BB	#003C9D	#003377
#CCCCFF	#9999FF	#5555FF	#0000FF	#0000CC	#0000AA	#000088
#CCBBFF	#9F88FF	#7744FF	#5500FF	#4400CC	#2200AA	#220088
#D1BBFF	#B088FF	#9955FF	#7700FF	#5500DD	#4400B3	#3A0088
#E8CCFF	#D28EFF	#B94FFF	#9900FF	#7700BB	#66009D	#550088
#F0BBFF	#E38EFF	#E93EFF	#CC00FF	#A500CC	#7A0099	#660077
#FFB3FF	#FF77FF	#FF3EFF	#FF0 0FF	#CC00CC	#990099	#770077
ggplot2::ggplot(ggdata) +
  geom_bar(aes(x = group,y=value,fill = species),stat = "identity",position = "fill") +
  scale_fill_manual(values  = c("#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",
                                "#CCEEFF","#0a0883","#49ddd0","#e0f8f7","#FF9800","#61A0FF","#0398A8",
                                "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#FFB7DD",
                                "#8f8f8f","#75a5c0","#b5181a","#b72d9e","#99FF33",
                                "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",
                                "#5599FF","#0a0883","#49ddd0","#e0f8f7",
                                "#651d91","#9d3eaf","#b9a3d0","#5b5a5b",
                                "#8f8f8f","#75a5c0","#b5181a","#b72d9e",
                                "#e4cccb","#f6f0eb","#e8c755","#d5d456",
                                "#cfe74f","#39cf40","#FFDDAA","#0a0883","#49ddd0","#e0f8f7",
                                "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f",
                                "#75a5c0","#b72d9e"))+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45, vjust =0.5),
        legend.text = element_text(size=12),legend.position = "right")+xlab("Group")+ylab("% of species")+ theme(legend.position="none")


marker_pub <- read.delim("F:/A_MetaboAnalysis/Microbiome/marker_taxa_in_phenotype_comparison_between_D006262_D015179.txt")
head(marker_pub)
select <- marker_pub[marker_pub$scientific_name%in%tmp$name,]
View(select$scientific_name)
dim(select)
tt <- merged_abundance_table_species[unique(select$scientific_name),]
tt <- test_species[unique(select$scientific_name),]
tt
n <- t(scale(t(tt)))
n[n>2]=2
n[n<-2]=-2
pheatmap::pheatmap(n,fontsize = 16,
                   border_color = "white",cluster_rows = F,cluster_cols = F,
                   color = colorRampPalette(c("#336699", "white", "#CC3333"))(50))

View(tmp[unique(select$scientific_name),1:4])
select$sig <- c(rep("down",10),rep("up",20))
dd <- select[!duplicated(select$scientific_name),]
rownames(dd) <- dd$scientific_name
View(dd)
ggplot(data=dd[unique(select$scientific_name),], aes(x=scientific_name, y=LDA, fill=sig))+
  geom_bar(stat = "identity",position = "identity",width=0.9)+theme_classic()+
  coord_flip()+scale_fill_manual(values = c("#9D22F0", "#1488F0"), guide=FALSE)+
  scale_x_discrete(limits=rev(unique(select$scientific_name)))
 geom_abline(linetype="dashed",intercept = 50, slope = 0,size=1,colour='gray')+
  geom_abline(linetype="dashed",intercept = -50, slope = 0,size=1,colour='gray')+
  theme(panel.grid=element_blank(),panel.border=element_blank(),
        xis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

#热图
 n <- t(scale(t(tmp2)))
 n[n>1]=1
 n[n<-1]=-1
 n

pheatmap::pheatmap(na.omit(n[,c(4,3,2,1)]),fontsize = 22,clustering_method = "ward.D2",
                    border_color = "white",cluster_rows = T,cluster_cols = F,angle_col = 45,
                    color = colorRampPalette(c("#336699", "white", "#CC3333"))(50))
 
pheatmap::pheatmap(tmp[1:30,],scale = "row",fontsize = 16,border_color = "white")

diff_apc <- merged_abundance_table_species[rownames(merged_abundance_table_species)%in%setdiff_apc_wt$V1,]



pheatmap::pheatmap(cor(mat1[,1:38]),fontsize = 20,
                   border_color = NA,cluster_rows = T,cluster_cols = T,
                   color = colorRampPalette(c("#336699", "white", "#CC3333"))(50))


## 稀释曲线

#library(microbiome)
#library(amplicon)
#merged_abundance_counts_species <- read.delim("F:/A_MetaboAnalysis/merged_abundance_counts_species.txt")
#rownames(merged_abundance_counts_species) <- merged_abundance_counts_species$clade_name
#merged_abundance_counts_species$clade_name <- NULL
#change_name <- read.delim("F:/A_MetaboAnalysis/Microbiome/change_name.txt", header=FALSE)
#colnames(merged_abundance_counts_species) <- strsplit2(colnames(merged_abundance_counts_species),"_")[,3]%>%as.character()
#View(change_name)
#change_name$V1 <- change_name$V1%>%as.character()
#merged_abundance_counts_species <- merged_abundance_counts_species[,change_name$V1]
#colnames(merged_abundance_counts_species) <- change_name$V2
#sampleInfo <- strsplit2(colnames(merged_abundance_counts_species),"_") %>%as.data.frame()
#colnames(sampleInfo) <- c("condidion","sex","replicate")
#sampleInfo$group <- paste0(sampleInfo$condidion,"_",sampleInfo$sex)
#rownames(sampleInfo) <- colnames(merged_abundance_counts_species)
#write.table(merged_abundance_counts_species,"F:/A_MetaboAnalysis/merged_counts.txt", sep = '\t', col.names = NA, quote = FALSE)



library(vegan)
BiocManager::install("picante")
library(picante)      
otu_count <- merged_abundance_counts_species
alpha_diversity <- function(x, tree = NULL) {
  otus <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
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
alpha$observed_species <- NULL
alpha$Simpson <- NULL
alpha$group <- sampleInfo$group
alpha
## 多重比较
BiocManager::install("EasyAovWlxPlot")
library(EasyAovWlxPlot)
library(ggplot2)
alpha <- data.frame(ID=rownames(alpha),alpha)%>%.[,c(1,6,2,3,4,5)]
alpha$Shannon <- as.numeric(alpha$Shannon)
alpha$Simpson <- as.numeric(alpha$Simpson)

# 正态检验
qqnorm(alpha$observed_species);qqline(alpha$observed_species)
qqnorm(alpha$Chao1);qqline(alpha$Chao1)
qqnorm(alpha$Shannon);qqline(alpha$Shannon)
qqnorm(alpha$Simpson);qqline(alpha$Simpson)
shapiro.test(alpha$observed_species)
shapiro.test(alpha$Chao1)
class(alpha$Shannon)
shapiro.test(as.numeric(alpha$Shannon))
shapiro.test(as.numeric(alpha$Simpson))

# 方差齐性检验
x <- na.omit(otu_diversity_long)
bartlett.test(value~index,data=x)
## 结论：既不正态又方差不齐，不能用annova

# 非参数检验 单个指标
res = KwWlx(data = alpha, i=6)
# 多个指标
alpha$group <- factor(alpha$group)
?MuiaovMcomper
result  = MuiKwWlx(data = alpha,num = c(3:6))
result
?FacetMuiPlotresultBox
?FacetMuiPlotresultBox
result1= FacetMuiPlotresultBox(data = alpha,num = c(3:6),result = result,sig_show ="abc",ncol = 2 )
result1[[1]]+
  theme_classic()+theme(axis.title =element_text(size = 18),
                        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text.x = element_text(size = 16,hjust = 1,angle = 45))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.text = element_text(size=14),
        legend.title  = element_text(size=14),strip.text = element_text(size=14))+
  scale_fill_manual(values=c("#CFB5C8","#B0C5C3","#9A668C","#51659C"))+scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))




#实现了多个指标批量整体运行；这个函数可以将我们的目标列做正态检验和方差齐性，
#然后根据结果选择方差检验或者多重比较方法，最后选择自己需要的出图方式和显著性
result1 = MuiStat(data = alpha,num = c(3:6),method_cv = "leveneTest",method_Mc = "Tukey",
                  sig_show  = "abc",ncol = 4,plot = "box",plottype = "mui")





?melt
head(alpha)
sampleInfo$group
otu_diversity_long <- alpha %>% gather(key = "index",value = "value", -"group")%>%as.data.frame(.)
class(otu_diversity_long$value)
otu_diversity_long$value <- as.numeric(otu_diversity_long$value)

na.omit(otu_diversity_long)
## violin
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)


Data_summary <- summarySE(otu_diversity_long, measurevar="value", groupvars=c("group","index"))
otu_diversity_long$group <- factor(otu_diversity_long$group,levels = c("WT-M","WT-F","Apc-M","Apc-F"))
ggplot(otu_diversity_long, aes(x=group, y=value,fill=group)) + 
  geom_violin(trim=FALSE,color="white") + 
  geom_boxplot(width=0.2,position=position_dodge(0.9))+
  theme_classic()+theme(axis.title =element_text(size = 18),
                        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text.x = element_text(size = 22,hjust = 1,angle = 45))+
  theme(axis.text.y = element_text(size = 22))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.text = element_text(size=14),
        legend.title  = element_text(size=14),strip.text = element_text(size=14))+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  geom_signif(comparisons = list(c("WT_M", "WT_F"),c("Apc_M", "Apc_F"),c("WT_M", "Apc_M"),c("WT_F","Apc_F")),
              map_signif_level = T,step_increase = 0.1)+ylab("")+
  stat_summary(fun=mean, geom="point", size=2) +
  facet_wrap(~index,scales = "free",nrow =1)+xlab("")

ggplot(otu_diversity_long[1:152,],aes(x = group,y = value,fill=group))+
  geom_boxplot()+
  theme_classic()+theme(axis.title =element_text(size = 18),
                        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text.x = element_text(size = 16,hjust = 1,angle = 45))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.text = element_text(size=14),
        legend.title  = element_text(size=14),strip.text = element_text(size=14))+
  scale_fill_manual(values=c("#CFB5C8","#B0C5C3","#9A668C","#51659C"))+
  geom_signif(comparisons = list(c("WT_M", "WT_F"),c("Apc_M", "Apc_F"),c("Apc_F", "WT_M")),
              map_signif_level = T )+ylab("")+
  stat_summary(fun=mean, geom="point", size=2) +
  facet_wrap(~index,scales = "free",nrow =2)+scale_x_discrete(limits=c("WT-M","WT-F","Apc-M","Apc-F"))

ps = phyloseq(otu_table(merged_abundance_counts_species, taxa_are_rows=TRUE), sample_data(sampleInfo))
View(merged_abundance_counts_species)
?alpha_rare_all
result = alpha_rare_all(otu = merged_abundance_counts_species,ps=ps, map = sampleInfo, group = "group", 
                        method = "chao1", start = 1000, step = 1000)
p = result[[4]]
# 结果返回列表，1为样本稀释曲线，2为数据表，3为按组均值的稀释曲线，4为组置信区间
# 样本稀释曲线
result[[1]]
#按组均值绘图
p = result[[3]]
# 按照分组绘制标准差稀释曲线
(p = result[[4]])
merged_abundance_counts_species
out= rarecurve(t(merged_abundance_counts_species), step=10000, sample=min(rowSums(t(merged_abundance_counts_species))), xlab="Sample size", ylab="Species",label=F)
out[[6]]

## 3 确定画图的横纵坐标的最大值
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)


## 4 plot绘图
out[[1]]
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size", ylab = "Species", type = "n")
N <- attr(out[[1]], "Subsample")
lines(N, out[[1]],col=color,lwd=1) 绘制线
legend(x,y, cex ,xpd=T,inset=0.2,lty=1,col=color,horiz=FALSE,legend=c(colnames(merged_abundance_counts_species)),bty="n")


# ABT analysis(Aggregated boosted tree（ABT）评估变量的相对重要性)
install.packages(c("dismo","gbm"))
library(dismo)
library(gbm)

#处理之前要保证两个数据文件均是行为样本。
#先计算OTU丰度表的Bray-Curtis距离。
library(vegan)
library(reshape2)
head(otu)
diff_microb
bray.otu <- vegdist(diff_microb,method = "bray")
bray.otu <- as.matrix(bray.otu)
bray.otu <- melt(bray.otu)
z <- c()
for (k in 1:ncol(diff_metab)) {
  x <- c()    
  for (i in 1:nrow(diff_metab)) {
    for (j in 1:nrow(diff_metab)) {
    y <- diff_metab[,k][i]-diff_metab[,k][j]
    x <- c(x,y)
  }
  }
  z <- cbind(z,x)
}
z
dim(diff_metab)
colnames(z) <- colnames(diff_metab)
ABT.data <- cbind(z,bray.otu$value)
ABT.data <- as.data.frame(ABT.data)
ABT.fit <- gbm.step(data = ABT.data,gbm.x = 1:57,gbm.y = 58,family = "laplace",
                    tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)
#gbm.x的数值为1至环境因子的数目，gbm.y的数值为环境因子数目加1

#将变量重要性的结果信息提取后输出
softcorals_var_influence <- summary(ABT.fit)
#write.csv(softcorals_var_influence, 'softcorals.var_influence.csv', row.names = FALSE, quote = FALSE)

#然后打开一个已经安装 ggplot2 的较新的 R 版本
#加载 ggplot2，读取数据后重新绘制变量重要性的柱形图
library(ggplot2)
softcorals_var_influence$
#softcorals_var_influence <- read.csv('softcorals.var_influence.csv', stringsAsFactors = FALSE)
softcorals_var_influence <- softcorals_var_influence[order(softcorals_var_influence$rel.inf,decreasing = T), ]
softcorals_var_influence$var <- factor(softcorals_var_influence$var, levels = softcorals_var_influence$var)

#颜色代表不同的变量
p <- ggplot(softcorals_var_influence, aes(var, rel.inf)) +
  coord_flip() +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '', y = 'Relative influence(%)', title = '')

p
p + geom_col(aes(fill = "#1E4161"), width = 0.7, show.legend = F) + theme(axis.text =element_text(size=14) )

p
#或者颜色以渐变色表示变量的重要性
p + geom_col(aes(fill = rel.inf), width = 0.7, show.legend = FALSE)+ theme(axis.text =element_text(size=14) )

# 普氏分析结果可视化--物种与环境、物种与物种、物种与功能关联分析
# https://mp.weixin.qq.com/s/CzlW_I8fcHrHdpJ_dqVV-A
# 在微生物群落研究的过程中，我们经常需要评估微生物群落结构与环境因子整体之间是否具有显著的相关性，此时，通常使用的方式是Mantel test和普氏分析。
# 当然除了分析群落结构与环境因子的相关性之外，这两个分析还可以用于分析同一样品不同类型微生物群落之间的相关性，比如同一样品的稀有和丰富物种或者同一样品细菌和真菌群落结构的相关性。
# 不同微生物群落之间，该分析更多的还是用于分析微生物群落组成结构与其它功能基因组之间的关系，比如细菌组成结构与抗生素抗性组的相关性。
# 最后还有一种不太常用的用法，就是分析配对的两种不同类型样品微生物群落的相关性，比如河流或海洋同一位置水体和沉积物细菌群落组成结构的相关性，从而分析这两种关联样品类型中微生物群落的转移。
library(vegan)

s.dist <- vegdist(diff_microb,method = "bray")
r.dist <- vegdist(diff_metab)
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
  labs(title="Correlation between microbiota and metabolome") + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\n    M2 = 0.8304, p-value = 0.004\nMantel test:\n    r = 0.1539, p-value = 0.017',
           x = -1.5, y = 1.2, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=16,colour = "black",hjust = 0,face = "bold"))
p

# Pathway humann3
humann_pathabundance <- read.delim("F:/A_MetaboAnalysis/humann/humann_pathabundance.tsv", header=T)
View(humann_pathabundance)
ss=humann_pathabundance$Pathway[grep(":",humann_pathabundance$Pathway)]%>%.[grep("[|]",.,invert = T)]
ss <- humann_pathabundance[humann_pathabundance$Pathway%in%ss,]
dim(ss)
ss
rownames(ss) <- strsplit2(ss$Pathway,": ")[,2]
ss$Pathway <- NULL
head(ss)
sampleInfo$group
colnames(ss)
mean_ss <-  as.data.frame(t(apply(as.data.frame(ss),1,function(a){tapply(a,sampleInfo$group,mean)})))
median_ss <-  as.data.frame(t(apply(as.data.frame(ss),1,function(a){tapply(a,sampleInfo$group,median)})))
n=t(scale(t(median_ss)))
View(mean_ss)

n=t(scale(t(mean_ss)))
n[n>2]=2
n[n<-2]=-2
rownames(n)
p=pheatmap::pheatmap(na.omit(n),fontsize = 20,show_rownames = T,kmeans_k = 4,border_color = "white")

cluster=p$kmeans$cluster%>%as.data.frame()
cluster$pathway <- rownames(cluster)
cluster1 <- cluster[cluster$.==1,]
dim(cluster3)
n1 <- n[rownames(n)%in%cluster1$pathway,]
p=pheatmap::pheatmap(na.omit(n1),fontsize = 20,show_rownames = T,border_color = "white")

# 脂肪酸氧化
median_path=as.data.frame(t(apply(as.data.frame(humann_pathabundance[,-1]),1,function(a){tapply(a,sampleInfo$group,median)})))
rownames(median_path) <- humann_pathabundance$Pathway
n=t(scale(t(median_path)))
n[n>2]=2
n[n<-2]=-2
p=pheatmap::pheatmap(na.omit(n),fontsize = 20,kmeans_k = 4,show_rownames = T,border_color = "white")
cc=as.data.frame(p$kmeans$cluster)
cc$pathway <- rownames(cc)
cluster3=cc[cc$`p$kmeans$cluster`==3,]
cluster3
View(cluster3)
write.table(cluster3,"F:/A_MetaboAnalysis/ApcM_pathway.txt", sep = '\t', col.names = NA, quote = FALSE)

cluster2=cc[cc$`p$kmeans$cluster`==2,]
write.table(cluster2,"F:/A_MetaboAnalysis/ApcF_pathway.txt", sep = '\t', col.names = NA, quote = FALSE)
ApcM_pathway <- read.delim("F:/A_MetaboAnalysis/ApcM_pathway.txt", row.names=1)
ApcF_pathway <- read.delim("F:/A_MetaboAnalysis/ApcF_pathway.txt", row.names=1)
ApcM_pathway$pathway <- strsplit2(ApcM_pathway$pathway,": ")[,2]
ApcF_pathway$pathway <- strsplit2(ApcF_pathway$pathway,": ")[,2]

akk <- ApcM_pathway[ApcM_pathway$species=="Akkermansia_muciniphila",]%>%rownames()
a <- na.omit(n[rownames(n)%in%akk,])%>%as.data.frame()
rownames(a)=strsplit2(rownames(a),": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]
pheatmap::pheatmap(na.omit(a),fontsize = 20,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)
a$species <- "Akkermansia_muciniphila"
a1=a
alis <- ApcM_pathway[ApcM_pathway$species=="Alistipes_inops",]%>%rownames()
a <- na.omit(n[rownames(n)%in%alis,])%>%as.data.frame()
rownames(a)=strsplit2(rownames(a),": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]
pheatmap::pheatmap(na.omit(a),fontsize = 20,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)
a$species <- "Alistipes_inops"
a2=a
firmm <- ApcM_pathway[ApcM_pathway$species=="Firmicutes_bacterium_ASF500",]%>%rownames()
a <- na.omit(n[rownames(n)%in%firm,])%>%as.data.frame()
rownames(a)=strsplit2(rownames(a),": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]
pheatmap::pheatmap(na.omit(a),fontsize = 20,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)
a$species <- "Firmicutes_bacterium_ASF500"
a3=a
rom <- ApcM_pathway[ApcM_pathway$species=="Romboutsia_ilealis",]%>%rownames()
a <- na.omit(n[rownames(n)%in%rom,])%>%as.data.frame()
rownames(a)=strsplit2(rownames(a),": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]
pheatmap::pheatmap(na.omit(a),fontsize = 20,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)
a$species <- "Romboutsia_ilealis"
a4=a

lac <- ApcM_pathway[ApcM_pathway$species=="Lactobacillus_johnsonii",]%>%rownames()
a <- na.omit(n[rownames(n)%in%lac,])%>%as.data.frame()
rownames(a)=strsplit2(rownames(a),": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]
pheatmap::pheatmap(na.omit(a)[-1,],fontsize = 20,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)

a$species <- "Lactobacillus_johnsonii"
a5=a

pm <- rbind(a1,a2)%>%rbind(.,a3)%>%rbind(.,a4)%>%rbind(.,a5)
pm$group <- "ApcM_Pathway"
## female 
unique(ApcF_pathway$species)
lact <- ApcF_pathway[ApcF_pathway$species=="Lactococcus_lactis",]%>%rownames()
a <- n[rownames(n)%in%lact,]%>%as.data.frame()%>%t()%>%as.data.frame()
rownames(a)=strsplit2(lact,": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]
pheatmap::pheatmap(na.omit(a),fontsize = 20,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)
a$species <- "Lactococcus_lactis"
a6=a

parad <- ApcF_pathway[ApcF_pathway$species=="Parabacteroides_distasonis",]%>%rownames()
a <- na.omit(n[rownames(n)%in%parad,])%>%as.data.frame()%>%.[c(-1,-2),]
dim(a)
rownames(a)=strsplit2(rownames(a),": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]

pheatmap::pheatmap(na.omit(a)[-1,],fontsize = 16,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)

a$species <- "Parabacteroides_distasonis"
pf=a
pf$group <- "ApcF_Pathway"

pathway <- rbind(pm,pf)


firm <- ApcF_pathway[ApcF_pathway$species=="Firmicutes_bacterium_ASF500",]%>%rownames()
a <- na.omit(n[rownames(n)%in%firm,])%>%as.data.frame()
rownames(a)=strsplit2(rownames(a),": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]
pheatmap::pheatmap(na.omit(a),fontsize = 20,show_rownames = T,
                   border_color = "white",cluster_cols = F,cluster_rows = F)


m <-c(akk,alis,firmm,rom,lac)
m <- strsplit2(m,": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]%>%.[c(-1,-50,-34)]
f <- c(lact,parad)
f <- strsplit2(f,": ")[,2]%>%strsplit2(.,"[|]")%>%as.data.frame(.)%>%.[,1]%>%.[c(-2,-3)]

am=setdiff(m,f)
bf=setdiff(f,m)

ll <- list(m,f)
names(ll) <- c("Apc_M","Apc_F")

p=ggvenn(ll,       
         show_percentage = F,show_elements=F,text_size = 8,
         stroke_color = "white",
         fill_color = c("#619CFF","#F8766D"),
         set_name_color = c("#619CFF","#F8766D"))  #最后一个颜色要改

p

select_path <- pathway[rownames(pathway)%in%c(am,bf),]
ann_colors = list(group = c(ApcM_Pathway="#619CFF",ApcF_Pathway= "#F8766D"), 
                  species = c(Akkermansia_muciniphila = "#68C3F0",Alistipes_inops = "#DCB0A7",
                              Firmicutes_bacterium_ASF500 = "#EAAD5A",
                              Romboutsia_ilealis = "#A8B192",
                              Parabacteroides_distasonis = "#DE8052") )
pheatmap::pheatmap(select_path[,1:4], annotation_row =select_path[,5:6],
                   annotation_colors =ann_colors,fontsize = 16,
                   cluster_cols = F,cluster_rows = F)



View(ApcM_pathway)
fatty <- median_path[grep("fatty acid",rownames(median_path)),]
rownames(fatty) <- strsplit2(rownames(fatty),": ")[,2]
View(fatty)
n=t(scale(t(fatty)))
n[n>2]=2
n[n<-2]=-2
p=pheatmap::pheatmap(na.omit(n),fontsize = 20,show_rownames = T,border_color = "white")

# 抗坏血酸降解
asco <- median_path[grep("L-asco",rownames(median_path)),]
rownames(asco) <- strsplit2(rownames(asco),": ")[,2]
View(asco)
n=t(scale(t(asco)))
n[n>2]=2
n[n<-2]=-2
p=pheatmap::pheatmap(na.omit(n),fontsize = 20,show_rownames = T,border_color = "white")

# 血红素合成
heme <- median_path[grep("heme",rownames(median_path)),]
rownames(heme) <- strsplit2(rownames(heme),": ")[,2]
View(heme)
n=t(scale(t(heme)))
n[n>2]=2
n[n<-2]=-2
p=pheatmap::pheatmap(na.omit(n),fontsize = 20,show_rownames = T,border_color = "white")

# 谷氨酸降解L-glu
L_glu <- median_path[grep("L-glu",rownames(median_path)),]
rownames(L_glu) <- strsplit2(rownames(L_glu),": ")[,2]
View(L_glu)
n=t(scale(t(L_glu)))
n[n>2]=2
n[n<-2]=-2
p=pheatmap::pheatmap(na.omit(n),fontsize = 20,show_rownames = T,border_color = "white")


# TCA

TCA <- median_path[grep("TCA",rownames(median_path)),]
rownames(TCA) <- strsplit2(rownames(TCA),": ")[,2]
View(TCA)
n=t(scale(t(TCA)))
n[n>2]=2
n[n<-2]=-2
p=pheatmap::pheatmap(na.omit(n),fontsize = 20,show_rownames = T,border_color = "white")
p
# 蛋氨酸合成









lefse_wt_apc_F.res <- read.delim("F:/A_MetaboAnalysis/lefse_wt_apc_F.res/lefse_wt_apc_F.res.txt")
lefse_wt_apc_M.res <- read.delim("F:/A_MetaboAnalysis/lefse_wt_apc_M_res/lefse_wt_apc_M.res.txt")
lefse_wt_apc_F.res$Log10AverageAbundance
ggplot(data = lefse_wt_apc_F.res, mapping = aes(x = LDA, y = Biomarkernames, fill = EnrichedGroups)) + geom_bar(stat = 'identity') +
  scale_y_discrete(limits=lefse_wt_apc_F.res$Biomarkernames)+
  theme_bw()+theme(axis.title =element_text(size = 18),
                   axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text.x = element_text(size = 18,hjust = 0.5,angle = 0))+
  theme(axis.text.y = element_text(size = 18))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),legend.text = element_text(size=12),strip.text = element_text(size=14))+
  scale_fill_manual(values=c("#DE5B0F","#0D5C98"))+
  geom_text(data = subset(lefse_wt_apc_F.res, LDA < 0),aes(y=Biomarkernames, x= 0, label= paste0(" ", Biomarkernames)),size = 6,hjust = "outward") +
  geom_text(data = subset(lefse_wt_apc_F.res, LDA > 0),aes(y=Biomarkernames, x= -0.2, label=Biomarkernames),size = 6, hjust = "outward")+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())  + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.8))+xlab("LDA_log10")+ylab("Species")

ggplot(data = lefse_wt_apc_M.res, mapping = aes(x = LDA, y = Biomarkernames, fill = EnrichedGroups)) + geom_bar(stat = 'identity') +
  scale_y_discrete(limits=lefse_wt_apc_M.res$Biomarkernames)+
  theme_bw()+theme(axis.title =element_text(size = 18),
                   axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text.x = element_text(size = 18,hjust = 0.5,angle = 0))+
  theme(axis.text.y = element_text(size = 18))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),legend.text = element_text(size=12),strip.text = element_text(size=14))+
  scale_fill_manual(values=c("#DE5B0F","#0D5C98"))+
  geom_text(data = subset(lefse_wt_apc_M.res, LDA < 0),aes(y=Biomarkernames, x= 0, label= paste0(" ", Biomarkernames)),size = 6,hjust = "outward") +
  geom_text(data = subset(lefse_wt_apc_M.res, LDA > 0),aes(y=Biomarkernames, x= -0.2, label=Biomarkernames),size = 6, hjust = "outward")+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())  + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.8))+xlab("LDA_log10")+ylab("Species")


marker_taxa<- read.delim("F:/A_MetaboAnalysis/public_marker.txt")
marker_taxa <- marker_taxa[order(marker_taxa$LDA),]
marker_taxa$scientific.name
p1<-ggplot(marker_taxa,aes(x=project.id,y=scientific.name)) #热图绘制
p2 <- p1+scale_color_gradientn(colours = c('#0D5C98','white','#DE5B0F'))+
  theme_bw()+
  geom_point(aes(size=8,
                 color=LDA))+
  theme(panel.grid = element_blank(),axis.text.x =element_text(size=16,angle =35,hjust =1))+
  theme(panel.grid = element_blank(),axis.text.y =element_text(size=16))+
  xlab(NULL) + ylab(NULL)+scale_y_discrete(limits=unique(marker_taxa$scientific.name))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+guides(size=FALSE)
p2



intersect(marker_taxa$scientific.name,gsub("_"," ",lefse_wt_apc_F.res$Biomarkernames))
intersect(marker_taxa$scientific.name,gsub("_"," ",lefse_wt_apc_M.res$Biomarkernames))
lefse_all <- read.delim("F:/A_MetaboAnalysis/lefse/lefse_wt_apc_all.txt")
lefse_all
aa <- intersect(marker_taxa$scientific.name,gsub("_"," ",lefse_all$Biomarkernames))
marker_taxa[marker_taxa$scientific.name%in%aa,]
lefse_all[lefse_all$Biomarkernames%in%gsub(" ","_",aa),]
View(lefse_all)
ggplot(data = lefse_all, mapping = aes(x = LDA, y = Biomarkernames, fill = EnrichedGroups)) + geom_bar(stat = 'identity') +
  scale_y_discrete(limits=lefse_all$Biomarkernames)+
  theme_bw()+theme(axis.title =element_text(size = 18),
                   axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text.x = element_text(size = 18,hjust = 0.5,angle = 0))+
  theme(axis.text.y = element_text(size = 18))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),legend.text = element_text(size=14),strip.text = element_text(size=14))+
  scale_fill_manual(values=c("#DE5B0F","#0D5C98"))+
  geom_text(data = subset(lefse_all, LDA < 0),aes(y=Biomarkernames, x= 0, label= paste0(" ", Biomarkernames)),size = 6,hjust = "inward") +
  geom_text(data = subset(lefse_all, LDA > 0),aes(y=Biomarkernames, x= -0.2, label=Biomarkernames),size = 6, hjust = "outward")+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())  + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.8))+xlab("LDA_log10")+ylab("Species")

diff_microb_new <- merged_abundance_table_species[rownames(merged_abundance_table_species)%in%rownames(diff_microb),]
write.table(diff_microb_new,"F:/A_MetaboAnalysis/diff_microb_new.txt", sep = '\t', col.names = NA, quote = FALSE)
View(diff_microb_new)
# sparCC
library(devtools)
devtools::install_github('zdk123/SpiecEasi')
library(SpiecEasi)
# 安装需要的包，默认不安装，没安装过的请取消如下注释
install.packages("igraph")
install.packages("psych")
# 加载包
library(igraph)
library(psych)

# 读取otu-sample矩阵，行为sample，列为otu
data= otu_order[rownames(otu_order)%in%setdiff,]


# 计算OTU间两两相关系数矩阵
# 数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
occor = corr.test(t(data),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r = occor$r # 取相关性矩阵R值
occor.p = occor$p # 取相关性矩阵p值

# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
# 构建igraph对象
igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph
# NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
# 可以按下面命令转换数据
# occor.r[occor.r!=0] = 1
# igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)

# 是否去掉孤立顶点，根据自己实验而定
# remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph

# 将igraph weight属性赋值到igraph.weight
igraph.weight = E(igraph)$weight

# 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight = NA


# 简单出图
# 设定随机种子数，后续出图都从同一随机种子数出发，保证前后出图形状相对应
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

# 如果构建网络时，weighted=NULL,此步骤不能统计
sum(igraph.weight>0)# number of postive correlation
sum(igraph.weight<0)# number of negative correlation

# set edge color，postive correlation 设定为red, negative correlation设定为blue
E.color = igraph.weight
E.color = ifelse(E.color>0, "red3",ifelse(E.color<0, "blue3","grey"))
E(igraph)$color = as.character(E.color)

# 改变edge颜色后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# 可以设定edge的宽 度set edge width，例如将相关系数与edge width关联
E(igraph)$width = abs(igraph.weight)*4

# 改变edge宽度后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


devtools::install_github('zdk123/SpiecEasi')
library(SpiecEasi)
library(ggClusterNet)
library(tidyverse)
library(phyloseq)

data(ps)

############################## micro_metab cor ######################################
data1 <- otu_order[rownames(otu_order)%in%setdiff,]
diff_ApcMF_micrometab<- read.delim("F:/A_MetaboAnalysis/Apc_metab/diff_ApcMF_micrometab20220729.txt")
colnames(diff_ApcMF_micrometab)
data2 <- diff_ApcMF_micrometab[,c(3,5:42)]
colnames(data2)
rownames(data2) <- data2$Metabolite
data2$Metabolite <- NULL
data2 <- data2[,colnames(data1)]
data1
library(psych)
library(reshape2)
n1 <- t(scale(t(data1)))
n2 <- t(scale(t(data2)))
colnames(n1)
cor <-corr.test(t(n1),t(n2), method = "spearman",adjust= "none")						
cor <-corr.test(t(data1[,22:37]),t(data2[,22:37]), method = "spearman",adjust= "none")						#cor <-corr.test(diff_microb, diff_metab, method = "pearson",adjust= "none")						
## 提取相关性、p值						
cor						
cmt <-cor$r						
View(cmt)						
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
head(df)						
#View(df)						
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
annotation_col <- data.frame(group=ApcMF[setdiff,]$group,row.names = setdiff)
#colnames(cmt)
#table(annotation_col$group)
#View(annotation_col)
#head(diff_microb)
#dim(diff_microb)
#View(cmt)
annotation_row <- data.frame(group=diff_ApcMF_micrometab$group,row.names = diff_ApcMF_micrometab$Metabolite)
anno_colors=list(group = c(Apc_F = "#CFB5C8", Apc_M = "#B0C5C3"))
p=pheatmap::pheatmap(t(cmt),fontsize = 16,show_colnames = T,legend_breaks = c(-1,-0.5,0,0.5,1),display_numbers = t(pmt),
                   color=mycol,annotation_row = annotation_row,annotation_colors = anno_colors,
                   annotation_col = annotation_col)			


	
p


setdiff

cmt						
ggdata <- t(rbind(data1[,22:37],data2[,22:37]))%>%as.data.frame(.)	

#dim(data)						
#dim(logh[,-78])						
colnames(ggdata)						
ggplot(data=ggdata, aes(x=`Akkermansia_muciniphila`, y=`LysoPC(16:1(9Z))`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					

cmt

ggplot(data=ggdata, aes(x=`Parabacteroides_goldsteinii`, y=`LysoPC(16:1(9Z))`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					
cor1=ggplot(data=ggdata, aes(x=`Akkermansia_muciniphila`, y=`L-Arginine`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					


ggplot(data=ggdata, aes(x=`Akkermansia_muciniphila`, y=`N-Acetyl-D-glucosamine`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					

ggplot(data=ggdata, aes(x=`Lactobacillus_fermentum`, y=`L-Arginine`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					

ggplot(data=ggdata, aes(x=`Lactobacillus_taiwanensis`, y=`L-Arginine`))+
  geom_point(color="black",size=5)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),
        legend.text = element_text(size=16,face = "bold"),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))					
