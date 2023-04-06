#####################
## alpha多样性分析
#####################

species_abundance <-read_csv("FMT/data/merge34_abundance.csv")
species_abundance<-column_to_rownames(species_abundance,var="...1")
#sample_id <- read.csv("FMT/data/sample_id2.csv",row.names = 1)
sample_id <- read.csv("FMT/data/micro_sample_id.csv",row.names = 1)

sample_id<-sample_id[order(sample_id$group),]
sample_id$Label<-paste(sample_id$host,"-FMT-C",sample_id$source,sep = "")


######计算α多样性#####
## Count file for alpha diversity
alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- vegan::diversity(x, index = 'shannon',base = 2)
  Simpson <- vegan::diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  
  
  result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson, goods_Coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_Coverage)
  }
  
  
  result
}

library(reshape2)
# install.packages("vegan")
library(vegan)
library(ggplot2)
# install.packages("ggdendro")
#library(ggdendro)
# install.packages("phyloseq")
library(phyloseq)
library(tidyverse)
library(ggfun)



## 只有算阿尔法的时候才会使用count矩阵
otu_count <- read_csv("FMT/data/merge34_count.csv")
otu_count<-column_to_rownames(otu_count,var="...1")

# library()
otu_t <- t(otu_count)
alpha <- alpha_diversity(otu_t)  ## 这里的alpha_diversity是自己定义的
alpha$ACE<-NULL
alpha$goods_Coverage <- NULL
#alpha$observed_species <- NULL
#alpha$Simpson <- NULL
rownames(sample_id)<-sample_id$micro_sample
alpha$group <-sample_id[match(rownames(alpha),sample_id$micro_sample),]$Label
alpha$Source <-sample_id[match(rownames(alpha),sample_id$micro_sample),]$source


head(alpha)
otu_diversity_long <- alpha %>% gather(key = "index",value = "value", -c("group","Source"))
otu_diversity_long$value <- as.numeric(otu_diversity_long$value)

###
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


Data_summary <- summarySE(otu_diversity_long, measurevar="value", groupvars=c("group","Source","index"))
factor(otu_diversity_long$group)
otu_diversity_long$group <- factor(otu_diversity_long$group,levels=c("M-FMT-CM","F-FMT-CM","M-FMT-CF","F-FMT-CF"))
otu_diversity_long$Source <- factor(otu_diversity_long$Source,levels=c("M","F"))


pdf("FMT/output_fig/fig/shanno.pdf",width = 7,height = 7)
ggplot(otu_diversity_long[otu_diversity_long$index=="Shannon",],
       aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  #scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  scale_fill_manual(values=c("#51659C","#4DA2AF","#9A668C","#F1C8C8"))+
  geom_signif(comparisons = list(c("M-FMT-CM", "F-FMT-CM"),c("M-FMT-CM", "M-FMT-CF"),
                                 c("F-FMT-CM", "F-FMT-CF"),c("M-FMT-CF", "F-FMT-CF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Shannon diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("M-FMT-CM","F-FMT-CM","M-FMT-CF","F-FMT-CF"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),
        axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),
        legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),
        axis.line.y=element_line(size=1.2))+
  theme(legend.position = 'none',
        plot.margin = margin(t = 0.5,  # 顶部边缘距离
                             # r = 4,  # 右边边缘距离
                             # b = 4,  # 底部边缘距离
                             l = 0.2,  # 左边边缘距离
                             unit = "cm"))
dev.off()

pdf("FMT/output_fig/fig/simpson.pdf",width = 7,height = 7)
ggplot(otu_diversity_long[otu_diversity_long$index=="Simpson",],
       aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+theme_classic()+
  #scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  scale_fill_manual(values=c("#51659C","#4DA2AF","#9A668C","#F1C8C8"))+
  geom_signif(comparisons = list(c("M-FMT-CM", "F-FMT-CM"),c("M-FMT-CM", "M-FMT-CF"),
                                 c("F-FMT-CM", "F-FMT-CF"),c("M-FMT-CF", "F-FMT-CF")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Simpson diversity indices")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  scale_x_discrete(limits=c("M-FMT-CM","F-FMT-CM","M-FMT-CF","F-FMT-CF"))+	
  theme(axis.title =element_text(size = 24,face = "bold"),
        axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust =1,face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),
        legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),
        axis.line.y=element_line(size=1.2))+
  theme(legend.position = 'none',
        plot.margin = margin(t = 0.5,  # 顶部边缘距离
                             # r = 4,  # 右边边缘距离
                             # b = 4,  # 底部边缘距离
                             l = 0.2,  # 左边边缘距离
                             unit = "cm"))
dev.off()



#####################
## beta多样性分析
#####################


##### beita diversity########
# install.packages("ape")
## 物种相对丰度表格otu_order
# install.packages("tidyverse")
library(tidyverse)
library(ape)
setwd("/home/zhangy/project/black_rice/")
species_abundance <-read_csv("FMT/data/merge34_abundance.csv")
species_abundance<-column_to_rownames(species_abundance,var="...1")
sample_id <- read.csv("FMT/data/sample_id2.csv",row.names = 1)

sample_id$group<-factor(sample_id$group,levels = c("M_FMT_CM","F_FMT_CM","M_FMT_CF","F_FMT_CF"))
sample_id<-sample_id[order(sample_id$group),]

species_abundance$sum <-apply(species_abundance,1,sum)
species_abundance <- species_abundance[order(species_abundance$sum, decreasing=TRUE), ]
otu_order<-species_abundance
otu_order<-otu_order[,-ncol(otu_order)]


## my_tree一个数据库里的数据
metaphlan_tree = ape::read.tree("/NAS/panxl/CRC_human_Metagenomic_Metabolomic/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")
metaphlan_tree$tip.label <- metaphlan_tree$tip.label %>% gsub(".+\\|s__","",.)
my_tree = ape::keep.tip(metaphlan_tree,intersect(rownames(otu_order),metaphlan_tree$tip.label))

# sample的metadata
metadata<-sample_id %>% column_to_rownames(var = "micro_sample")
# metadata<-sample_id %>% column_to_rownames(var = "sample_id")

# metaphlan_anno 微生物信息表
library(limma)
metaphlan_tree2 = ape::read.tree("/NAS/panxl/CRC_human_Metagenomic_Metabolomic/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")
metaphlan_anno <- as.data.frame(strsplit2(metaphlan_tree2$tip.label,"\\|"))
metaphlan_anno$V1 <- NULL
unique(metaphlan_anno$V3)
metaphlan_anno <- apply(metaphlan_anno,2,function(x){gsub(".__","",x)})
colnames(metaphlan_anno) <- c("Kindom","Phylum","Class","Order","Family",	"Genus","Species")
metaphlan_anno <- as.data.frame(metaphlan_anno)
rownames(metaphlan_anno) <- metaphlan_anno$Species
metaphlan_anno <- metaphlan_anno[intersect(rownames(otu_order),rownames(metaphlan_anno)),]

library(phyloseq)
my_phy = phyloseq(otu_table(as.matrix(otu_order),
                            taxa_are_rows = TRUE),
                  sample_data(metadata),
                  phy_tree(my_tree),
                  tax_table(metaphlan_anno %>% as.matrix()))  #这里是否要as.matrix很严格

########

####unweight_unifrac###### 
unweight_unifrac = phyloseq::distance(my_phy,method = "unifrac") %>% as.matrix()

pcoa.au <- cmdscale(unweight_unifrac, k=3, eig=T)
points.au <- as.data.frame(pcoa.au$points)
colnames(points.au) <- c("x", "y", "z") 
eig.au <- pcoa.au$eig
points.au <-  cbind(points.au, metadata[match(rownames(points.au), rownames(metadata)), ])

### 绘pCOA图
points.au$source<-factor(points.au$source,levels = c("M","F"))
mi=c("#51659C","#9A668C")
points.au$group<-factor(points.au$group,levels = c("M_FMT_CM","F_FMT_CM","M_FMT_CF","F_FMT_CF"))
# mi=c("#51659C","#4DA2AF","#9A668C","#F1C8C8")
mi=c("#51659C","#9A668C","#51659C","#9A668C")


pdf("FMT/output_fig/fig/pcoa.pdf",width = 8.15,height = 7.8)
p2 <- ggplot(points.au,aes(x = x,y = y, fill = group,shape = group)) + 
  theme_bw() + 
  stat_ellipse(type = "t", 
               linetype = 1, 
               show.legend = FALSE,
               aes(color = group))+
  geom_point(alpha = 1, size =6) + 
  labs(x=paste("PCoA 1 (", format(100 * eig.au[1] / sum(eig.au), digits=4), "%)", sep=""), 
       y=paste("PCoA 2 (", format(100 * eig.au[2] / sum(eig.au), digits=4), "%)", sep=""),
       title = "PCoA plot of unweighted unifrac distance") + 
  theme_classic()+
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  scale_shape_manual(values = c(21,21,22,22))+ #自定义图形
  scale_fill_manual(values = mi)+
  scale_color_manual(values = mi)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( #调整图例
    legend.position = "top",
    legend.background = element_rect(fill = "transparent"),
    #legend.title  = element_blank(),
    title = element_text(size=22,face = "bold") ,
    text = element_text(size=24),
    axis.text.x=element_text(size=24,face = "bold",color = 'black'),
    axis.title = element_text(size=24,face = "bold"),
    axis.text.y = element_text(size =24,face = "bold",color = 'black'),
    axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2),
    legend.text=element_text(size=18,face = "bold"),
    #legend.key.size = unit(spaceLegend, "lines"),
    legend.title=element_text(size=0))+
  theme(#legend.position = 'none',
    plot.margin = margin(t = 0.5,  # 顶部边缘距离
                         r = 0.5,  # 右边边缘距离
                         # b = 4,  # 底部边缘距离
                         l = 0.1,  # 左边边缘距离
                         unit = "cm"))
p2
dev.off()

##########重写一个函数##########
rewrite_pairwise.adonis <- function(x,   ## 矩阵
                                    factors,     ## 分组向量
                                    sim.function = 'vegdist',     ## daisy，vegdist，distance(使用distance时只能用wunifrac，或者unifrac)
                                    sim.method = 'bray',   ## wunifrac，unifrac，bray
                                    p.adjust.m ='bonferroni',
                                    phy_obj=my_phy){
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); 
      x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else if(sim.function == 'vegdist'){
      x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)
    } else if(sim.function == 'distance'){
      library(phyloseq)
      ## 指定是phyloseq包里的distance函数
      x1 =phyloseq::distance(phy_obj,method = sim.method) %>% as.matrix()
      x1=x1[factors %in% c(co[1,elem],co[2,elem]),factors %in% c(co[1,elem],co[2,elem])]
      x1=as.dist(x1)
    }
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
} 


## 计算显著性：
Group=metadata[colnames(otu_order),]$Group
unif_pairwise <- rewrite_pairwise.adonis(x=t(otu_order), 
                                         factors=Group, 
                                         sim.function = "distance",
                                         sim.method = "unifrac",
                                         p.adjust.m = "BH",
                                         phy_obj=my_phy)

unif_pairwise

Group=metadata[colnames(otu_order),]$source
unif_pairwise <- rewrite_pairwise.adonis(x=t(otu_order), 
                                         factors=Group, 
                                         sim.function = "distance",
                                         sim.method = "unifrac",
                                         p.adjust.m = "BH",
                                         phy_obj=my_phy)
unif_pairwise


#####################
## 微生物的差异分析
#####################


species_abundance <-read_csv("FMT/data/merge34_abundance.csv")
slpecies_abundance<-column_to_rownames(species_abundance,var="...1")
sample_id <- read.csv("FMT/data/sample_id2.csv",row.names = 1)


mat<-species_abundance
FC_u_test<-function(groupA,groupB,data){
  Pvalue<-c(rep(0,nrow(data))) 
  names(Pvalue)=rownames(data)
  FC<-c(rep(0,nrow(data))) 
  names(FC)=rownames(data)
  for(i in rownames(data)){
    ## u 检验
    print(i)
    u=wilcox.test(as.numeric(data[i,groupA]),as.numeric(data[i,groupB]))
    Pvalue[i]<-u$p.value
    ## 
    FC[i]<-mean(as.numeric(data[i,groupA]))/mean(as.numeric(data[i,groupB]))
    #}
  }
  
  res=data.frame(Pvalue=Pvalue,
                 FC=FC,
                 row.names = rownames(data))
  return(res)
}

group_a=sample_id$micro_sample[grep(sample_id$Group,pattern = "_M")]#16
group_b=sample_id$micro_sample[grep(sample_id$Group,pattern = "_F")]#18

micro_fVSm_u<-FC_u_test(groupA=group_a,
                        groupB=group_b,
                        data =mat)

micro_fVSm_u$FDR=p.adjust(micro_fVSm_u$Pvalue, "BH")
micro_fVSm_u$log2FC=log2(as.numeric(micro_fVSm_u$FC))
micro_wt_fVSm_diff_u=micro_fVSm_u[micro_fVSm_u$FDR<=0.05,]#67
micro_wt_fVSm_diff_u$diff<-ifelse(micro_wt_fVSm_diff_u$log2FC>0,"Up","Down")
table(micro_wt_fVSm_diff_u$diff)
# Down   Up 
# 16   51 

# save(micro_fVSm_u,micro_wt_fVSm_diff_u,file = "FMT/wt_diff.RData")


######
## 画热图
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(-2, 0, 2), c("#336699", "white", "#CC3333"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#51659C","#4DA2AF","#9A668C","#F1C8C8")),
                       labels = c("M_M","F_M","M_F","F_F"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Group = factor(substr(colnames(mat),1,3),levels = c("M_M","F_M","M_F","F_F"))
## 看数据整体分布情况

n <- t(scale(t(mat)))
round(sum(n[1,]),2)
range(n)

pdf(file = "FMT/output_fig/Diff_micro_heatmap.pdf",width = 10,height = 13)
Heatmap(n[rownames(micro_wt_fVSm_diff_u),],name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = Group,
        show_heatmap_legend = T,
        border = "black",
        show_column_names = T,column_names_gp = gpar(fontsize=14),
        show_row_names = T,row_names_gp = gpar(fontsize=14),
        column_title = NULL,cluster_rows = T,cluster_columns = F,
        row_names_max_width = max_text_width(
          rownames(mat), 
          gp = gpar(fontsize = 14)
        ))
dev.off()


### 雄性中：
group_c=sample_id$micro_sample[grep(sample_id$Group,pattern = "M_M")]#9
group_d=sample_id$micro_sample[grep(sample_id$Group,pattern = "M_F")]#10

M_fVSm_u<-FC_u_test(groupA=group_c,
                    groupB=group_d,
                    data =mat)
M_fVSm_u$FDR=p.adjust(M_fVSm_u$Pvalue, "BH")
M_fVSm_u$log2FC=log2(as.numeric(M_fVSm_u$FC))
M_fVSm_diff_u=M_fVSm_u[M_fVSm_u$FDR<=0.05,]
M_fVSm_diff_u<-na.omit(M_fVSm_diff_u)#24
M_fVSm_diff_u$diff<-ifelse(M_fVSm_diff_u$log2FC>0,"Up","Down")
table(M_fVSm_diff_u$diff)   ## up 23  down 11

pdf(file = "FMT/output_fig/DiffM_micro_heatmap.pdf",width = 10,height = 8)
Heatmap(n[rownames(M_fVSm_diff_u),],name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = Group,
        show_heatmap_legend = T,
        border = "black",
        show_column_names = T,column_names_gp = gpar(fontsize=14),
        show_row_names = T,row_names_gp = gpar(fontsize=14),
        column_title = NULL,cluster_rows = T,cluster_columns = F,
        row_names_max_width = max_text_width(
          rownames(n[rownames(M_fVSm_diff_u),]), 
          gp = gpar(fontsize = 14)
        ))
dev.off()

### 雌性中：
group_e=sample_id$micro_sample[grep(sample_id$Group,pattern = "F_M")]#16
group_f=sample_id$micro_sample[grep(sample_id$Group,pattern = "F_F")]#12


F_fVSm_u<-FC_u_test(groupA=group_e,
                    groupB=group_f,
                    data =mat)
F_fVSm_u$FDR=p.adjust(F_fVSm_u$Pvalue, "BH")
F_fVSm_u$log2FC=log2(as.numeric(F_fVSm_u$FC))
F_fVSm_diff_u=F_fVSm_u[F_fVSm_u$FDR<=0.05,]
F_fVSm_diff_u<-na.omit(F_fVSm_diff_u)#8
F_fVSm_diff_u$diff<-ifelse(F_fVSm_diff_u$log2FC>0,"Up","Down")
table(F_fVSm_diff_u$diff)  # up 18   down 8

pdf(file = "FMT/output_fig/DiffF_micro_heatmap.pdf",width = 10,height = 8)
Heatmap(n[rownames(F_fVSm_diff_u),],name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = Group,
        show_heatmap_legend = T,
        border = "black",
        show_column_names = T,column_names_gp = gpar(fontsize=14),
        show_row_names = T,row_names_gp = gpar(fontsize=14),
        column_title = NULL,cluster_rows = T,cluster_columns = F,
        row_names_max_width = max_text_width(
          rownames(n[rownames(F_fVSm_diff_u),]), 
          gp = gpar(fontsize = 14)
        ))
dev.off()

save(micro_fVSm_u,micro_wt_fVSm_diff_u,   ##  M_M&F_M vs M_F&F_F 都是雄VS雌
     F_fVSm_u,F_fVSm_diff_u,    ##  M_Mvs M_F 都是雄VS雌
     M_fVSm_u,M_fVSm_diff_u,     ## F_M vs F_F 都是雄VS雌
     file = "/home/zhangy/project/black_rice/FMT/34_mirco_diff.RData")


########微生物的热图###########
## 用这里面的metadata表格
load("~/project/black_rice/FMT/34_mirco_diff.RData")


## 方向都是一样的
wt_and_F<-data.frame(wt=micro_wt_fVSm_diff_u[intersect(rownames(F_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     Fhost=F_fVSm_diff_u[intersect(rownames(F_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     row.names = intersect(rownames(F_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)))

wt_and_M<-data.frame(wt=micro_wt_fVSm_diff_u[intersect(rownames(M_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     Mhost=M_fVSm_diff_u[intersect(rownames(M_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     row.names = intersect(rownames(M_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)))

micro_33<-micro_wt_fVSm_diff_u[unique(c(rownames(wt_and_F),rownames(wt_and_M))),]

#######差异微生物并集的热图#############
# micro_83<-c(rownames(F_fVSm_diff_u),rownames(M_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u))
micro<-c(rownames(micro_33),
         c("Bacteroides_fragilis",
           "Bacteroides_intestinalis",
           "Bifidobacterium_animalis",
           "Bifidobacterium_breve",
           "Bifidobacterium_dentium",
           "Bifidobacterium_longum",
           "Clostridioides_difficile",
           "Clostridium_perfringens",
           "Dialister_invisus",
           "Dubosiella_newyorkensis",
           "Eubacterium_rectale",
           "Helicobacter_japonicus",
           "Klebsiella_pneumoniae",
           "Lactobacillus_salivarius",
           "Megamonas_funiformis",
           "Megamonas_hypermegale",
           "Methanobrevibacter_smithii",
           "Ruminococcus_bromii",
           "Streptococcus_salivarius",
           "Streptococcus_vestibularis",
           "Veillonella_parvula"
         ))


species_abundance <-read_csv("FMT/data/merge34_abundance.csv")
species_abundance<-column_to_rownames(species_abundance,var="...1")

species_abundance$sum <-apply(species_abundance,1,sum)
species_abundance <- species_abundance[order(species_abundance$sum, decreasing=TRUE), ]
otu_order<-species_abundance[,-ncol(species_abundance)]

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333"))
mat<-otu_order[unique(micro_83),]



top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#51659C","#4DA2AF","#9A668C","#F1C8C8")),
                       labels =c("M_M","F_M","M_F","F_F"),
                       labels_gp = gpar(col = "white", fontsize = 18)))
Group = factor(metadata[colnames(mat),]$group,levels = c("M_M","F_M","M_F","F_F"))


n <- t(scale(t(mat)))
round(sum(n[1,]),2)
range(n)

lgd1 = Legend(title = "", 
              at = c(-1,-0.5,0,0.5,1),
              col_fun = colorRamp2(c(-1, 0, 1), c("#336699", "white", "#CC3333")),
              direction = "vertical",
              labels_gp = gpar(fontsize = 14,fontface = "bold"),
              grid_width = unit(5, "mm"))
draw(lgd1, x = unit(0.94, "npc"), y = unit(0.65, "npc"))

pdf("FMT/output_fig/fig/micro_heatmap.pdf",width = 12,height = 8)
#pdf(file = "FMT/output_fig/DiffF_micro33_heatmap.pdf",width = 10,height = 8)
Heatmap(n[unique(micro),],name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = Group,
        show_heatmap_legend = T,
        border = "black",
        show_column_names = F,  ## 检查一下是否有误，再将其隐藏
        column_names_gp = gpar(fontsize=14),
        show_row_names = T,row_names_gp = gpar(fontsize=14,fontface = "bold"),
        column_title = NULL,cluster_rows = T,
        cluster_columns = F,
        #cluster_column_slices = FALSE,
        row_names_max_width = max_text_width(
          rownames(n[unique(micro),]), 
          gp = gpar(fontsize = 14)
        ))
draw(lgd1, x = unit(0.94, "npc"), y = unit(0.65, "npc"))
dev.off()



#####################################
## 差异微生物和差异代谢物相关性分析
#####################################

setwd("/home/zhangy/project/black_rice/")
# 信息表
sample_id <- read.csv("FMT/data/sample_id2.csv",row.names = 1)
#sample_id<-sample_id[sample_id$micro_sample %in% colnames(species_abundance),]


## 差异微生物 
load("~/project/black_rice/FMT/34_mirco_diff.RData")

## 微生物的相对丰度表格：
library(readr)
library(tidyverse)

species_abundance <-read_csv("~/project/black_rice/FMT/data/merge34_abundance.csv")
species_abundance<-column_to_rownames(species_abundance,var="...1")

species_abundance$sum <-apply(species_abundance,1,sum)
species_abundance <- species_abundance[order(species_abundance$sum, decreasing=TRUE), ]
otu_order<-species_abundance[,-ncol(species_abundance)]

## otu_order就是一会要分析的表格

## 差异代谢物：
library(readr)
#wt_fvsm_diff<- read.csv("output/001_fvsm_diffMetab.csv",row.names = 1)
res_MF_281<-read.csv("../meta_analyse/2023_output/res_MF_281.txt",sep = "\t")
res_MF_281[which(res_MF_281$metab_id=="metab_6784"),]$KEGG.Compound.ID="C04230"

# ## 代谢物和样本的表达矩阵：
load("~/project/black_rice/output/001_F_FvsF_M.Rdata")  
rm(oplsda,oplsda_score,vipVn,a,fVSm)  ## 保留c，因为c是填充好后的表达矩阵

######  wt的差异代谢物和差异微生物相关性########
wtDiff_Metabolites_Origin <- read_csv("net_analyse_res/2023_4_5/Metabolites_Origin.csv")

wt_id<-wtDiff_Metabolites_Origin[which(wtDiff_Metabolites_Origin$Origin %in% c("Microbiota","Co-Metabolism")),]
wt_id<-res_MF_281[which(res_MF_281$KEGG.Compound.ID %in% wt_id$KEGGID),]
rownames(wt_id)<-wt_id$metab_id
wt_metab<-c[rownames(wt_id),]
rownames(wt_metab)<-wt_id[rownames(wt_metab),]$Metabolite

sample_id$metab_sample<-str_replace_all(sample_id$metab_sample,"-","_")
wt_metab<-wt_metab[,sample_id$metab_sample]
#colnames(wt_metab)<-sample_id[which(sample_id$metab_sample %in% colnames(wt_metab)),]$micro_sample
colnames(wt_metab)<-sample_id[which(sample_id$metab_sample %in% colnames(wt_metab)),]$sample_id

#######20个微生物和代谢物相关性分析######
load("~/project/black_rice/FMT/34_mirco_diff.RData")

wt_and_F<-data.frame(wt=micro_wt_fVSm_diff_u[intersect(rownames(F_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     Fhost=F_fVSm_diff_u[intersect(rownames(F_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     row.names = intersect(rownames(F_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)))
## 5个  17个

wt_and_M<-data.frame(wt=micro_wt_fVSm_diff_u[intersect(rownames(M_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     Mhost=M_fVSm_diff_u[intersect(rownames(M_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)),]$diff,
                     row.names = intersect(rownames(M_fVSm_diff_u),rownames(micro_wt_fVSm_diff_u)))

## 26个

micro_33<-micro_wt_fVSm_diff_u[unique(c(rownames(wt_and_F),rownames(wt_and_M))),]
select_micro<-c(rownames(micro_33),c("Bacteroides_fragilis",
                                     "Bacteroides_intestinalis",
                                     "Bifidobacterium_animalis",
                                     "Bifidobacterium_breve",
                                     "Bifidobacterium_dentium",
                                     "Bifidobacterium_longum",
                                     "Clostridioides_difficile",
                                     "Clostridium_perfringens",
                                     "Dialister_invisus",
                                     "Dubosiella_newyorkensis",
                                     "Eubacterium_rectale",
                                     "Helicobacter_japonicus",
                                     "Klebsiella_pneumoniae",
                                     "Lactobacillus_salivarius",
                                     "Megamonas_funiformis",
                                     "Megamonas_hypermegale",
                                     "Methanobrevibacter_smithii",
                                     "Ruminococcus_bromii",
                                     "Streptococcus_salivarius",
                                     "Streptococcus_vestibularis",
                                     "Veillonella_parvula"
))



micro_39<-otu_order[unique(select_micro),
                    sample_id$micro_sample]

#wt_micro<-otu_order[rownames(micro_wt_fVSm_diff_u),sample_id$micro_sample]
colnames(micro_39)<-sample_id[which(sample_id$micro_sample %in% colnames(micro_39)),]$sample_id



## 相关性：
library(psych)
cor <-corr.test(t(wt_metab[,colnames(micro_39)]),t(micro_39), method = "spearman",adjust= "none")

cmt <-cor$r   ## 相关性系数
head(cmt)
library(reshape2)
pmt <- cor$p

## cor,pvalue,padj
df <-melt(cmt,value.name= "cor")  ## 180个   297
df$pvalue <- as.vector(pmt)
head(df)
df$fdr <- p.adjust(df$pvalue,method = "fdr")

pmt_adj<-pmt
pmt_adj<-matrix(p.adjust(as.numeric(pmt_adj),method = "fdr"),nrow = nrow(pmt_adj))
rownames(pmt_adj)<-rownames(pmt)
colnames(pmt_adj)<-colnames(pmt)

if(!is.null(pmt_adj)){
  ssmt <- pmt_adj< 0.01
  pmt_adj[ssmt] <- '**'
  smt <- pmt_adj > 0.01& pmt_adj < 0.05
  pmt_adj[smt] <- '*'
  pmt_adj[!ssmt&!smt]<- ''
} else{
  pmt_adj <- F
}
pmt_adj

##### 绘制相关性热图#######

##自定义颜色范围
mycol <-  colorRampPalette(c("#336699", "white", "#CC3333"))(200)


load("~/project/black_rice/output/001_wt_Fvswt_M.Rdata")  
rm(oplsda,oplsda_score,vipVn,a) 

metab_group<-res_MF_281[which(res_MF_281$Metabolite %in% c(rownames(cmt),"1-Pyrroline-4-hydroxy-2-carboxylate")),]
# metab_group[which(metab_group$Metabolite=="1-Pyrroline-4-hydroxy-2-carboxylate"),]$Metabolite<-"1-P-4-h-2-carboxylate"

annotation_col<-data.frame(Sex=if_else(metab_group$group=="Source_M","Male","Female"),
                           row.names = metab_group$Metabolite)



micro_group<-micro_wt_fVSm_diff_u[colnames(cmt),]
annotation_row <- data.frame(Sex=ifelse(micro_group$diff=="Up","Male","Female"),
                             row.names = rownames(micro_group))

head(wt_fvsm_diff)
anno_colors=list(Sex = c(Female = "#CFB5C8", Male = "#B0C5C3"))
pdf("../black_rice/net_analyse_res/2023_4_5/corr_comm_metab_micro.pdf",height = 9.2,width =13)
pheatmap::pheatmap((cmt),fontsize = 14,
                   show_colnames = T,
                   legend_breaks = c(-1,-0.5,0,0.5,1),
                   display_numbers = (pmt_adj),
                   color=mycol,
                   annotation_colors = anno_colors,
                   #annotation_legend = F,  ## 不展示分组legend
                   annotation_row = annotation_col,
                   row_names = gpar(fontsize=24 ,fontface = "bold.italic"),
                   annotation_col = annotation_row)
dev.off()
max(cmt)
