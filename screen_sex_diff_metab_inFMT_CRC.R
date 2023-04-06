#####################
## 代谢物的PCA分析
#####################

### 代谢组预处理：缺失值填充，以及PCA，t-test分析等等#####
setwd("/home/zhangy/project/black_rice/")
#BiocManager::install("ropls")
# library(mixOmics)
library(ropls)
library(tidyverse)

name <- read.delim2("data/metab_desc.txt")
name=name[grep("metab_",name$Metabolite,invert = T),]  ## 获取那些不是metab_开头的元素,只有1123个
intersect(name$metab_id,rownames(c))


###数据导入###
## 代谢物浓度表
org_mix<-read.table("data/metab_abund.txt",header = T,row.names = 1)

# 保留样本分组表
host<-sapply(strsplit(colnames(org_mix),split = "_|0"), function(x){
  return(x[[1]])
})
source<-sapply(strsplit(colnames(org_mix),split = "_|0"), function(x){
  return(x[[2]])
})

meta_data<-data.frame(sample = colnames(org_mix),
                      host=host,
                      source=source)
meta_data<-mutate(meta_data,Label=paste(host,source,sep = "_")) %>%
  filter(!str_detect(sample,"QC"))%>%
  column_to_rownames(var="sample")

#### 缺失值处理#####
#####将缺失值替换为数值0，方便后续我们进行数值判断
## 如果没有缺失值，就把0换成组内最小值或者是组内最小值的一半
# 计算组间最小值
Replace_Min<-function(mt_group,mt_df){
  mt_group=unique(mt_group)  ## 不同组
  mt_list=list()
  for (i in mt_group) {
    ## 不同组的矩阵
    group_df=mt_df[,str_detect(colnames(mt_df),i)]  
    kk = apply(group_df,1,function(x){
      if(sum(is.na(x))<length(x)/2){   ## 如果组内NA值小于组内样本数量一半，那么就用组内最小值填充
        min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
        # half_replace=min_replace/2  ## 也可以用最小值的一半替代
        x[is.na(x)] = min_replace   ## 该行中NA值替换成最小值 
      }else{      ## 如果组内NA值大于等于组内样本数量一半，那么就用组内最小值一半填充
        min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
        half_replace=min_replace/2  ## 也可以用最小值的一半替代
        x[is.na(x)] = half_replace   ## 该行中NA值替换成最小值,如果都是NA，则最小值是Inf  
      }
      return(x)
    })%>%t()
    
    mt_list[[i]]<-kk
  }
  return(mt_list)
}
a <- org_mix[,!str_detect(colnames(org_mix),pattern = "^QC")]
a[a==0]<-NA  ## 把0换成NA值


######缺失值过多处理#######
sum(is.na(a))

result <- data.frame("rowmname"=rownames(a), "missing"=rowSums(is.na(a)))
result$missing_percent<-result$missing/length(colnames(a)[grep("QC",colnames(a),invert = T)]) ## 34个样本
# res_filter<-filter(result,missing_percent>0.5 | missing_percent==0.5)  ## 0.5筛选
res_keep<-filter(result,missing_percent<0.5)  ## 0.5筛选

a=a[rownames(res_keep),]
sum(is.na(a))

a["metab_12706",1:10]

c<-Replace_Min(unique(meta_data$Label),a)
c<-do.call(cbind,c)
c["metab_12706",1:10]
# c[c==Inf]=0  ## 这步先不处理


######RSD筛选#########
# RSD
QC<- org_mix[,str_detect(colnames(org_mix),pattern = "^QC")]
head(QC)
RSD <- data.frame(rsd=apply(QC,1,sd)/apply(QC,1,mean))
head(RSD,35)
RSD0.2 <- subset(RSD,RSD$rsd< 0.2)  ## 因为这个数据是液相的，所以筛选值是0.2，气相的是0.3的筛选值
# RSD0.3 <- subset(RSD,RSD$rsd< 0.3)
dim(RSD0.2)
dim(c)
c=c[rownames(c) %in% rownames(RSD0.2),]
## 确实会有很多在组内代谢物浓度为0，但是QC样本RSD没有问题的组，是不是真的有这么大的差异？？
# c[c==Inf]=0  ## 这步先不处理

c<-apply(c,1,function(x){
  min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
  # half_replace=min_replace/2  ## 也可以用最小值的一半替代
  x[x==Inf] = min_replace   ## 该行中INF值用所有样本中的最小值填充
  return(x)
})%>%t()


#### log 转换#####
dim(c)
logc <- log10(c+1)   
head(logc)
name <- read.delim2("data/metab_desc.txt")
name=name[grep("metab_",name$Metabolite,invert = T),]  
intersect(name$metab_id,rownames(c))

dim(logc)
logc <- logc[rownames(logc)%in%name$metab_id,]     
dim(name)
dim(logc)

pca_count<-cbind(c,org_mix[rownames(c),str_detect(colnames(org_mix),"QC")])
logpca_count <- log10(pca_count+1)
logpca_count <- logpca_count[rownames(logpca_count)%in%name$metab_id,]   ##963

library(PCAtools)
library(patchwork)
scale_logc<- t(scale(t(logpca_count), center = T, scale = T) ) 
round(sum(scale_logc[1,]),2)  ## 要看小数点，不然sum结果不是0

plot(scale_logc, main="scaled data")

pca<-prcomp(t(scale_logc), center = F,scale. = F)  
plot(pca$x, main="after PCA")

test<-as.data.frame(pca$x[,1:2])
library(ggplot2)
test$Label=factor(substr(rownames(test),1,3),
                  levels = c("M_M","F_M","M_F","F_F","QC0"))

summary(pca) 
# Proportion of Variance  0.2584  0.1265


ggplot(test,aes(x = PC1 , y = PC2,shape=Label,fill=Label))+
  #geom_point(size = 8)+
  geom_jitter(size = 8)+
  scale_shape_manual(values = c(21,22,21,22,23))+ #自定义图形
  #scale_color_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  scale_color_manual(values=c("#51659C","#4DA2AF","#9A668C","#F1C8C8","#F67866"))+
  # scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  scale_fill_manual(values=c("#51659C","#4DA2AF","#9A668C","#F1C8C8","#F67866"))+
  theme_classic()+
  guides(colour=guide_legend(title = NULL),
         fill=guide_legend(title = NULL),
         shape=guide_legend(title = NULL))+
  theme( #调整图例
    legend.position = c(0.16,0.3),
    #legend.direction = "horizontal",
    #legend.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "transparent"),
    text = element_text(size=26,face = "bold"),
    axis.text.x=element_text(hjust=1, size=26,face = "bold"),
    axis.title = element_text(size=26,face = "bold"),
    axis.text.y = element_text(size =26,face = "bold"),
    axis.text =element_text(size = 26, color = 'black'),
    axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2),
    legend.text=element_text(size=15,face = "bold"),
    legend.title=element_text(size=15,face = "bold"))+
  xlab("PC1(25.8%)")+
  ylab("PC2(12.7%)")




#####################
## 代谢物的差异分析
#####################

### 代谢组预处理：缺失值填充，以及PCA，t-test分析等等#####
setwd("/home/zhangy/project/black_rice/")
#BiocManager::install("ropls")
# library(mixOmics)
library(ropls)
library(tidyverse)

###数据导入###
## 代谢物浓度表
org_mix<-read.table("data/metab_abund.txt",header = T,row.names = 1)

# 保留样本分组表
host<-sapply(strsplit(colnames(org_mix),split = "_|0"), function(x){
  return(x[[1]])
})
source<-sapply(strsplit(colnames(org_mix),split = "_|0"), function(x){
  return(x[[2]])
})

meta_data<-data.frame(sample = colnames(org_mix),
                      host=host,
                      source=source)
meta_data<-mutate(meta_data,Label=paste(host,source,sep = "_")) %>%
  filter(!str_detect(sample,"QC"))%>%
  column_to_rownames(var="sample")

meta_data$group=paste(meta_data$source,meta_data$host,sep = "→")
meta_data$group<-paste0("Human_",meta_data$group)

write.csv(meta_data,"data/sample_id.csv")
#### 缺失值处理#####
#####将缺失值替换为数值0，方便后续我们进行数值判断
## 如果没有缺失值，就把0换成组内最小值或者是组内最小值的一半
# 计算组间最小值
Replace_Min<-function(mt_group,mt_df){
  mt_group=unique(mt_group)  ## 不同组
  mt_list=list()
  for (i in mt_group) {
    ## 不同组的矩阵
    group_df=mt_df[,str_detect(colnames(mt_df),i)]  
    kk = apply(group_df,1,function(x){
      if(sum(is.na(x))<length(x)/2){   ## 如果组内NA值小于组内样本数量一半，那么就用组内最小值填充
        min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
        # half_replace=min_replace/2  ## 也可以用最小值的一半替代
        x[is.na(x)] = min_replace   ## 该行中NA值替换成最小值 
      }else{      ## 如果组内NA值大于等于组内样本数量一半，那么就用组内最小值一半填充
        min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
        half_replace=min_replace/2  ## 也可以用最小值的一半替代
        x[is.na(x)] = half_replace   ## 该行中NA值替换成最小值,如果都是NA，则最小值是Inf  
      }
      return(x)
    })%>%t()
    
    mt_list[[i]]<-kk
  }
  return(mt_list)
}
a <- org_mix[,!str_detect(colnames(org_mix),pattern = "^QC")]
a[a==0]<-NA  ## 把0换成NA值


######缺失值过多处理#######
sum(is.na(a))

result <- data.frame("rowmname"=rownames(a), "missing"=rowSums(is.na(a)))
result$missing_percent<-result$missing/length(colnames(a)[grep("QC",colnames(a),invert = T)]) ## 34个样本
# res_filter<-filter(result,missing_percent>0.5 | missing_percent==0.5)  ## 0.5筛选
res_keep<-filter(result,missing_percent<0.5)  ## 0.5筛选

a=a[rownames(res_keep),]
sum(is.na(a))

a["metab_12706",1:10]

c<-Replace_Min(unique(meta_data$Label),a)
c<-do.call(cbind,c)
c["metab_12706",1:10]
# c[c==Inf]=0  ## 这步先不处理


######RSD筛选#########
# RSD
QC<- org_mix[,str_detect(colnames(org_mix),pattern = "^QC")]
head(QC)
RSD <- data.frame(rsd=apply(QC,1,sd)/apply(QC,1,mean))
head(RSD,35)
RSD0.2 <- subset(RSD,RSD$rsd< 0.2)  ## 因为这个数据是液相的，所以筛选值是0.2，气相的是0.3的筛选值
# RSD0.3 <- subset(RSD,RSD$rsd< 0.3)
dim(RSD0.2)
dim(c)
c=c[rownames(c) %in% rownames(RSD0.2),]
## 确实会有很多在组内代谢物浓度为0，但是QC样本RSD没有问题的组，是不是真的有这么大的差异？？
# c[c==Inf]=0  ## 这步先不处理

c<-apply(c,1,function(x){
  min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
  # half_replace=min_replace/2  ## 也可以用最小值的一半替代
  x[x==Inf] = min_replace   ## 该行中INF值用所有样本中的最小值填充
  return(x)
})%>%t()


#### log 转换#####
dim(c)
logc <- log10(c+1)   
head(logc)
name <- read.delim2("data/metab_desc.txt")
name=name[grep("metab_",name$Metabolite,invert = T),] 
intersect(name$metab_id,rownames(c))

dim(logc)
logc <- logc[rownames(logc)%in%name$metab_id,]    
dim(name)
dim(logc)

######T-test检验和FC计算##########
# mt_df=logc
FC_Ttests_Anal<-function(groupA,groupB,mt_df){
  Pvalue<-c(rep(0,nrow(mt_df))) 
  names(Pvalue)=rownames(mt_df)
  FC<-c(rep(0,nrow(mt_df))) 
  names(FC)=rownames(mt_df)
  
  for(i in rownames(mt_df)){
    if(sd(mt_df[i,groupA])==0&&sd(mt_df[i,groupB])==0){
      Pvalue[i] <-NA
      FC[i]<-NA
    }else{
      ## test检验
      print(i)
      y=t.test(as.numeric(mt_df[i,groupA]),as.numeric(mt_df[i,groupB]))
      Pvalue[i]<-y$p.value
      ## 下面传入的是缺失值填充完之后的矩阵c计算FC
      FC[i]<-mean(as.numeric(c[i,groupA]))/mean(as.numeric(c[i,groupB]))
    }
  }
  
  res=data.frame(Pvalue=Pvalue,
                 FC=FC,
                 row.names = rownames(mt_df))
  return(res)
}
groupA=rownames(meta_data)[grep(meta_data$group,pattern = "M→")]  ## 男性来源
groupB=rownames(meta_data)[grep(meta_data$group,pattern = "F→")]  ## 女性来源

## 结果是mVSf
mVSf<-FC_Ttests_Anal(groupA=groupA,
                     groupB =groupB,
                     mt_df = logc)
# mVSf<-na.omit(mVSf)
mVSf$FDR=p.adjust(mVSf$Pvalue, "BH")
mVSf$log2FC=log2(mVSf$FC)

## 因为log10(0+1)还是=0，所以才会出现0/x的情况,所以在前面做好0的填充
## 这里的FC倍数，可以用来在后面画火山图，或者筛选差异代谢物FC>1.2|FC<0.8,但是我们统一使用vip值筛选，就用不到FC了


library(stringr)
metba_metadata<-data.frame(metba_sample=rownames(test),
                           host=str_sub(rownames(test),1,1),
                           source=str_sub(rownames(test),3,3))
metba_metadata$Label=paste(metba_metadata$host,"-FMT-H",metba_metadata$source,sep = "")
# F_FMT_HM
head(metba_metadata)
rownames(metba_metadata)<-metba_metadata$metba_sample


###OPLS-DA分析####
# 只能用在两组之间，不同组之间的比较vip值不同
# VIP value
table(name$Mode)

scale_logc<-t(scale(t(logc),center = T,scale = T))

group_1 <-factor(meta_data[rownames(t(scale_logc)),]$source)

oplsda<-ropls::opls(t(scale_logc),    ## 样本在行，代谢物在列
                    group_1, ## 分组
                    predI = 1,# OPLS设为1
                    orthoI = 1 # 正交分量数 pls=0,opls=1
                    # 当设置为NA时，执行OPLS，并使用交叉验证(最多9个正交分量)自动计算正交分量的数量。
)    ## 
library(ropls)
vipVn<-data.frame(vip=getVipVn(oplsda))
vipVn$vip_1<-if_else(vipVn$vip>1,"T","F")
table(vipVn$vip_1)

oplsda_score <- as.data.frame(oplsda@scoreMN)
# oplsda_score$group <-group_1
oplsda_score$group<-factor(meta_data[rownames((oplsda_score)),]$group,
                           levels =c("Human_M→M","Human_F→M","Human_M→F","Human_F→F") )
oplsda_score$o1 <- oplsda@orthoScoreMN[,1]

oplsda@summaryDF

library(ggrepel)
library(ggsci)
ggplot(oplsda_score, aes(p1, o1, color = group,shape=group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point(size=2) +
  #scale_shape_manual(values = c(21,21,22,22))+ #自定义图形
  scale_color_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  scale_fill_manual(values=c("#51659C","#9A668C","#B0C5C3","#CFB5C8"))+
  labs(title = "",
       x = paste0("T score[1](",round(oplsda@modelDF[1,1]*100,2),"%)"),
       y = paste0("Orthogonal T score[1](",round(oplsda@modelDF[2,1]*100,2),"%)")) +
  stat_ellipse(aes(p1, o1, fill = group), geom = "polygon",alpha = 1/4,level = 0.95)+
  # geom_label_repel(aes(label = rownames(t(scale_logc))),colour ="black")+
  theme_bw()+
  guides(colour=guide_legend(title = NULL),
         fill=guide_legend(title = NULL),
         shape=guide_legend(title = NULL))+
  theme(  
    legend.position = c(0.19,0.807),
    #         #legend.direction = "horizontal",
    legend.background = element_rect(#fill ="transparent",
      color="white"),
    text = element_text(size=26),
    axis.text.x=element_text(hjust=1, size=26),
    axis.title = element_text(size=26),
    axis.text.y = element_text(size =26),
    legend.text=element_text(size=18),
    legend.key.size = unit(1.5,"line"),  ## 图例元素之间的距离
    legend.title=element_text(size=18))
## 773 666

######差异代谢物筛选与火山图########
mVSf$vip<-vipVn[rownames(mVSf),]$vip
mVSf$type<-if_else(mVSf$Pvalue<0.05 & mVSf$vip>1,"Diff","noSig")
mVSf[which(mVSf$type=="Diff"),]$type<-ifelse(mVSf[which(mVSf$type=="Diff"),]$log2FC>0,
                                             "Male",
                                             "Female")
table(mVSf$type)
# Female   Male  noSig 
# 147    134    682 
# 
mVSf[which(mVSf$vip<1),]$type <- "noSig"
table(mVSf$type)

res_MF<- as.data.frame(mVSf)
res_MF$type <- ifelse(res_MF$Pvalue < 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')
table(res_MF$type)

res_MF[which(res_MF$vip<1),]$type <- "noSig"
table(res_MF$type)
# down noSig    up 
# 147   682   134 

res_MF$metab_id <- rownames(res_MF)
res_MF<- merge(res_MF,name,all.x=T,by="metab_id")

#########合并导出代谢物表格##########
# FC	FDR	log2FC	VIP	Metabolite	ID	KEGG.Compound.ID	Library.ID	group
res_MF_281<-res_MF[res_MF$type!="noSig",c("metab_id","FC","FDR","log2FC",
                                          "vip","Metabolite","ID",
                                          "KEGG.Compound.ID",
                                          "Library.ID","type")]

res_MF_281$group<-ifelse(res_MF_281$type=="up","Source_M","Source_F")
res_MF_281$group<-factor(res_MF_281$group,levels = c("Source_M","Source_F"))

res_MF_281<-group_by(res_MF_281,group) %>%
  arrange(group,desc(log2FC))

write.table(res_MF_281,file = "../meta_analyse/2023_output/res_MF_281.txt",sep = "\t",
            quote = F,
            row.names = T)


p=ggplot(res_MF,aes(x = log2FC,y = -log10(FDR))) +
  geom_point(aes(color = type,size=VIP),alpha = 1,size = 3) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1,
        # 标题居中
        #plot.title = element_text(),
        plot.title = element_text(hjust = 0.5,face = "bold",size = 24),
        axis.text = element_text(color = 'black',size=24,face = "bold"),
        axis.title.x  = element_text(color = 'black',size=24,face = "bold"),
        axis.title = element_text(color = 'black'),
        axis.title.y   = element_text(color = 'black',size=24,face = "bold"),
        legend.text = element_text(size=16,face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "top"
        #legend.position = c(1.3,0.5)
  ) +
  scale_color_manual(name = '',
                     # color or three types
                     # c('Female'='#DA1212','noSig'='grey','Male'='#3E7C17')
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='SourceM (n=134)','noSig'='noSig (n=682)','down'='SourceF (n=147)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Source_M vs Source_F')
p

# library(ggrepel)
# pdf("FMT/output_fig/fig/metab_volcal_2.pdf",width = 8,height = 7)
# p + geom_text_repel(data = meta,aes(x = log2FC,y = -log10(FDR),label = meta$Metabolite),
#                     force=20,
#                     color="black",
#                     size=4,point.padding = 0.5,hjust = 0.5,
#                     arrow = arrow(length = unit(0.01, "npc"), 
#                                   type = "open", ends = "last"),
#                     segment.color="black",segment.size=0.2,nudge_y=1)
# dev.off()



#########Female##############
######T-test检验和FC计算##########
# mt_df=logc
FC_Ttests_Anal<-function(groupA,groupB,mt_df){
  Pvalue<-c(rep(0,nrow(mt_df))) 
  names(Pvalue)=rownames(mt_df)
  FC<-c(rep(0,nrow(mt_df))) 
  names(FC)=rownames(mt_df)
  
  for(i in rownames(mt_df)){
    if(sd(mt_df[i,groupA])==0&&sd(mt_df[i,groupB])==0){
      Pvalue[i] <-NA
      FC[i]<-NA
    }else{
      ## test检验
      print(i)
      y=t.test(as.numeric(mt_df[i,groupA]),as.numeric(mt_df[i,groupB]))
      Pvalue[i]<-y$p.value
      ## 下面传入的是缺失值填充完之后的矩阵c计算FC
      FC[i]<-mean(as.numeric(c[i,groupA]))/mean(as.numeric(c[i,groupB]))
    }
  }
  
  res=data.frame(Pvalue=Pvalue,
                 FC=FC,
                 row.names = rownames(mt_df))
  return(res)
}

## 宿主为雌性F时
groupA=rownames(meta_data)[grep(meta_data$group,pattern = "M→F")]
groupB=rownames(meta_data)[grep(meta_data$group,pattern = "F→F")]

mVSf<-FC_Ttests_Anal(groupA=groupA,
                     groupB =groupB,
                     mt_df = logc)
# mVSf<-na.omit(mVSf)
mVSf$FDR=p.adjust(mVSf$Pvalue, "BH")
mVSf$log2FC=log2(mVSf$FC)

## 因为log10(0+1)还是=0，所以才会出现0/x的情况
## 这里的FC倍数，可以用


meta_data=meta_data[c(groupA,groupB),]
logc=logc[,rownames(meta_data)]
#######PCA看样本聚类########
#PCA（主成分分析）
library(PCAtools)
library(patchwork)
scale_logc<- t(scale(t(logc), center = T, scale = T) )
scale_logc<-na.omit(scale_logc)  ## 有一个都是0的，scale之后变成了NA值

round(sum(scale_logc[1,]),2)  ## 要看小数点，不然sum结果不是0

###OPLS-DA分析####
# 只能用在两组之间，不同组之间的比较vip值不同
# VIP value
table(name$Mode)

scale_logc<-t(scale(t(logc),center = T,scale = T))
scale_logc<-na.omit(scale_logc)  ## 有一个都是0的，scale之后变成了NA值
group_1 <-factor(meta_data[rownames(t(scale_logc)),]$source)

oplsda<-ropls::opls(t(scale_logc),    ## 样本在行，代谢物在列
                    group_1, ## 分组
                    predI = 1,# OPLS设为1
                    orthoI = 1 # 正交分量数 pls=0,opls=1
                    # 当设置为NA时，执行OPLS，并使用交叉验证(最多9个正交分量)自动计算正交分量的数量。
)    ## 
library(ropls)
vipVn<-data.frame(vip=getVipVn(oplsda))
vipVn$vip_1<-if_else(vipVn$vip>1,"T","F")
table(vipVn$vip_1)

oplsda_score <- as.data.frame(oplsda@scoreMN)
# oplsda_score$group <-group_1
oplsda_score$group<-factor(meta_data[rownames((oplsda_score)),]$group)
oplsda_score$o1 <- oplsda@orthoScoreMN[,1]

oplsda@summaryDF

######差异代谢物筛选与火山图########
mVSf$vip<-vipVn[rownames(mVSf),]$vip
mVSf$type<-if_else(mVSf$Pvalue<0.05 & mVSf$vip>1,"Diff","noSig")
mVSf[which(mVSf$type=="Diff"),]$type<-ifelse(mVSf[which(mVSf$type=="Diff"),]$log2FC>0,
                                             "Male",
                                             "Female")
table(mVSf$type)
# Female   Male  noSig 
# 65     84    814 



p2<-ggplot(mVSf,aes(x = log2FC,y = -log10(Pvalue))) +
  geom_point(aes(color = type,size=VIP),alpha = 1,size = 3)+
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1,
        # 标题居中
        #plot.title = element_text(),
        plot.title = element_text(hjust = 0.5,face = "bold",size = 24),
        axis.text = element_text(color = 'black',size=24,face = "bold"),
        axis.title.x  = element_text(color = 'black',size=24,face = "bold"),
        axis.title = element_text(color = 'black'),
        axis.title.y   = element_text(color = 'black',size=24,face = "bold"),
        legend.text = element_text(size=16,face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "top"
        #legend.position = c(1.3,0.5)
  ) +
  scale_color_manual(name = '',
                     # color or three types
                     # c('Female'='#DA1212','noSig'='grey','Male'='#3E7C17')
                     values = c('Male'='#DA1212','noSig'='grey','Female'='#3E7C17'),
                     # legend labels
                     label = c('Male'='SourceM (n=84)','noSig'='noSig (n=814)','Female'='SourceF (n=65)'))  +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Source_M vs Source_F (Female mouse)')
p2

F_res_MF<-mVSf
F_res_MF$metab_id <- rownames(F_res_MF)
F_res_MF<- merge(F_res_MF,name,all.x=T,by="metab_id")

F_res_MF_149<-F_res_MF[F_res_MF$type!="noSig",c("metab_id","FC","FDR","log2FC",
                                                "vip","Metabolite","ID",
                                                "KEGG.Compound.ID",
                                                "Library.ID","type")]

F_res_MF_149$group<-F_res_MF_149$type

F_res_MF_149<-group_by(F_res_MF_149,group) %>%
  arrange(group,desc(log2FC))

write.table(F_res_MF_149,file = "../meta_analyse/2023_output/F_res_MF_149.txt",sep = "\t",
            quote = F,
            row.names = T)
p|p2



#########宿主为Male##############
######T-test检验和FC计算##########
# mt_df=logc
FC_Ttests_Anal<-function(groupA,groupB,mt_df){
  Pvalue<-c(rep(0,nrow(mt_df))) 
  names(Pvalue)=rownames(mt_df)
  FC<-c(rep(0,nrow(mt_df))) 
  names(FC)=rownames(mt_df)
  
  for(i in rownames(mt_df)){
    if(sd(mt_df[i,groupA])==0&&sd(mt_df[i,groupB])==0){
      Pvalue[i] <-NA
      FC[i]<-NA
    }else{
      ## test检验
      print(i)
      y=t.test(as.numeric(mt_df[i,groupA]),as.numeric(mt_df[i,groupB]))
      Pvalue[i]<-y$p.value
      ## 下面传入的是缺失值填充完之后的矩阵c计算FC
      FC[i]<-mean(as.numeric(c[i,groupA]))/mean(as.numeric(c[i,groupB]))
    }
  }
  
  res=data.frame(Pvalue=Pvalue,
                 FC=FC,
                 row.names = rownames(mt_df))
  return(res)
}

## 宿主为雌性F时
groupA=rownames(meta_data)[grep(meta_data$group,pattern = "M→M")]
groupB=rownames(meta_data)[grep(meta_data$group,pattern = "F→M")]

mVSf<-FC_Ttests_Anal(groupA=groupA,
                     groupB =groupB,
                     mt_df = logc)
# mVSf<-na.omit(mVSf)
mVSf$FDR=p.adjust(mVSf$Pvalue, "BH")
mVSf$log2FC=log2(mVSf$FC)

## 因为log10(0+1)还是=0，所以才会出现0/x的情况
## 这里的FC倍数，可以用


meta_data=meta_data[c(groupA,groupB),]
logc=logc[,rownames(meta_data)]
#######PCA看样本聚类########
#PCA（主成分分析）
library(PCAtools)
library(patchwork)
scale_logc<- t(scale(t(logc), center = T, scale = T) )
scale_logc<-na.omit(scale_logc)  ## 有一个都是0的，scale之后变成了NA值

round(sum(scale_logc[1,]),2)  ## 要看小数点，不然sum结果不是0

###OPLS-DA分析####
# 只能用在两组之间，不同组之间的比较vip值不同
# VIP value
table(name$Mode)

scale_logc<-t(scale(t(logc),center = T,scale = T))
scale_logc<-na.omit(scale_logc)  ## 有一个都是0的，scale之后变成了NA值
group_1 <-factor(meta_data[rownames(t(scale_logc)),]$source)

oplsda<-ropls::opls(t(scale_logc),    ## 样本在行，代谢物在列
                    group_1, ## 分组
                    predI = 1,# OPLS设为1
                    orthoI = 1 # 正交分量数 pls=0,opls=1
                    # 当设置为NA时，执行OPLS，并使用交叉验证(最多9个正交分量)自动计算正交分量的数量。
)    ## 
library(ropls)
vipVn<-data.frame(vip=getVipVn(oplsda))
vipVn$vip_1<-if_else(vipVn$vip>1,"T","F")
table(vipVn$vip_1)

oplsda_score <- as.data.frame(oplsda@scoreMN)
# oplsda_score$group <-group_1
oplsda_score$group<-factor(meta_data[rownames((oplsda_score)),]$group)
oplsda_score$o1 <- oplsda@orthoScoreMN[,1]

oplsda@summaryDF

######差异代谢物筛选与火山图########
mVSf$vip<-vipVn[rownames(mVSf),]$vip
mVSf$type<-if_else(mVSf$Pvalue<0.05 & mVSf$vip>1,"Diff","noSig")
mVSf[which(mVSf$type=="Diff"),]$type<-ifelse(mVSf[which(mVSf$type=="Diff"),]$log2FC>0,
                                             "Male",
                                             "Female")
table(mVSf$type)
# Female   Male  noSig 
# 121    140    702



p3<-ggplot(mVSf,aes(x = log2FC,y = -log10(Pvalue))) +
  geom_point(aes(color = type,size=VIP),alpha = 1,size = 3)+
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1,
        # 标题居中
        #plot.title = element_text(),
        plot.title = element_text(hjust = 0.5,face = "bold",size = 24),
        axis.text = element_text(color = 'black',size=24,face = "bold"),
        axis.title.x  = element_text(color = 'black',size=24,face = "bold"),
        axis.title = element_text(color = 'black'),
        axis.title.y   = element_text(color = 'black',size=24,face = "bold"),
        legend.text = element_text(size=16,face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "top"
        #legend.position = c(1.3,0.5)
  ) +
  scale_color_manual(name = '',
                     # color or three types
                     # c('Female'='#DA1212','noSig'='grey','Male'='#3E7C17')
                     values = c('Male'='#DA1212','noSig'='grey','Female'='#3E7C17'),
                     # legend labels
                     label = c('Male'='SourceM (n=140)','noSig'='noSig (n=702)','Female'='SourceF (n=121)'))  +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Source_M vs Source_F (Male mouse)')
p3

M_res_MF<-mVSf
M_res_MF$metab_id <- rownames(M_res_MF)
M_res_MF<- merge(M_res_MF,name,all.x=T,by="metab_id")

M_res_MF_261<-M_res_MF[M_res_MF$type!="noSig",c("metab_id","FC","FDR","log2FC",
                                                "vip","Metabolite","ID",
                                                "KEGG.Compound.ID",
                                                "Library.ID","type")]

M_res_MF_261$group<-M_res_MF_261$type

M_res_MF_261<-group_by(M_res_MF_261,group) %>%
  arrange(group,desc(log2FC))

#write.csv(M_res_MF_261,file = "../meta_analyse/2023_output/M_res_MF_261.csv")
p|p2|p3
write.table(M_res_MF_261,file = "../meta_analyse/2023_output/res_MF_281.txt",sep = "\t",
            quote = F,
            row.names = T)

#############
M_res_MF_261
F_res_MF_149
res_MF_281

library(venn)
venn(list(M=M_res_MF_261[which(M_res_MF_261$type!="nosig"),]$metab_id,
          FMT=res_MF_281[which(res_MF_281$type!="nosig"),]$metab_id,
          F=F_res_MF_149[which(F_res_MF_149$type!="nosig"),]$metab_id),
     # zcolor=c("#B2182B","#2166AC"),
     zcolor=c("#51659C","#9A668C","#F1C8C8"),
     ilcs = 2,
     opacity = 0.85,
     sncs=2,borders = FALSE, box = FALSE)



ll<-read.csv("../meta_analyse/2023_output/res_MF_281.txt",sep = "\t")
ll1<-read.csv("../meta_analyse/2023_output/F_res_MF_149.txt",sep = "\t")

res_MF_281<-read.csv("../meta_analyse/2023_output/res_MF_281.txt",sep = "\t")
F_res_MF_149<-read.csv("../meta_analyse/2023_output/F_res_MF_149.txt",sep = "\t")
M_res_MF_261<-read.csv("../meta_analyse/2023_output/M_res_MF_261.txt",sep = "\t")

rownames(res_MF_281)<-res_MF_281$metab_id

inter_F<-res_MF_281[intersect(res_MF_281$metab_id,F_res_MF_149$metab_id),]
inter_M<-res_MF_281[intersect(res_MF_281$metab_id,M_res_MF_261$metab_id),]

inter_217_metab<-res_MF_281[unique(inter_F$metab_id,inter_M$metab_id),]
write.table(inter_217_metab,file = "../meta_analyse/2023_output/inter_217_metab.txt",sep = "\t",
            quote = F,
            row.names = T)


#####################
## 代谢物不同来源条形图和富集通路图
#####################

### 代谢物来源条形图#############
######  wt的差异代谢物和差异微生物相关性########
wtDiff_Metabolites_Origin <- read_csv("net_analyse_res/wtDiff_Metabolites_Origin.csv")
tttt<-wtDiff_Metabolites_Origin
tttt$sample="Human_metab"
table(tttt$Origin)

tttt$Origin<-ifelse(tttt$Origin=="Drug related"|tttt$Origin=="Unknown",
                    "Others",tttt$Origin)
tttt$Origin<-factor(tttt$Origin,levels = c(
  "Microbiota",
  "Co-Metabolism",
  "Host",
  "Food related",
  # "Drug related",
  # "Unknown"
  "Others"))
table(tttt$Origin)
tttt$sample<-"SourceMF"

pdf("FMT/output_fig/fig/metab_origin.pdf",width = 6,height = 6)
ggplot2::ggplot(tttt,aes(x=sample,fill=Origin)) +
  geom_bar(stat = "count",position = "fill") +
  xlab("")+ylab("% of Origin")+
  #scale_fill_manual(values  = c("#4D7587","#85AEC9","#976276","#E1AB8D","#AE6944","#C6BEB0"))+
  theme_classic()+
  theme(axis.title =element_text(size = 24,face = "bold", color = 'black'),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "right",legend.title = element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+
  #scale_y_continuous(labels = scales::percent)+
  
  scale_fill_manual(name = '',
                    values = c('Microbiota'='#4D7587','Co-Metabolism'='#85AEC9',"Host"="#A0B3A1",'Food related'='#976276','Others'='#E1AB8D'),
                    label = c('Microbiota'='Microbiota (num=3)','Co-Metabolism'='Co-Metabolism (num=6)',
                              "Host"="Host (num=2)",
                              'Food related'='Food related (num=38)','Others'='Others (num=45)'))
dev.off()




#######wt微生物和宿主代谢物富集的通路图############
wt_MPEA_Results_Microbiota <- read.csv("net_analyse_res/wt_MPEA_Results_Microbiota.csv",row.names = 1)
wt_MPEA_Results_Co_Metabolism <- read.csv("net_analyse_res/wt_MPEA_Results_Co-Metabolism.csv",row.names = 1)
wt_MPEA_Results_Host <- read.csv("net_analyse_res/MPEA_Results_Host.csv",row.names = 1)


wt_MPEA_Results<-do.call(rbind,list(wt_MPEA_Results_Microbiota,
                                    wt_MPEA_Results_Co_Metabolism,
                                    wt_MPEA_Results_Host)
)
wt_MPEA_Results<-wt_MPEA_Results[order(wt_MPEA_Results$Pvalue,decreasing = T),]
#wt_MPEA_Results$Name<-factor(wt_MPEA_Results$Name,levels = wt_MPEA_Results$Name)


wt_MPEA_Results$Metabolites
wt_MPEA_Results$micro_source<-if_else(wt_MPEA_Results$ID %in% wt_MPEA_Results_Co_Metabolism$ID,
                                      "Co_Metabolism",if_else(wt_MPEA_Results$ID %in% wt_MPEA_Results_Microbiota$ID,
                                                              "Micro_Metabolism","H_Metabolism"))


wt_MPEA_Results[wt_MPEA_Results$Name=="Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",]$Name<-"Glycosaminoglycan biosynthesis-\nchondroitin sulfate / dermatan sulfate"
wt_MPEA_Results$Name<-factor(wt_MPEA_Results$Name,levels = wt_MPEA_Results$Name)
wt_MPEA_Results$micro_source<-factor(wt_MPEA_Results$micro_source,
                                     levels =c("Micro_Metabolism","H_Metabolism","Co_Metabolism") )

pdf("FMT/output_fig/fig/metab_pathway.pdf",width = 14,height = 12)
ggplot(wt_MPEA_Results,aes(x=-log(Pvalue),y=Name,fill=micro_source))+
  #geom_col(fill="#03989E")+
  geom_col()+
  scale_fill_manual(values = c(
    "Micro_Metabolism"="#03989E",
    "H_Metabolism"="#D0CEC1",
    "Co_Metabolism"="#F4C492"))+   
  theme_classic()+
  theme(axis.title =element_text(size = 24,face = "bold", color = 'black'),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),
        axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),
        legend.position = "top",
        #legend.position = c(0.1,0.9),
        #legend.position = c(0.47,0.98),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+ylab("")
dev.off()


