setwd("~/FMT_metab/")
org_mix <- read.delim("score50/metab_abund.txt",row.names = 1)
colnames(org_mix)
#############################################################################################################
# 保留样本分组表
Label<-sapply(strsplit(colnames(org_mix),split = "G|0"), function(x){
  return(x[[1]])
})
meta_data<-data.frame(Label=Label,row.names = colnames(org_mix))
# RSD
QC<- org_mix[,str_detect(colnames(org_mix),pattern = "^QC")]
head(QC)
RSD <- data.frame(rsd=apply(QC,1,sd)/apply(QC,1,mean))
head(RSD,35)
RSD0.3 <- subset(RSD,RSD$rsd< 0.3)
dim(RSD0.3)
# 缺失值处理
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
      min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
      # half_replace=min_replace/2  ## 也可以用最小值的一半替代
      x[is.na(x)] = min_replace   ## 该行中NA值替换成最小值,如果都是NA，则最小值是Inf
      return(x)
    })%>%t()
    
    mt_list[[i]]<-kk
  }
  return(mt_list)
}
a <- org_mix[,!str_detect(colnames(org_mix),pattern = "^QC")]
a[a==0]<-NA  ## 假设这里的0还没有处理好，还都是NA值

c<-Replace_Min(unique(meta_data$Label)[1:4],a)
c<-do.call(cbind,c)
c[c==Inf]=0



# score50
setwd("~/FMT_metab/")
library(mixOmics)
org_mix <- read.delim("score50/metab_abund.txt",row.names = 1)
colnames(org_mix)
group <- substr(colnames(org_mix),1,3)
group[39:43] <- "QC" 
Group <- group[1:38]
colnames(org_mix)
table(group)
colnames(org_mix) <- c(paste0("F_FMT_AM_",1:10),paste0("M_FMT_AM_",1:10),paste0("F_FMT_AF_",1:9),paste0("M_FMT_AF_",1:9),paste0("QC_",1:5))
dim(org_mix)
meta_data <- data.frame(group=substr(colnames(org_mix),1,8)[1:38])
rownames(meta_data) <- colnames(org_mix)[1:38]
meta_data$source <- substr(meta_data$group,8,8)
meta_data$target <- substr(meta_data$group,1,1)
meta_data$group <- gsub("_","-",meta_data$group)

# 缺失值处理
#####将缺失值替换为数值0，方便后续我们进行数值判断
str(org_mix)
a <- org_mix[,1:38]
str(a)
a[a==0] <- NA
a
dim(a)
#View(a)
result <- data.frame(rowmname=rownames(a),missing=rowSums(is.na(a)))
head(result)
result$missing_percent<-result$missing/length(colnames(a)) ## 34个样本
res_filter<-filter(result,missing_percent>0.5 | missing_percent==0.5)  ## 0.5筛选
dim(res_filter)
res_keep<-filter(result,missing_percent<0.5)  ## 0.5筛选
dim(res_keep)
c=a[rownames(a)%in%res_keep$rowmname,]
colnames(c)
#colnames(c) <- paste0("Apc_",substr(colnames(c),3,3),"→",substr(colnames(c),1,1))

# min 填充
min <- as.data.frame(t(apply(c,1,function(x){
  tapply(x,substr(colnames(c),1,8),min,na.rm=T)})))
head(min)
colnames(c)
for (i in 1:nrow(a)){
  c[i,1:10][is.na(c[i,1:10])]=min[i,2]}
for (i in 1:nrow(a)){
  c[i,11:20][is.na(c[i,11:20])]=min[i,4]}
for (i in 1:nrow(a)){
  c[i,21:29][is.na(c[i,21:29])]=min[i,1]}
for (i in 1:nrow(a)){
  c[i,30:38][is.na(c[i,30:38])]=min[i,3]}
c[c==Inf]=0
# RSD
QC<- org_mix[,39:43]
head(QC)
RSD <- data.frame(rsd=apply(QC,1,sd)/apply(QC,1,mean))
head(RSD,35)
RSD0.2 <- subset(RSD,RSD$rsd< 0.2)
dim(RSD0.2)
RSD0.2
c <- c[rownames(c)%in%rownames(RSD0.2),]
dim(c)
# log 转换
logc <- log10(c+1)
head(logc)

#View(logc)
name <- read.delim("score50/metab_desc.txt",encoding = "utf-8")
rownames(names) <- name$X
name=name[grep("metab_",name$Metabolite,invert = T),]
#View(name)
dim(name)
name <- name[,c(1:3,12:13)]
head(name)
#write.table(name,"score50/name.txt", sep = '\t', col.names = NA, quote = FALSE)
#BiocManager::install("ropls")

library(ropls)
dim(logc)
intersect(name$metab_id,rownames(logc))
logc <- logc[rownames(logc)%in%name$metab_id,]
logc
colnames(logc)
#View(logc)
#pheatmap(cor(logc))
dim(name)
dim(logc)
dim(c)
#logc$`Apc_F→F_4` <- NULL
#logc$`Apc_F→F_7` <- NULL
#logc$`Apc_M→M_4` <- NULL
e <- c[rownames(c)%in%rownames(logc),]
colnames(logc)
##########################source M vs sourceF##############################
# pvalue FC
dim(logc)
dim(name)
logc <- logc[rownames(logc)%in%name$metab_id,]
colnames(logc)

Pvalue<-c(rep(0,nrow(logc))) 
FC<-c(rep(0,nrow(logc))) 

for(i in 1:nrow(logc)){
  if(sd(logc[i,1:20])==0&&sd(logc[i,21:38])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,1:20]),as.numeric(logc[i,21:38]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(e[i,1:20]))/mean(as.numeric(e[i,21:38])) 
    }
}

colnames(logc)

#for(i in 1:nrow(logc)){
#  if(sd(logc[i,1:19])==0&&sd(logc[i,20:35])==0){
#    Pvalue[i] <-"NA"
#    FC[i]<-"NA"
#  }else{
#    y=t.test(as.numeric(logc[i,1:19]),as.numeric(logc[i,20:35]))
#    Pvalue[i]<-y$p.value
#    FC[i]<-mean(as.numeric(e[i,1:19]))/mean(as.numeric(e[i,20:35])) 
#  }
#}

FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)
head(mix)

# z-score
dim(mix)
head(mix)
logc_scale <- t(scale(t(mix[,1:38])))
dim(logc_scale)
#PCA（主成分分析）
#BiocManager::install("PCAtools")

library(PCAtools)
#meta_data$group <- substr(colnames(logc),1,7)
rownames(meta_data) <- colnames(c)
length(colnames(logc_scale))
dim(meta_data)

colnames(logc)==rownames(meta_data)
pca <- pca(logc, metadata=meta_data)

biplot(pca, x = 'PC1', y = 'PC2') #可以看到这两个成分对样品的解释度
#screen(p) #可以看到所有成分对样品的解释
##将pca与关联样本信息表 
pca$variance
pca_rotated_plus <- cbind(pca$rotated,meta_data)
library(ggsci)
pca_rotated_plus$group <- factor(pca_rotated_plus$group,levels=c("M-FMT-AM","F-FMT-AM","M-FMT-AF","F-FMT-AF"))
ggplot(pca_rotated_plus,aes(x = PC1 , y = PC2))+
  geom_point(size = 8,aes(shape = group, fill = group)) +
  stat_ellipse(aes(color = group,
                   fill=pca_rotated_plus$group),
               linetype = 'dashed',
               size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'PC1 (32.22%)',y = 'PC2 (12.90%)') + 
  scale_shape_manual(values = c(21,22,23,24))+ 
  scale_fill_manual(values = c("#194F90","#4E8DC3","#C5000C","#EE7762"))+
  scale_color_manual(values = c("#194F90","#4E8DC3","#C5000C","#EE7762"))+
  theme_classic()+theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "top",legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))


# VIP value
colnames(logc_scale)
pos <- logc_scale[rownames(logc_scale)%in%name[grep("pos",name$ID),]$metab_id,]
neg <- logc_scale[rownames(logc_scale)%in%name[grep("neg",name$ID),]$metab_id,]
dim(pos)
dim(neg)
l <- list(pos,neg)
colnames(l[[1]])
dat <- list()
group <- list()
data <- list()
out <- list()
meta_data$source <- substr(meta_data$group,8,8)
meta_data$target <- substr(meta_data$group,1,1)
meta_data
colnames(pos)
oplsda <- purrr::map(1:2,function(x){
  dat[[x]] <- as.data.frame(t(l[[x]]),stringsAsFactors = F)
  group[[x]] <-substr(colnames(pos),8,8) %>% factor()
  data[[x]]= as.data.frame(apply(dat[[x]],2,as.numeric))
  rownames(data[[x]]) <- rownames(dat[[x]])
  out[[x]] <- opls(data[[x]], group[[x]], orthoI = 1) 
})
vipVn <- list()
vipVn_select <- list()
vip<- do.call(rbind,purrr::map(1:2,function(x){
  vipVn[[x]] <- getVipVn(oplsda[[x]])%>%as.data.frame(.)
}))
dim(vip)
mix <- cbind(mix,vip)
head(vip)
dim(mix)
mix$log2FC <- log2(mix$FC)
mix_select1=subset(mix,FDR<0.05 & .>1 )
dim(mix_select1)
#View(mix_select1)
res_MF<- as.data.frame(mix)
res_MF$type <- ifelse(res_MF$FDR < 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

res_MF[which(res_MF$.<1),]$type <- "noSig"
res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
# check
head(res_MF,5)
res_MF$metab_id <- rownames(res_MF)
res_MF<- merge(res_MF,name,all.x=T,by="metab_id")
#View(res_MF)
# filter gene

length(Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name)
Metabolites_Origin <- read.csv("score50/source_MetOrigin_20220805160748/02_Origin_Analysis/Metabolites_Origin.csv")

meta <- res_MF[res_MF$Metabolite%in%c(Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name,
               "Oxoadipic acid",
               "Folinic acid",
               "Deoxycholic acid 3-glucuronide",
               "4-Hydroxycyclohexylcarboxylic acid",
               "Hippuric acid",
               "Chlorogenic Acid",
               "N-Acetyl-L-glutamic acid",
               "2-Hydroxyphenylacetic Acid",
               "5-(hydroxymethyl)-2-Furancarboxylic acid",
               "5'-Deoxy-5'-(methylthio)adenosine",
               "N-Acetyl-L-glutamate 5-semialdehyde",
               "N-Acetylneuraminic acid",
               "Ansamitocin P3",
               "Indolelactic acid",
               "Tetrahydrofolyl-[Glu](2)",
               "Sinapic acid",
               "7,8-Dihydropteroic acid",
               "5-Aminoimidazole ribonucleotide"),]
dim(meta)
meta <- meta[meta$Metabolite%in%c("Sinapic acid","Levan","Indole","ChlorogenicAcid",
                                  "Biliverdin","DTMP","Folinic acid"),]
head(res_MF)
dim(meta)
p=ggplot(res_MF,aes(x = log2FC,y = -log10(FDR))) +
  geom_point(aes(color = type,size=VIP),alpha = 0.5,size = 3) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1,
        # 标题居中
        plot.title = element_text(),
        axis.text = element_text(color = 'black',size=24),axis.title.x  = element_text(color = 'black',size=24),
        axis.title = element_text(color = 'black'),axis.title.y   = element_text(color = 'black',size=24),
        legend.text = element_text(size=14),legend.position = "right") +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='up (SourceM, num=162)','noSig'='noSig (FDR>0.05 & VIP<1, num=728)','down'='down (SourceF, num=219)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Source_M vs Source_F')
p + geom_text_repel(data = meta,aes(x = log2FC,y = -log10(FDR),label = meta$Metabolite),
                      force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                      arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                      segment.color="black",segment.size=0.2,nudge_y=1)+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))

# fdr vip mix
mix_select1=subset(mix,FDR<0.05 & .>1 )
# fc
#mix_select2=subset(mix_select1,FC< 0.7| FC >1.5 )
subset(mix_select1,FC< 1)%>%dim()
subset(mix_select1,FC> 1)%>%dim()

mix_select1$metab_id <- rownames(mix_select1)
mix_select <- merge(mix_select1,name,all.x=T,by="metab_id")
mix_select <- mix_select[order(mix_select$FC,decreasing = T),]

#View(mix_select)
mix_select$group <- c(rep("Source_M",162),rep("Source_F",219))
View(mix_select)

write.table(mix_select,"score50/Source_Mvs_F_20220805.txt", sep = '\t', col.names = NA, quote = FALSE)

cc <- logc
cc$metab_id <- rownames(cc)
colnames(mix_select)
diff_sourceMF_metab <- merge(mix_select[,c("metab_id","group","Metabolite")],cc,all.x=T,by="metab_id")
dim(diff_sourceMF_metab)

#View(Metabolites_Origin)

diff_sourceMF_micrometab <- diff_sourceMF_metab[diff_sourceMF_metab$Metabolite%in%c(
  Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name,
  "Oxoadipic acid",
  "Folinic acid",
  "Deoxycholic acid 3-glucuronide",
  "4-Hydroxycyclohexylcarboxylic acid",
  "Hippuric acid",
  "Chlorogenic Acid",
  "N-Acetyl-L-glutamic acid",
  "2-Hydroxyphenylacetic Acid",
  "5-(hydroxymethyl)-2-Furancarboxylic acid",
  #"5'-Deoxy-5'-(methylthio)adenosine",
  "N-Acetyl-L-glutamate 5-semialdehyde",
  "N-Acetyl-L-glutamicacid",
  "N-Acetylneuraminic acid",
  "Tyramine glucuronide",
  "Ansamitocin P3",
  "Indolelactic acid",
  "Tetrahydrofolyl-[Glu](2)",
  "Sinapic acid",
  "7,8-Dihydropteroic acid",
  "5-Aminoimidazole ribonucleotide"),]
dim(diff_sourceMF_micrometab)
write.table(diff_sourceMF_micrometab,"score50/diff_sourceMF_micrometab20220805.txt", sep = '\t', col.names = NA, quote = FALSE)


Source_Mvs_F_select <- mix_select
dim(Source_Mvs_F_select)
Apc_MvsApc_F_select <- read.delim("~/Apc_metab/score50/Apc_MvsApc_F_select20220725.txt")
Source_Mvs_F_select <- read.delim("~/FMT_metab/score50/Source_Mvs_F_select20220725.txt")
intersect(Source_Mvs_F_select[Source_Mvs_F_select$group=="Source_M",]$Metabolite,Apc_MvsApc_F_select[Apc_MvsApc_F_select$group=="Apc_M",]$Metabolite)

intersect(Source_Mvs_F_select[Source_Mvs_F_select$group=="Source_F",]$Metabolite,Apc_MvsApc_F_select[Apc_MvsApc_F_select$group=="Apc_F",]$Metabolite)
intersect(Source_Mvs_F_select[Source_Mvs_F_select$group=="Source_M",]$Metabolite,Apc_MvsApc_F_select[Apc_MvsApc_F_select$group=="Apc_M",]$Metabolite)
#intersect(M_MvsM_F_select[M_MvsM_F_select$group=="Source_F",]$Metabolite,M_MvsM_F_select[M_MvsM_F_select$group=="Apc_F",]$Metabolite)
#intersect(M_MvsM_F_select[M_MvsM_F_select$group=="Source_M",]$Metabolite,M_MvsM_F_select[M_MvsM_F_select$group=="Apc_M",]$Metabolite)


#  M_M vs M_F
colnames(logc)
d <- logc[,c(11:20,30:38)]
e <- c[,c(11:20,30:38)][rownames(c[,c(11:20,30:38)])%in%rownames(logc),]

colnames(e)
dim(e)
Pvalue<-c(rep(0,nrow(d))) 
FC<-c(rep(0,nrow(d))) 
colnames(d)

for(i in 1:nrow(d)){
  if(sd(d[i,1:10])==0&&sd(d[i,11:19])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(d[i,1:10]),as.numeric(d[i,11:19]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(e[i,1:10]))/mean(as.numeric(e[i,11:19])) 
  }
}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(d,FC)%>%cbind(.,FDR)
colnames(pos)
l <- list(pos[,c(11:20,30:38)],neg[,c(11:20,30:38)])
dat <- list()
group <- list()
data <- list()
out <- list()
oplsda <- purrr::map(1:2,function(x){
  dat[[x]] <- as.data.frame(t(l[[x]]),stringsAsFactors = F)
  group[[x]] <- substr(colnames(d),8,8) %>% factor()
  data[[x]]= as.data.frame(apply(dat[[x]],2,as.numeric))
  rownames(data[[x]]) <- rownames(dat[[x]])
  out[[x]] <- opls(data[[x]], group[[x]], orthoI = 1) 
})
vipVn <- list()
vipVn_select <- list()
vip<- do.call(rbind,purrr::map(1:2,function(x){
  vipVn[[x]] <- getVipVn(oplsda[[x]])%>%as.data.frame(.)
}))
head(vip)
mix <- cbind(mix,vip)
#View(mix)
# fdr vipmix
mix_select1=subset(mix,FDR<0.05 & .>1 )
# fc
#mix_select2=subset(mix_select1,FC< 0.7| FC >1.5 )
subset(mix_select1,FC< 1)%>%dim()
subset(mix_select1,FC> 1)%>%dim()

dim(mix_select)
mix_select1$metab_id <- rownames(mix_select1)
mix_select <- merge(mix_select1,name,all.x=T,by="metab_id")
mix_select <- mix_select[order(mix_select$FC,decreasing = T),]
#View(mix_select)
mix_select$group <- c(rep("Source_M",73),rep("Source_F",108))
#View(mix_select)
#write.table(mix_select,"score50/M_MvsM_F_select20220805.txt", sep = '\t', col.names = NA, quote = FALSE)

dim(mix_select)
M_MvsM_F_select <- mix_select
dim(mix)
mix$log2FC <- log2(mix$FC)
res_MF<- subset(mix)
res_MF$type <- ifelse(res_MF$FDR < 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

res_MF[which(res_MF$.<1),]$type <- "noSig"
res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
# check
head(res_MF,5)
res_MF$metab_id <- rownames(res_MF)
res_MF<- merge(res_MF,name,all.x=T,by="metab_id")
#View(res_MF)

# filter gene
Metabolites_Origin <- read.csv("~/FMT_metab/score50/MM_FM_MetOrigin_20220730212209/02_Origin_Analysis/Metabolites_Origin.csv")
View(res_MF)
View(Metabolites_Origin)
meta <- res_MF[res_MF$Metabolite%in%c(Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name,
                                      "Deoxycholic acid 3-glucuronide","Cholic acid","5-Aminoimidazole ribonucleotide","Methyl n-acetylanthranilate") ,]
dim(meta)
p=ggplot(res_MF,aes(x = log2FC,y = -log10(FDR))) +
  geom_point(aes(color = type,size=VIP),alpha = 0.5,size = 3) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1,
        # 标题居中
        plot.title = element_text(),
        axis.text = element_text(color = 'black',size=24),axis.title.x  = element_text(color = 'black',size=24),
        axis.title = element_text(color = 'black'),axis.title.y   = element_text(color = 'black',size=24),
        legend.text = element_text(size=14),legend.position = "right") +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='up (M→M, num=65)','noSig'='noSig (FDR>0.05 & VIP<1, num=942)','down'='down (F→M, num=102)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('M→M vs F→M')
p + geom_text_repel(data = meta,aes(x = log2FC,y = -log10(FDR),label = meta$Metabolite),
                    force=20,color="black",size=4,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)


Metabolites_Origin <- read.csv("score50/MM_FM_MetOrigin_20220730212209/02_Origin_Analysis/Metabolites_Origin.csv")
#Metabolites_Origin <- read.csv("score50/MetOrigin_20220522192255/02_Origin_Analysis/Metabolites_Origin.csv")

View(Metabolites_Origin)
#Metabolites_Origin$diff <- gsub(1,"SigDiff",Metabolites_Origin$Diff)%>%gsub("0","UnSigDiff",.)
dim(Metabolites_Origin)
Metabolites_Origin$value <- 1
Metabolites_Origin$Origin <- gsub("Drug related", "Others",Metabolites_Origin$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin$Origin <- factor(Metabolites_Origin$Origin,levels = c("Microbiota","Co-Metabolism","Host","Food related","Others"  ))
unique(Metabolites_Origin$Origin)
Metabolites_Origin$sample <- "SourceMF"

ggplot(Metabolites_Origin) +
  geom_bar(aes(x =sample,y=value,fill = Origin),stat = "identity",position = "fill") +
  scale_fill_manual(values  = c("#4D7587","#85AEC9","#A0B3A1","#976276","#E1AB8D"))+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text = element_text(size = 20),
        legend.text = element_text(size=14),legend.title  = element_text(size=14),legend.position = "right")+
  xlab("")+ylab("% of Origin")


ggdata = list(c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism")]),c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Host","Co-Metabolism")]))
names(ggdata) <- c("Microbiota","Host")
ggdata
#BiocManager::install("ggvenn")
#library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 8,
         stroke_color = "white",
         fill_color = c("#85AEC9","#A0B3A1"),
         set_name_color = c("#85AEC9","#A0B3A1")) 
p

MPEA_Results_Co.Metabolism <- read.csv("score50/MM_FM_MetOrigin_20220730212209//03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
MPEA_Results_Mico.Metabolism <- read.csv("score50/MM_FM_MetOrigin_20220730212209/03_Function_Analysis/MPEA_Results_Microbiota.csv")
MPEA_Results_H.Metabolism <- read.csv("score50/MM_FM_MetOrigin_20220730212209//03_Function_Analysis/MPEA_Results_Host.csv")
MPEA_Results_Co.Metabolism$group <- "Co_Metabolism"
MPEA_Results_Mico.Metabolism$group <- "Micro_Metabolism"
MPEA_Results_H.Metabolism$group <- "H_Metabolism"

MPEA_Results.Metabolism <- rbind(MPEA_Results_Co.Metabolism,MPEA_Results_Mico.Metabolism)%>%rbind(.,MPEA_Results_H.Metabolism)
MPEA_Results.Metabolism$logPvalue <- -log(MPEA_Results.Metabolism$Pvalue)
data <- MPEA_Results.Metabolism[order(MPEA_Results.Metabolism$Pvalue,decreasing = T),]
data$group <- factor(data$group,levels = c("Micro_Metabolism","H_Metabolism","Co_Metabolism"))
ggplot(data, aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#03989E","#D0CEC1","#FF9800"))+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=12),legend.position = "right")+
  ylab("")+xlab("log Pvalue")+
  #geom_text(aes(label=Sig..Metabolites),hjust=1,colour="white")+
  scale_y_discrete(limits=unique(data$Name))













Source_Mvs_F_select <- read.delim("~/FMT_metab/score50/Source_Mvs_F_select.txt")
M_MvsM_F_select <- read.delim("~/FMT_metab/score50/M_MvsM_F_select.txt")
Apc_MvsApc_F_select <- read.delim("~/Apc_metab/score50/Apc_MvsApc_F_select.txt")
#intersect(Source_Mvs_F_select[Source_Mvs_F_select$group=="Source_M",]$Metabolite,M_MvsM_F_select[M_MvsM_F_select$group=="Source_M",]$Metabolite)
#intersect(Source_Mvs_F_select[Source_Mvs_F_select$group=="Source_F",]$Metabolite,M_MvsM_F_select[M_MvsM_F_select$group=="Source_F",]$Metabolite)
intersect(Apc_MvsApc_F_select[Apc_MvsApc_F_select$group=="Apc_M",]$Metabolite,mix_select[mix_select$group=="Source_M",]$Metabolite)
intersect(Apc_MvsApc_F_select[Apc_MvsApc_F_select$group=="Apc_F",]$Metabolite,mix_select[mix_select$group=="Source_F",]$Metabolite)

intersect(Metabolites_Origin[Metabolites_Origin$Origin%in%c("Co-Metabolism","Microbiota"),]$Name,mix_select[mix_select$group=="Source_M",]$Metabolite)

data1 <- mix_select[mix_select$Metabolite%in%c("LysoPE(16:1(9Z)/0:0)","LysoPC(14:0/0:0)","L-Glutamine","Dehydrosoyasaponin I"),][,2:20]
rownames(data1) <- c("LysoPE(16:1(9Z)/0:0)","LysoPC(14:0/0:0)","L-Glutamine","Dehydrosoyasaponin I")
View(mix_select)

data1 <- M_MvsM_F_select[M_MvsM_F_select$Metabolite%in%unique(c(ggdata$Microbiota,
                                                         "5-Aminoimidazole ribonucleotide","Cholic acid","Methyl n-acetylanthranilate",
                                                         "7,8-Dihydropteroic acid","Sinapic acid","Ansamitocin P3")),]

rownames(data1) <- data1$Metabolite
data1 <- data1[,2:20]

New_FMT <- read.delim("~/FMT_metab/score50/New_FMT.txt")
New_FMT <- na.omit(New_FMT)
rownames(New_FMT) <- New_FMT$X
New_FMT$X <- NULL
data2 <- data1[,colnames(data1)%in%colnames(New_FMT)]
colnames(data2)
data3 <- New_FMT[,colnames(data2)]


cor <-corr.test(t(data2),t(data3), method = "spearman",adjust= "none")
#cor <-corr.test(diff_microb, diff_metab, method = "pearson",adjust= "none")
## 提取相关性、p值
cor
cmt <-cor$r
head(cmt)
c=melt(cmt)
pmt <- cor$p
#p=melt(pmt)
#cp <- cbind(c,p[,3])

## 输出相关系数表格,第一行为代谢物信息，第一列为物种信息
cmt.out<-cbind(rownames(cmt),cmt)
#write.table(cmt.out,file= "cor.txt",sep= "t",row.names=F)
## 输出p值表格，第一行为代谢物信息，第一列为物种信息
pmt.out<-cbind(rownames(pmt),pmt)
#write.table(pmt.out,file= "pvalue.txt",sep= "t",row.names=F)
## 第一列为物种名，第二列为代谢物名，第三、第四列对应显示相关系数与p值
df <-melt(cmt,value.name= "cor")
df$pvalue <- as.vector(pmt)
head(df)
df$fdr <- p.adjust(df$pvalue,method = "fdr")
df <- subset(df,abs(cor)>0.6& fdr<0.05)
#View(df)
write.table(df,file= "cor-p.txt",sep= "t")
# 绘制显著性标记
## 对所有p值进行判断，p< 0.01的以“*”标注，p值 0.01<p< 0.05的以“”标注
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


# 绘制相关性热图

##自定义颜色范围
mycol <-  colorRampPalette(c("#336699", "white", "#CC3333"))(200)

#绘制热图,可根据个人需求调整对应参数
#scale=”none” 不对数据进行均一化处理 可选 "row", "column"对行、列数据进行均一化
#cluster_row/col=T 对行或列数据进行聚类处理，可选F为不聚类
#border=NA 各自边框是否显示、颜色，可选“white”等增加边框颜色
#number_color=”white” 格子填入的显著性标记颜色
#cellwidth/height=12 格子宽度、高度信息
head(diff_metab)
dim(diff_metab)
data2
#annotation_col <- data.frame(group=diff_metab[diff_metab$Metabolite%in%rownames(data2),][,27],row.names = rownames(diff_metab[diff_metab$Metabolite%in%rownames(data2),]))
colnames(cmt)

table(annotation_col$group)
View(annotation_col)
head(diff_microb)
dim(diff_microb)
View(cmt)
#annotation_row <- data.frame(group=diff_microb[,20],row.names = rownames(diff_microb))

#anno_colors=list(group = c(Apc_F = "#CFB5C8", Apc_M = "#B0C5C3"))
pheatmap::pheatmap(t(cmt),fontsize = 20,show_colnames = T,legend_breaks = c(-1,-0.5,0,0.5,1),display_numbers = t(pmt),
                   color=mycol)

























data1 <- mix_select[mix_select$Metabolite%in%ggdata$Microbiota,][,2:20]
rownames(data1) <- c("LysoPE(16:1(9Z)/0:0)","LysoPC(14:0/0:0)","L-Glutamine","Dehydrosoyasaponin I")
View(mix_select)


New_FMT <- read.delim("~/FMT_metab/score50/New_FMT.txt",col.names = 1)

New_FMT <- na.omit(New_FMT)
rownames(New_FMT) <- New_FMT$X
New_FMT$X <- NULL
data2 <- data1[,colnames(data1)%in%colnames(New_FMT)]
colnames(data2)
data3 <- New_FMT[,colnames(data2)]
data <- rbind(data2,data3) 


colnames(pos)

colnames(neg)




# F_MvsF_F
colnames(logc)
d <- logc[,c(1:10,21:29)]
e <- c[,c(1:10,21:29)]
colnames(e)
Pvalue<-c(rep(0,nrow(d))) 
FC<-c(rep(0,nrow(d))) 
for(i in 1:nrow(d)){
  if(sd(d[i,1:10])==0&&sd(d[i,11:19])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(d[i,1:10]),as.numeric(d[i,11:19]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(e[i,1:10]))/mean(as.numeric(e[i,11:19])) 
  }
}

FDR=p.adjust(Pvalue, "BH")
mix <- cbind(d,FC)%>%cbind(.,FDR)
l <- list(pos[,c(1:10,21:29)],neg[,c(1:10,21:29)])
dat <- list()
group <- list()
data <- list()
out <- list()
oplsda <- purrr::map(1:2,function(x){
  dat[[x]] <- as.data.frame(t(l[[x]]),stringsAsFactors = F)
  group[[x]] <- substr(colnames(d),8,8) %>% factor()
  data[[x]]= as.data.frame(apply(dat[[x]],2,as.numeric))
  rownames(data[[x]]) <- rownames(dat[[x]])
  out[[x]] <- opls(data[[x]], group[[x]], orthoI = 1) 
})
vipVn <- list()
vipVn_select <- list()
vip<- do.call(rbind,purrr::map(1:2,function(x){
  vipVn[[x]] <- getVipVn(oplsda[[x]])%>%as.data.frame(.)
}))
dim(vip)
dim(mix)
mix <- cbind(mix[rownames(mix)!="metab_20800",],vip)
#View(mix)
# fdr vipmix
mix_select1=subset(mix,FDR<0.05 & .>1 )
# fc
#mix_select2=subset(mix_select1,FC< 0.7| FC >1.5 )
subset(mix_select1,FC< 1)%>%dim()
subset(mix_select1,FC> 1)%>%dim()
dim(mix_select)
mix_select1$metab_id <- rownames(mix_select1)
mix_select <- merge(mix_select1,name,all.x=T,by="metab_id")
mix_select <- mix_select[order(mix_select$FC,decreasing = T),]
#View(mix_select)
mix_select$group <- c(rep("Source_M",199),rep("Source_F",218))
write.table(mix_select,"score50/F_MvsF_F_select20220805.txt", sep = '\t', col.names = NA, quote = FALSE)
M_MvsM_F <- read.delim("~/FMT_metab/score50/M_MvsM_F_select20220805.txt")
F_MvsF_F<- read.delim("~/FMT_metab/score50/F_MvsF_F_select20220805.txt")
Source_MvsF <- read.delim("~/FMT_metab/score50/Source_Mvs_F_20220805.txt")
intersect(Source_MvsF$Metabolite,c(M_MvsM_F$Metabolite,F_MvsF_F$Metabolite))
intersect(M_MvsM_F$Metabolite,F_MvsF_F$Metabolite)
a=intersect(M_MvsM_F$Metabolite,F_MvsF_F$Metabolite)%>%intersect(.,Source_MvsF$Metabolite)

all_diff=name[name$Metabolite%in%c(Source_MvsF$Metabolite,M_MvsM_F$Metabolite,F_MvsF_F$Metabolite),]
dim(all_diff)
write.table(all_diff,"score50/all_diff20220805.txt", sep = '\t', col.names = NA, quote = FALSE)
ggdata = list(M_MvsM_F$Metabolite,F_MvsF_F$Metabolite,Source_MvsF$Metabolite)
names(ggdata) <- c("M-FMT(MvsF)","F-FMT(MvsF)","FMT-MvsFMT-F")
ggdata
#BiocManager::install("ggvenn")
#library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 10,
         stroke_color = "black",
         fill_color = c("#85AEC9","#A0B3A1","#B2662A"),
         set_name_color = c("black","black","black")) 
p






## fece
# 载入 XML 包
gc()
library("XML")
xmldataframe <- xmlToDataFrame("../feces_metabolites/feces_metabolites.xml")
#print(xmldataframe)
#View(xmldataframe)
CRC <- xmldataframe[grep("Colorectal",xmldataframe$diseases,ignore.case=T),]
dim(CRC)
View(CRC)
CRC_feces <- CRC[grep("feces",CRC$abnormal_concentrations,ignore.case=T),]
write.table(CRC_feces,"../feces_metabolites/CRC_feces.txt", sep = '\t', col.names = NA, quote = FALSE)
Source_Mvs_F_select <- read.delim("~/FMT_metab/score50/Source_Mvs_F_select.txt")
View()
intersect(CRC_feces$traditional_iupac,Source_Mvs_F_select$Metabolite)










#metab_origin
Metabolites_Origin <- read.csv("score50/source_MetOrigin_20220805160748//02_Origin_Analysis/Metabolites_Origin.csv")
#Metabolites_Origin <- read.csv("score50/source_MFMetOrigin_20220804220343/02_Origin_Analysis/Metabolites_Origin.csv")

View(Metabolites_Origin)
#Metabolites_Origin$diff <- gsub(1,"SigDiff",Metabolites_Origin$Diff)%>%gsub("0","UnSigDiff",.)
dim(Metabolites_Origin)
Metabolites_Origin$value <- 1
Metabolites_Origin$Origin <- gsub("Drug related", "Others",Metabolites_Origin$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin$Origin <- factor(Metabolites_Origin$Origin,levels = c("Microbiota","Co-Metabolism","Host","Food related","Others"  ))
unique(Metabolites_Origin$Origin)
Metabolites_Origin$sample <- "SourceMF"
table(Metabolites_Origin$Origin)
ggplot(Metabolites_Origin) +
  geom_bar(aes(x =sample,y=value,fill = Origin),stat = "identity",position = "fill") +
  scale_fill_manual(values  = c("#4D7587","#85AEC9","#A0B3A1","#976276","#E1AB8D"))+theme_classic()+
  xlab("")+ylab("% of Origin")+theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "right",legend.title = element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+xlab("")+ylab("% of Origin")+
  scale_fill_manual(name = '',
                    values = c('Microbiota'='#4D7587','Co-Metabolism'='#85AEC9','Host'='#A0B3A1','Food related'='#976276','Others'='#E1AB8D'),
                    label = c('Microbiota'='Microbiota (num=14)','Co-Metabolism'='Co-Metabolism (num=27)','Host'='Host (num=2)',
                              'Food related'='Food related (num=143)','Others'='Others (num=194)'))



ggdata = list(c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism")]),c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Host","Co-Metabolism")]))
names(ggdata) <- c("Microbiota","Host")
ggdata
#BiocManager::install("ggvenn")
#library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 8,
         stroke_color = "black",
         fill_color = c("#85AEC9","#A0B3A1"),
         set_name_color = c("#85AEC9","#A0B3A1")) 
p

MPEA_Results_Co.Metabolism <- read.csv("score50/source_MetOrigin_20220805160748//03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
MPEA_Results_Mico.Metabolism <- read.csv("score50/source_MetOrigin_20220805160748//03_Function_Analysis/MPEA_Results_Microbiota.csv")
MPEA_Results_H.Metabolism <- read.csv("score50/source_MetOrigin_20220805160748//03_Function_Analysis/MPEA_Results_Host.csv")
MPEA_Results_Co.Metabolism$group <- "Co_Metabolism"
MPEA_Results_Mico.Metabolism$group <- "Micro_Metabolism"
MPEA_Results_H.Metabolism$group <- "H_Metabolism"

MPEA_Results.Metabolism <- rbind(MPEA_Results_Co.Metabolism,MPEA_Results_Mico.Metabolism)%>%rbind(.,MPEA_Results_H.Metabolism)
MPEA_Results.Metabolism$logPvalue <- -log(MPEA_Results.Metabolism$Pvalue)
data <- MPEA_Results.Metabolism[order(MPEA_Results.Metabolism$Pvalue,decreasing = T),]
data$group <- factor(data$group,levels = c("Micro_Metabolism","H_Metabolism","Co_Metabolism"))

ggplot(data, aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#03989E","#D0CEC1","#FF9800"))+
  ylab("")+xlab("log Pvalue")+
  #geom_text(aes(label=Sig..Metabolites),hjust=1,colour="white")+
  scale_y_discrete(limits=unique(data$Name))+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "top",legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))



#metab_origin

Metabolites_Origin <- read.csv("score50/Metabolites_Origin.csv")
View(Metabolites_Origin)
Metabolites_Origin$value <- 1
Metabolites_Origin <- Metabolites_Origin[Metabolites_Origin$Name%in%mix_select$Metabolite,]
Metabolites_Origin$Origin <- gsub("Drug related", "Others",Metabolites_Origin$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin$Origin <- factor(Metabolites_Origin$Origin,levels = c("Co-Metabolism","Microbiota","Host","Food related","Others"  ))
unique(Metabolites_Origin$Origin)
Metabolites_Origin$sample <- "SourceMF_metab"

ggplot(Metabolites_Origin) +
  geom_bar(aes(x =sample,y=value,fill = Origin),stat = "identity",position = "fill") +
  scale_fill_manual(values  = c("#4D7587","#85AEC9","#A0B3A1","#976276","#E1AB8D"))+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text = element_text(size = 20),
        legend.text = element_text(size=14),legend.title  = element_text(size=14),legend.position = "right")+
  xlab("")+ylab("% of Origin")

ggdata = list(c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism")]),c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Host","Co-Metabolism")]))
names(ggdata) <- c("Microbiota","Host")
ggdata
#BiocManager::install("ggvenn")
#library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 8,
         stroke_color = "white",
         fill_color = c("#85AEC9","#A0B3A1"),
         set_name_color = c("#85AEC9","#A0B3A1")) 
p

MPEA_Results_Co.Metabolism <- read.csv("score50/MetOrigin_20220502172449/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
MPEA_Results_Mico.Metabolism <- read.csv("score50/MetOrigin_20220502172449/03_Function_Analysis/MPEA_Results_Microbiota.csv")
MPEA_Results_H.Metabolism <- read.csv("score50/MetOrigin_20220502172449/03_Function_Analysis/MPEA_Results_Host.csv")
MPEA_Results.Metabolism <- rbind(MPEA_Results_Co.Metabolism,MPEA_Results_Mico.Metabolism)%>%rbind(.,MPEA_Results_H.Metabolism)
MPEA_Results.Metabolism$logPvalue <- -log(MPEA_Results.Metabolism$Pvalue)
ggplot(MPEA_Results.Metabolism[order(MPEA_Results.Metabolism$Pvalue,decreasing = T),], aes(x=logPvalue,y=Name)) +
  geom_bar(stat="identity",fill="#03989E") +theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text = element_text(size = 20),
        legend.text = element_text(size=12),legend.position = "right")+
  ylab("")+xlab("log Pvalue")+
  geom_text(aes(label=Sig..Metabolites),hjust=1.2,colour="white")+
  scale_y_discrete(limits=MPEA_Results_Co.Metabolism[order(MPEA_Results_Co.Metabolism$Pvalue,decreasing = T),]$Name)


