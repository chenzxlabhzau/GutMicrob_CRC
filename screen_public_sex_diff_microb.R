#####################################
## 公共数据的下载和样本的筛选
#####################################

#### step1:下载数据 ######
library(curatedMetagenomicData)
??curatedMetagenomicData

## 查看有哪些文件
curatedMetagenomicData("YuJ_20.+")
curatedMetagenomicData("YachidaS_20.+")

count=list()
abundance=list()

study=c("FengQ_2015",
        "GuptaA_2019",
        "HanniganGD_2017",
        "ThomasAM_2018a",
        "ThomasAM_2018b",
        "ThomasAM_2019_c",
        "WirbelJ_2018",
        "YachidaS_2019",
        "ZellerG_2014",
        "YuJ_2015",
        "VogtmannE_2016"
)


for (i in study) {
  count[[i]]<-curatedMetagenomicData(paste0(i,".relative_abundance"), 
                                     dryrun = FALSE,
                                     ## download count data
                                     counts = TRUE, 
                                     ## use the long name of the species
                                     rownames = "long")
  
  count[[i]] <- assay(count[[i]][[1]]) |>
    as.matrix() |>
    as.data.frame()
  
  
  abundance[[i]]<-curatedMetagenomicData(paste0(i,".relative_abundance"), 
                                         dryrun = FALSE,
                                         ## Download relative abundance
                                         counts = F, 
                                         ## use the short name of the species
                                         rownames = "short")
  
  abundance[[i]] <- assay(abundance[[i]][[1]]) |>
    as.matrix() |>
    as.data.frame()
  
}


all_study_abundance<-abundance[-1]
all_study_count<-count[-1]



all_sample<-sampleMetadata %>%
  filter(study_name %in% c("FengQ_2015",
                           "GuptaA_2019",
                           "HanniganGD_2017",
                           "ThomasAM_2018a",
                           "ThomasAM_2018b",
                           "ThomasAM_2019_c",
                           "WirbelJ_2018",
                           "YachidaS_2019",
                           "ZellerG_2014",
                           "YuJ_2015",
                           "VogtmannE_2016"
  ))

save(all_study_abundance,all_study_count,all_sample,file = "/home/zhangy/project/meta_analyse/data_11study/data_11study.RData")


#############################step1 sample_info####################################
# [1] "X"                   "study_name"          "sample_id"           "subject_id"          "study_condition"    
# [6] "disease"             "age"                 "age_category"        "gender"              "BMI"                
# [11] "country"             "PMID"                "number_reads"        "number_bases"        "minimum_read_length"
# [16] "median_read_length"  "tnm"                 "disease_location"    "curator" 
load("/home/zhangy/project/meta_analyse/data_11study/data_11study.RData")
colnames(all_sample)
all_sample$gender[is.na(all_sample$study_condition)]
all_sample$study_condition[is.na(all_sample$study_condition)]<-"control"
all_sample<-all_sample%>%
  filter(!(study_condition=="carcinoma_surgery_history"))%>%
  filter(!(is.na(gender)))
#write.csv(all_sample,"/home/panxl/Metagenomic_data/sample_info/all_sampleinfo_metadata.csv")

##################################step2 筛选数据 #####################
################################2.1除去CRC手术史样本和没有性别或者样本47#######
###提取文件
load("/home/zhangy/project/meta_analyse/data_11study/data_11study.RData")
df<-names(all_study_count)
for(i in df){
  data<-all_study_count[[i]]%>%as.data.frame()
  write.table(data,paste0("/home/panxl/Metagenomic_data/Count_data/",i,"_merged_count_table.txt",spe=""),sep="\t",quote=F)
}
###合并文件
filename<-list.files(path ="/home/panxl/Metagenomic_data/Count_data/",pattern = "table.txt",full.names = TRUE )
for (i in 1:length(filename)) {
  assign(as.character(filename[[i]]),read.csv(filename[i],sep="\t",header = TRUE))
}
filename_species<-list.files(path ="/home/panxl/Metagenomic_data/Count_data/",pattern = "table.txt",all.files  = TRUE )
filename_species<-filename_species[-1]
tmp<-read.csv("/home/panxl/Metagenomic_data/Count_data/FengQ_2015_merged_count_table.txt",sep="\t")
tmp$clade_name<-rownames(tmp)
for(i in filename_species){
  data<-read.csv(paste("/home/panxl/Metagenomic_data/Count_data/",i,sep=""),sep="\t")
  data$clade_name<-rownames(data)
  tmp<-full_join(tmp,data,by="clade_name")
}
data<-tmp%>%as.data.frame()
rownames(data)<-data$clade_name
data<-data[,-155]
for (i in 1:ncol(data)) {
  data[,i][is.na(data[,i])]<-0
}
all_sample<-read.csv("/home/panxl/Metagenomic_data/sample_info/all_sampleinfo_metadata.csv")
which(is.na(all_sample$age))

all_sample=all_sample[-419,]
table(all_sample$study_condition,all_sample$gender)
#          female male
# adenoma     77  132
# control    311  383
# CRC        257  444
data1<-data
colnames(data1)<-gsub("[.]","-",colnames(data1))
all_samp1<-all_sample
all_samp1$sample_id<-gsub("[.]","-",all_samp1$sample_id)
all_sample_metadata<-all_samp1
all_count<-data1
all_count<-all_count[,colnames(all_count)%in%all_sample_metadata$sample_id]
#save(all_sample_metadata,all_count,file ="/home/panxl/Metagenomic_data/result/all_sample_count.RData")

#####加重庆数据
load("/home/panxl/Metagenomic_data/result/all_sample_count.RData")
cq_count<-read.csv("/home/panxl/CRC_Metagenomic_Metab/chongqing_data/merged_abundance_counts_species.txt",sep="\t")
colnames(cq_count)<-gsub("_modify2_counts","",colnames(cq_count))
all_count$clade_name<-rownames(all_count)

data<-full_join(all_count,cq_count,by="clade_name")
rownames(data)<-data$clade_name
load("/home/panxl/CRC_Metagenomic_Metab/chongqing_data/metaphlan_species_sample_info.RData")
colnames(sample_info)
colnames(all_sample_metadata)
all_sample_metadata<-all_sample_metadata[,-c(1:2)]
sample_info$study_name<-"YangJ_2020"
sample_info$sample_id<-sample_info$Run
sample_info$subject_id<-sample_info$Run
sample_info$study_condition<-ifelse(sample_info$State=="Case","CRC","control")
#sample_info$disease<-"CRC"
sample_info$disease<-ifelse(sample_info$State=="Case","CRC","healthy")
sample_info$age<-sample_info$Age
sample_info$age_category<-ifelse(sample_info$Age<=65,"adult","senior")
sample_info$gender<-ifelse(sample_info$Gender=="Female","female","male")
sample_info$country<-"CHN"
sample_info$PMID<-31971861
cq_count_number_reads<-as.data.frame(apply(cq_count[,-1],2,sum))
cq_count_number_reads$number_reads<-rownames(cq_count_number_reads)
colnames(cq_count_number_reads)<-c("number_reads","sample_id")
sample_info$number_reads<-cq_count_number_reads[match(sample_info$sample_id,cq_count_number_reads$sample_id),]$number_reads
sample_info$number_bases<-"NA"
sample_info$minimum_read_length<-"NA"
sample_info$median_read_length<-"NA"
sample_info$curator<-"NA"
sample_info$BMI<-sample_info$`BMI(kg/m²)`
sample_info$tnm<-sample_info$`TNM stage`
sample_info$disease_location<-"NA"
sample_info1<-sample_info[,c(21:38)]
colnames(sample_info1)

all_sample_metadata<-all_sample[,colnames(sample_info1)]  
all_sample_metadata1<-rbind(all_sample_metadata,sample_info1)
data1<-data[,colnames(data)%in%all_sample_metadata1$sample_id]
all_sample_metadata<-all_sample_metadata1
all_count<-data1
dim(all_sample_metadata)
dim(all_count)
save(all_sample_metadata,all_count,file ="/home/panxl/Metagenomic_data/result/all_sample_count.RData")

######################2.2除去测序深度小于100万的样本 4个样本##############
load("/home/panxl/Metagenomic_data/result/all_sample_count.RData")
which(is.na(all_sample_metadata$age))#853
#all_sample_metadata=all_sample_metadata[-853,]
all_sample_metadata1<-all_sample_metadata[all_sample_metadata$number_reads>=1000000,]
all_count<-all_count[,colnames(all_count)%in%all_sample_metadata1$sample_id]
######################2.3 除去细菌含量大于总量50%和小于1/n * 1/100 的样本  2个样本###########
#all_sample_sum<-tapply(all_sample_metadata1$number_reads, all_sample_metadata1$study_name, sum) %>% as.data.frame()
df<-unique(all_sample_metadata1$study_name)
study_condition<-unique(all_sample_metadata1$study_condition)
all_sample_metadata2<-purrr::map(study_condition,function(y){
  tmp<-data_frame()
  for(x in df){
    data<-all_sample_metadata1[all_sample_metadata1$study_name==x&all_sample_metadata1$study_condition==y,]
    sum<-sum(data$number_reads)
    data1<-data[data$number_reads<=0.5*sum&data$number_reads>=1/nrow(data)*sum*0.1,]
    tmp<-rbind(tmp,data1)
  } 
  return(tmp)
})%>%dplyr::bind_rows()
all_count<-all_count[,colnames(all_count)%in%all_sample_metadata2$sample_id]
all_sample_metadata<-all_sample_metadata2
save(all_sample_metadata,all_count,file ="/home/panxl/Metagenomic_data/result/all_sample_count.RData")

###################################2.4筛选菌  ######################################
load("/home/panxl/Metagenomic_data/result/all_sample_count.RData")
for(i in 1:ncol(all_count)){
  all_count[,i][is.na(all_count[,i])]=0
  
  
}
# convert to relative abundance
all_merged_model_metaphlan_species1=as.matrix(all_count)
all_rel_metaphlan_species = prop.table(all_merged_model_metaphlan_species1, 2)#将feat.ab.crc按列求百分比
lib.size = 5000
all_rar_metaphlan_species= t(rrarefy(floor(t(all_merged_model_metaphlan_species1)), sample=lib.size))#绘制稀释曲线
dim(all_rar_metaphlan_species)

temp.max.ab = t(sapply(row.names(all_rel_metaphlan_species),
                       FUN=function(marker){sapply(unique(all_sample_metadata$study_name),
                                                   FUN=function(study_name, marker){
                                                     max.ab = max(all_rel_metaphlan_species[marker, which(all_sample_metadata$study_name==study_name)])
                                                   },
                                                   marker=marker)}))
dim(temp.max.ab)
### low abundance filter
ab.cutoff=0.001
f.idx = rowSums(temp.max.ab >= ab.cutoff) >= 6 &
  row.names(all_rel_metaphlan_species) != '-1'

all_ab.red =all_count[f.idx,]
all_rel.red = all_rel_metaphlan_species[f.idx,]
all_rar.red = all_rar_metaphlan_species[f.idx,]#绘制稀释曲线
cat('Retaining', sum(f.idx), 'features after low-abundance filtering...\n')#251 
all_sample_metadata$Label=case_when(all_sample_metadata$study_condition=="control"&all_sample_metadata$gender=="female"~"H-F",
                                    all_sample_metadata$study_condition=="control"&all_sample_metadata$gender=="male"~"H-M",
                                    all_sample_metadata$study_condition=="CRC"&all_sample_metadata$gender=="female"~"C-F",
                                    all_sample_metadata$study_condition=="CRC"&all_sample_metadata$gender=="male"~"C-M",
                                    all_sample_metadata$study_condition=="adenoma"&all_sample_metadata$gender=="female"~"A-F",
                                    all_sample_metadata$study_condition=="adenoma"&all_sample_metadata$gender=="male"~"A-M",
)
dim(all_ab.red)
save(all_ab.red,all_rel.red,all_rar.red,all_sample_metadata,all_count,file ="/home/panxl/Metagenomic_data/result/all_data.RData")



#############################
## 公共数据样本数量的展示图
#############################
load("~/project/meta_analyse/data_11study/data_11study.RData")
load("/home/panxl/Metagenomic_data/result/all_data.RData")
sample_meta<-all_sample_metadata[which(all_sample_metadata$disease %in% c("CRC","healthy","adenoma")),]

ttt<-all_sample[all_sample$sample_id %in% sample_meta$sample_id,]
ttt<-ttt[,c("sample_id","alcohol_numeric","smoker")]

meta<-left_join(sample_meta,ttt)
table(meta$disease)

# Get Data
meta=meta[which(meta$disease!="adenoma"),]
meta <- meta[!(meta$country=="CAN" & meta$study_name=="HanniganGD_2017"),]

#####样本数量统计#######
t<-as.data.frame(table(meta$study_name)) %>% arrange(Freq)
meta$study_name<-factor(meta$study_name,levels = t$Var1)

df1<-meta[which(meta$study_condition=="CRC"),]
df1<-as.data.frame(table(df1$study_name,df1$gender))
colnames(df1)<-c("study_name","gender","number")
df1$gender<-ifelse(df1$gender=="male","M-C","F-C")

df2<-meta[which(meta$study_condition=="control"),]
df2<-as.data.frame(table(df2$study_name,df2$gender))
colnames(df2)<-c("study_name","gender","number")
df2$gender<-ifelse(df2$gender=="male","M-H","F-H")

df2$number<-df2$number*(-1)

COLOR=c(
  "F-H"="#FFD9C0",
  # "F-A"="#CC9C75",
  "F-C"="#A25B5B",
  "M-H"="#8AB6D6",
  # "M-A"="#84A7E6",
  "M-C"="#23689B"
)


ggplot() +
  geom_col(data = df2,aes(x = study_name,y=number,fill=gender),position = 'stack',width = 0.7) +
  geom_col(data=df1,aes(x=study_name,y=number,fill=gender),position = 'stack',width = 0.7)+
  scale_fill_manual(values = COLOR)+
  theme_classic()+
  coord_flip()+
  geom_hline(yintercept=0)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,hjust = 1,angle = 45))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size=14),
        legend.title  = element_text(size=14),
        #legend.position = "top",
        strip.text = element_text(size=14))


#####################################
## Health中影响H-F和H-M比较的混杂因素分析
#####################################

load("~/project/meta_analyse/data_11study/data_11study.RData")
load("/home/panxl/Metagenomic_data/result/all_data.RData")
sample_meta<-all_sample_metadata[which(all_sample_metadata$disease %in% c("CRC","healthy","adenoma")),]

ttt<-all_sample[all_sample$sample_id %in% sample_meta$sample_id,]
ttt<-ttt[,c("sample_id","alcohol_numeric","smoker")]

meta<-left_join(sample_meta,ttt)
meta=meta[which(meta$disease!="adenoma"),]
meta <- meta[!(meta$country=="CAN" & meta$study_name=="HanniganGD_2017"),]
## 有二型糖尿病、脂肪肝的都不要了

count<-all_ab.red[,meta$sample_id]
abundce<-all_rel.red[,meta$sample_id]



##################H-FM#####################
library(coin)
# load("/home/panxl/CRC_Metagenomic_Metab/meta_count.RData")
# meta=all_sample_info

meta_TFM=meta[meta$Label%in%c("H-F","H-M"),]
studies <- meta %>% pull(study_name) %>% unique
meta_TFM$block<-meta_TFM$study_name
colnames(meta_TFM)[2]<-"sample_id"

all_rel.red_TFM=abundce[,meta_TFM$sample_id]

# all_rel.red_TFM=all_rel.red[,colnames(all_rel.red)%in%meta_TFM$sample_id]

p.val <- matrix(NA, nrow=nrow(all_rel.red_TFM), ncol=length(studies)+1, 
                dimnames=list(row.names(all_rel.red_TFM), c(studies, 'all')))
fc <- p.val
aucs.mat <- p.val
aucs.all  <- vector('list', nrow(all_rel.red_TFM))

cat("Calculating effect size for every feature...\n")
pb <- txtProgressBar(max=nrow(all_rel.red_TFM), style=3)
colnames(meta_TFM)[2]="sample_id"
# caluclate wilcoxon test and effect size for each feature and study
for (f in row.names(all_rel.red_TFM)) {
  
  # for each study
  for (s in studies) {
    # s="FengQ_2015" f="Faecalibacterium_prausnitzii"
    
    x <- all_rel.red_TFM[f, meta_TFM %>% filter(study_name==s) %>% 
                           filter(Label=="H-F") %>% pull(sample_id)]
    
    
    y <- all_rel.red_TFM[f, meta_TFM %>% filter(study_name==s) %>% 
                           filter(Label=="H-M") %>% pull(sample_id)]
    # Wilcoxon
    p.val[f,s] <- wilcox.test(x, y, exact=FALSE)$p.value
    # FC
    q.p <- quantile(log10(x+1e-05), probs=seq(.1, .9, .05))
    q.n <- quantile(log10(y+1e-05), probs=seq(.1, .9, .05))
    fc[f, s] <- sum(q.p - q.n)/length(q.p)
  }
  
  # calculate effect size for all studies combined
  # Wilcoxon + blocking factor
  d <- data.frame(y=all_rel.red_TFM[f,], 
                  x=meta_TFM$Label, block=meta_TFM$block)
  #d$x <-as.factor(d$x)
  d$x <-factor(d$x) 
  d$block <-factor(d$block) 
  b<-wilcox_test(y ~ x|block, data=d)
  p.val[f,'all'] <- pvalue(wilcox_test(y ~ x | block, data=d)$p.value)
  
  x1 <- all_rel.red_TFM[f, meta_TFM %>% 
                          filter(Label=="H-F") %>% pull(sample_id)]
  y1 <- all_rel.red_TFM[f, meta_TFM %>% 
                          filter(Label=="H-M") %>% pull(sample_id)]
  #p.val[f,'all'] <- wilcox.test(x1, y1, exact=FALSE)$p.value
  
  # other metrics   
  # x <- all_rel.red_TFM[f, meta_TFM %>%  filter(Label=="C-F") %>% pull(sample_id)]
  # y <- all_rel.red_TFM[f, meta_TFM %>% filter(Label=="H-F") %>% pull(sample_id)]
  # FC
  fc[f, 'all'] <- mean(fc[f, studies])
  # AUC
  # aucs.mat[f,'all'] <- c(roc(controls=y, cases=x, 
  #                            direction='<', ci=TRUE, auc=TRUE)$ci)[2]
  
  # progressbar
  setTxtProgressBar(pb, (pb$getVal()+1))
}


# multiple hypothesis correction
p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method="fdr"),
                    check.names = FALSE)
T_FM_p.adj=p.adj

#save(T_FM_p.adj,TH_p.adj,file = "/home/panxl/CRC_Metagenomic_Metab/meta_padj.RData")



##################################
meta_TFM$Library_Size<-meta_TFM$number_reads


meta_1 <- meta_TFM %>%
  # age
  mutate(age_factor=as.factor(
    cut(meta_TFM$age, breaks = quantile(meta_TFM$age), labels=c(1,2,3,4)))) %>%
  # bmi
  mutate(bmi_factor=as.factor(
    cut(meta_TFM$BMI, breaks = c(0, 25, 30, 100),
        labels=c('lean', 'overweight', 'obese')))) %>%
  # library size
  mutate(lib_size_factor=as.factor(
    cut(meta_TFM$Library_Size, breaks = quantile(meta_TFM$Library_Size),
        labels=c(1,2,3,4))))

meta_1$alcohol_factor<-ifelse(meta_1$alcohol_numeric>0,1,0)
meta_1$alcohol_factor<-factor(meta_1$alcohol_factor)
unique(meta_1$country)
meta_1$continent_factor<-if_else(meta_1$country %in% c("CHN","IND","JPN"),
                                 "Asian",ifelse(meta_1$country %in% c("USA","CAN"),
                                                "North American","European"))
meta_1$continent_factor<-factor(meta_1$continent_factor)

unique(meta_1$disease_location)



meta_1$smoker<-ifelse(meta_1$smoker=="yes",1,
                      ifelse(meta_1$smoker=="no",0,
                             meta_1$smoker))
meta_1$smoker<-factor(meta_1$smoker)

####
#  variance explained by disease status
ss.disease <- apply(all_rel.red_TFM, 1, FUN=function(x, Label){
  rank.x <- rank(x)/length(x)#百分比排名
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)#所有样本方差
  ss.o.i <- sum(vapply(unique(Label), function(l){
    sum((rank.x[Label==l] - mean(rank.x[Label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)#每个分组样本方差
  return(1-ss.o.i/ss.tot)
}, Label=meta_1 %>% pull(Label))   ## 这个地方不用study_condition，但是现在是性别

# calculate trimmed mean abundance
t.mean <- apply(all_rel.red_TFM, 1, mean, trim=0.1)

df.plot.all <- tibble(
  species=rownames(all_rel.red_TFM),
  disease=ss.disease,
  t.mean=t.mean,
  adj.p.val=p.adj[rownames(all_rel.red_TFM), 'all'],
  meta.significance=p.val[rownames(all_rel.red_TFM), 'all'] < 0.05)


df.list <- list()
for (meta.var in c('age_factor',"bmi_factor","study_name", 
                   "continent_factor","lib_size_factor",
                   #"disease_location",
                   "smoker",
                   "alcohol_factor")){
  cat('###############################\n', meta.var, '\n')
  meta.c <- meta_1 %>%
    filter(!is.na(eval(parse(text=meta.var))))
  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$gender, meta.c %>% pull(meta.var)))
  print(table(meta.c$study_name))
  feat.red <- all_rel.red_TFM[,meta.c$sample_id]
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, Label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(Label), function(l){
      sum((rank.x[Label==l] - mean(rank.x[Label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
  }, Label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}


#########################

data=df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  mutate(facet=case_when( #type=='Gender' ~ 'Sex',
    type=='age_factor' ~ 'Age',
    type=='study_name' ~ 'Study',
    type=='bmi_factor' ~ 'BMI',
    # type=='study_condition' ~ 'study_condition',
    type=='continent_factor' ~ 'continent',
    type=='lib_size_factor' ~ 'library size',
    type=='alcohol_factor' ~ 'alcohol',
    #type=='disease_location' ~ 'location',
    TRUE ~ type))
colnames(data)  
alpha.meta <- 0.05

load("~/project/meta_analyse/output/HFM_result.Rdata")
alpha.meta <- 0.05
data<-TFM_data

p2<-ggplot(data[which(data$facet!="continent"),], aes(x=disease, y=meta, size=t.mean+1e-08, 
                                                      col=meta.significance)) +
  geom_point(shape=19) +
  
  xlab('Variance by gender status(Health)') +
  ylab('Variance by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=2) +
  theme(
    axis.text = element_text(size=24, color="black"),
    #text = element_text(size=15),
    axis.text.x=element_text(size=24,face = "bold"),
    axis.title = element_text(size=24,face = "bold"),
    axis.text.y = element_text(size =24,face = "bold"),
    legend.text=element_text(size=24,face = "bold"),
    legend.title=element_text(size=24,face = "bold"),
    strip.text = element_text(size=24,face="bold",color="black"),    #控制分面标题的文字
    strip.background = element_blank(),###分页标题背景
    panel.grid.minor = element_blank()##网格线
  )+
  
  #scale_x_continuous(breaks = seq(from=0, to=0.32, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
  scale_colour_manual(
    # values = alpha(c('black', '#CC071E'),
    #                                  alpha=c(0.4, .75)  ),
    values = c("#0000001A","#CC071EBF"),
    name=paste0('Significance\n(', alpha.meta, ' Pvalue)')) +
  scale_size_area(name='Trimmed mean\nabundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')
p2

HFM_df.plot.all=df.plot.all
HFM_data<-data
table(HFM_data$meta.significance)   ## 50个微生物
save(HFM_data,HFM_df.plot.all,file = "../meta_analyse/output/HFM_result.Rdata")


ggplot(data[data$facet=="Study",], aes(x=disease, 
                                       y=meta, 
                                       size=t.mean+1e-08, 
                                       col=meta.significance)) +
  geom_point(shape=19) +
  
  xlab('Variance by gender status') +
  ylab('Variance by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=2) +
  theme(
    axis.text = element_text(size=15, color="black"),
    #text = element_text(size=15),
    axis.text.x=element_text(hjust=1, size=18,face = "bold"),
    axis.title = element_text(size=24,face = "bold"),
    axis.text.y = element_text(size =24,face = "bold"),
    legend.text=element_text(size=24,face = "bold"),
    legend.title=element_text(size=24,face = "bold"),
    strip.text = element_text(size=24,face="bold",color="black"),    #控制分面标题的文字
    strip.background = element_blank(),###分页标题背景
    panel.grid.minor = element_blank()##网格线
  )+
  #scale_x_continuous(breaks = seq(from=0, to=0.32, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
  scale_colour_manual(values = alpha(c('black', '#CC071E'),
                                     alpha=c(0.4, .75)
  ),
  name=paste0('Significance\n(', alpha.meta, ' Pvalue)')) +
  scale_size_area(name='Trimmed mean\nabundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')

p2<-p2+
  theme(legend.position = 'none')

library(patchwork)
p2|p1 +plot_layout(guides = 'collect') 



#####################################
## CRC中影响C-F和C-M比较的混杂因素分析
#####################################

load("~/project/meta_analyse/data_11study/data_11study.RData")
load("/home/panxl/Metagenomic_data/result/all_data.RData")
sample_meta<-all_sample_metadata[which(all_sample_metadata$disease %in% c("CRC","healthy","adenoma")),]

ttt<-all_sample[all_sample$sample_id %in% sample_meta$sample_id,]
ttt<-ttt[,c("sample_id","alcohol_numeric","smoker")]

meta<-left_join(sample_meta,ttt)
## 有二型糖尿病、脂肪肝的都不要了

count<-all_ab.red[,meta$sample_id]
abundce<-all_rel.red[,meta$sample_id]



##################T-FM#####################
library(coin)
# load("/home/panxl/CRC_Metagenomic_Metab/meta_count.RData")
# meta=all_sample_info

meta_TFM=meta[meta$Label%in%c("C-F","C-M"),]
studies <- meta %>% pull(study_name) %>% unique
meta_TFM$block<-meta_TFM$study_name
colnames(meta_TFM)[2]<-"sample_id"

all_rel.red_TFM=abundce[,meta_TFM$sample_id]

# all_rel.red_TFM=all_rel.red[,colnames(all_rel.red)%in%meta_TFM$sample_id]

p.val <- matrix(NA, nrow=nrow(all_rel.red_TFM), ncol=length(studies)+1, 
                dimnames=list(row.names(all_rel.red_TFM), c(studies, 'all')))
fc <- p.val
aucs.mat <- p.val
aucs.all  <- vector('list', nrow(all_rel.red_TFM))

cat("Calculating effect size for every feature...\n")
pb <- txtProgressBar(max=nrow(all_rel.red_TFM), style=3)
colnames(meta_TFM)[2]="sample_id"
# caluclate wilcoxon test and effect size for each feature and study
for (f in row.names(all_rel.red_TFM)) {
  
  # for each study
  for (s in studies) {
    # s="FengQ_2015" f="Faecalibacterium_prausnitzii"
    
    x <- all_rel.red_TFM[f, meta_TFM %>% filter(study_name==s) %>% 
                           filter(Label=="C-F") %>% pull(sample_id)]
    
    
    y <- all_rel.red_TFM[f, meta_TFM %>% filter(study_name==s) %>% 
                           filter(Label=="C-M") %>% pull(sample_id)]
    # Wilcoxon
    p.val[f,s] <- wilcox.test(x, y, exact=FALSE)$p.value
    # # AUC
    # aucs.all[[f]][[s]] <- c(roc(controls=y, cases=x, 
    #                             direction='<', ci=TRUE, auc=TRUE)$ci)
    # aucs.mat[f,s] <- c(roc(controls=y, cases=x, 
    #                        direction='<', ci=TRUE, auc=TRUE)$ci)[2]
    # FC
    q.p <- quantile(log10(x+1e-05), probs=seq(.1, .9, .05))
    q.n <- quantile(log10(y+1e-05), probs=seq(.1, .9, .05))
    fc[f, s] <- sum(q.p - q.n)/length(q.p)
  }
  
  # calculate effect size for all studies combined
  # Wilcoxon + blocking factor
  d <- data.frame(y=all_rel.red_TFM[f,], 
                  x=meta_TFM$Label, block=meta_TFM$block)
  #d$x <-as.factor(d$x)
  d$x <-factor(d$x) 
  d$block <-factor(d$block) 
  b<- 
    p.val[f,'all'] <- pvalue(wilcox_test(y ~ x | block, data=d))
  
  x1 <- all_rel.red_TFM[f, meta_TFM %>% 
                          filter(Label=="C-F") %>% pull(sample_id)]
  y1 <- all_rel.red_TFM[f, meta_TFM %>% 
                          filter(Label=="C-M") %>% pull(sample_id)]
  #p.val[f,'all'] <- wilcox.test(x1, y1, exact=FALSE)$p.value
  
  # other metrics   
  # x <- all_rel.red_TFM[f, meta_TFM %>%  filter(Label=="C-F") %>% pull(sample_id)]
  # y <- all_rel.red_TFM[f, meta_TFM %>% filter(Label=="H-F") %>% pull(sample_id)]
  # FC
  fc[f, 'all'] <- mean(fc[f, studies])
  # AUC
  # aucs.mat[f,'all'] <- c(roc(controls=y, cases=x, 
  #                            direction='<', ci=TRUE, auc=TRUE)$ci)[2]
  
  # progressbar
  setTxtProgressBar(pb, (pb$getVal()+1))
}


# multiple hypothesis correction
p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method="fdr"),
                    check.names = FALSE)
T_FM_p.adj=p.adj

#save(T_FM_p.adj,TH_p.adj,file = "/home/panxl/CRC_Metagenomic_Metab/meta_padj.RData")



##################################
meta_TFM$Library_Size<-meta_TFM$number_reads


meta_1 <- meta_TFM %>%
  # age
  mutate(age_factor=as.factor(
    cut(meta_TFM$age, breaks = quantile(meta_TFM$age), labels=c(1,2,3,4)))) %>%
  # bmi
  mutate(bmi_factor=as.factor(
    cut(meta_TFM$BMI, breaks = c(0, 25, 30, 100),
        labels=c('lean', 'overweight', 'obese')))) %>%
  # library size
  mutate(lib_size_factor=as.factor(
    cut(meta_TFM$Library_Size, breaks = quantile(meta_TFM$Library_Size),
        labels=c(1,2,3,4))))

meta_1$alcohol_factor<-ifelse(meta_1$alcohol_numeric>0,1,0)
meta_1$alcohol_factor<-factor(meta_1$alcohol_factor)
unique(meta_1$country)
meta_1$continent_factor<-if_else(meta_1$country %in% c("CHN","IND","JPN"),
                                 "Asian",ifelse(meta_1$country %in% c("USA","CAN"),
                                                "North American","European"))
meta_1$continent_factor<-factor(meta_1$continent_factor)

unique(meta_1$disease_location)
meta_1$disease_location<-if_else(meta_1$disease_location=="lc",
                                 "left_colon",
                                 meta_1$disease_location)
meta_1$disease_location<-if_else(meta_1$disease_location=="rc",
                                 "right_colon",
                                 meta_1$disease_location)
meta_1$disease_location<-if_else(meta_1$disease_location=="left_colon_and_right_colon",
                                 "colon",
                                 meta_1$disease_location)
meta_1$disease_location<-if_else(meta_1$disease_location %in% c("rectum_and_left_colon",
                                                                "rectum,_left_colon_and_right_colon",
                                                                "left_colon_and_rectum"),
                                 "rectum-colon",
                                 meta_1$disease_location)

meta_1[which(meta_1$disease_location=="NA"),]$disease_location<-NA
meta_1$disease_location<-factor(meta_1$disease_location)

meta_1$smoker<-ifelse(meta_1$smoker=="yes",1,
                      ifelse(meta_1$smoker=="no",0,
                             meta_1$smoker))
meta_1$smoker<-factor(meta_1$smoker)

####
#  variance explained by disease status
ss.disease <- apply(all_rel.red_TFM, 1, FUN=function(x, Label){
  rank.x <- rank(x)/length(x)#百分比排名
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)#所有样本方差
  ss.o.i <- sum(vapply(unique(Label), function(l){
    sum((rank.x[Label==l] - mean(rank.x[Label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)#每个分组样本方差
  return(1-ss.o.i/ss.tot)
}, Label=meta_1 %>% pull(Label))   ## 这个地方不用study_condition，但是现在是性别

# calculate trimmed mean abundance
t.mean <- apply(all_rel.red_TFM, 1, mean, trim=0.1)

df.plot.all <- tibble(
  species=rownames(all_rel.red_TFM),
  disease=ss.disease,
  t.mean=t.mean,
  adj.p.val=p.adj[rownames(all_rel.red_TFM), 'all'],
  meta.significance=p.val[rownames(all_rel.red_TFM), 'all'] < 0.05)


df.list <- list()
for (meta.var in c('age_factor',"bmi_factor","study_name", 
                   "continent_factor","lib_size_factor",
                   "disease_location",
                   "smoker",
                   "alcohol_factor")){
  cat('###############################\n', meta.var, '\n')
  meta.c <- meta_1 %>%
    filter(!is.na(eval(parse(text=meta.var))))
  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$gender, meta.c %>% pull(meta.var)))
  print(table(meta.c$study_name))
  feat.red <- all_rel.red_TFM[,meta.c$sample_id]
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, Label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(Label), function(l){
      sum((rank.x[Label==l] - mean(rank.x[Label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
  }, Label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}

CFM_df.plot.all=df.plot.all
CFM_data<-data
table(CFM_data$meta.significance)   
save(CFM_data,CFM_df.plot.all,file = "../meta_analyse/output/CFM_result.Rdata")

#########################
data=df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  mutate(facet=case_when( #type=='Gender' ~ 'Sex',
    type=='age_factor' ~ 'Age',
    type=='study_name' ~ 'Study',
    type=='bmi_factor' ~ 'BMI',
    # type=='study_condition' ~ 'study_condition',
    type=='continent_factor' ~ 'continent',
    type=='lib_size_factor' ~ 'library size',
    type=='alcohol_factor' ~ 'alcohol',
    type=='disease_location' ~ 'location',
    TRUE ~ type))



load("~/project/meta_analyse/output/CFM_result.Rdata")
colnames(data)  
alpha.meta <- 0.05

data<-CFM_data

p1<-ggplot(data[which(data$facet!="continent"),], 
           aes(x=disease, y=meta, size=t.mean+1e-08, 
               col=meta.significance)) +
  geom_point(shape=19) +
  
  xlab('Variance by gender status(Tumor)') +
  ylab('Variance by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=2) +
  theme(axis.text = element_text(size=18, color="black"),
        #text = element_text(size=15),
        axis.text.x=element_text(size=18,face = "bold"),
        axis.title = element_text(size=18,face = "bold"),
        axis.text.y = element_text(size =18,face = "bold"),
        legend.text=element_text(size=24,face = "bold"),
        legend.title=element_text(size=24,face = "bold"),
        strip.text = element_text(size=24,face="bold",color="black"),    #控制分面标题的文字
        strip.background = element_blank(),###分页标题背景
        panel.grid.minor = element_blank()##网格线
  )+
  
  #scale_x_continuous(breaks = seq(from=0, to=0.32, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
  scale_colour_manual(
    # values = alpha(c('black', '#CC071E'),
    #                alpha=c(0.1, .75)),
    values = c("#0000001A","#CC071EBF"),
    name=paste0('Significance\n(', alpha.meta, ' FDR)')) +
  scale_size_area(name='Trimmed mean\nabundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')

p1

data<-CFM_data
alpha.meta=0.05

ggplot(data[data$facet=="Study",], aes(x=disease, 
                                       y=meta, 
                                       size=t.mean+1e-08, 
                                       col=meta.significance)) +
  geom_point(shape=19) +
  
  xlab('Variance by gender status(Tumor)') +
  ylab('Variance by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=2) +
  theme(
    axis.text = element_text(size=15, color="black"),
    #text = element_text(size=15),
    axis.text.x=element_text(size=18,face = "bold"),
    axis.title = element_text(size=18,face = "bold"),
    axis.text.y = element_text(size =18,face = "bold"),
    legend.text=element_text(size=18,face = "bold"),
    legend.title=element_text(size=18,face = "bold"),
    strip.text = element_text(size=18,face="bold",color="black"),    #控制分面标题的文字
    strip.background = element_blank(),###分页标题背景
    panel.grid.minor = element_blank()##网格线
  )+
  #scale_x_continuous(breaks = seq(from=0, to=0.32, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
  scale_colour_manual(values = alpha(c('black', '#CC071E'),
                                     alpha=c(0.4, .75)
  ),
  name=paste0('Significance\n(', alpha.meta, ' Pvalue)')) +
  scale_size_area(name='Trimmed mean\nabundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')



#####################
## alpha多样性分析
#####################

# # ##############################################################################
## 多样性的计算
load("~/project/meta_analyse/data_11study/data_11study.RData")
load("/home/panxl/Metagenomic_data/result/all_data.RData")
sample_meta<-all_sample_metadata[which(all_sample_metadata$disease %in% c("CRC","healthy","adenoma")),]

ttt<-all_sample[all_sample$sample_id %in% sample_meta$sample_id,]
ttt<-ttt[,c("sample_id","alcohol_numeric","smoker")]

meta<-left_join(sample_meta,ttt)
table(meta$disease)

# Get Data
meta=meta[which(meta$disease!="adenoma"),]
meta <- meta[!(meta$country=="CAN" & meta$study_name=="HanniganGD_2017"),]
# count<-all_ab.red[,meta$sample_id]
# abundce<-all_rel.red[,meta$sample_id]

feat.all<-all_rar.red

feat.all <- prop.table(feat.all)

meta <- meta[!(meta$country=="CAN" & meta$study_name=="HanniganGD_2017"),]
feat.all<-feat.all[,meta$sample_id]
feat.all_count=all_rar.red[,meta$sample_id]
dim(feat.all)
dim(meta)


# ##############################################################################
# Compute Alpha diversity

df.div <- meta %>% 
  ### ### OTU
  mutate(`All species` =estimateR(t(feat.all_count))[1, ]) %>%   #observed_species
  
  select(Label, study_name, `All species`) %>%
  gather(key=type, value=diversity, -c(study_name, Label)) 

df.div$study_name<-factor(df.div$study_name)
df.div$Label<-factor(df.div$Label)




###########整体P值的计算###############
## 下面的可以作为考虑study因素之后的组间的P值，可以在合并的时候，标上去
# blocked wilcoxon
library(coin)
df.div %>% 
  filter(type=='All species') %>% 
  filter(Label %in%c("C-M","C-F")) %>%
  wilcox_test(diversity~Label|study_name, data=.)

# oberserved: Z = -2.0449, p-value = 0.04087  ## 显著
# data:  diversity by
# Label (C-F, C-M) 
# stratified by study_name
# Z = -2.0449, p-value = 0.04087

df.div[which(df.div$Label %in% c("H-M","H-F")),] %>% 
  filter(type=='All species') %>% 
  wilcox_test(diversity~Label|study_name, data=.)

# oberserved:
# data:  diversity by
# Label (H-F, H-M) 
# stratified by study_name
# Z = 0.86047, p-value = 0.3895
# alternative hypothesis: true mu is not equal to 0


df.div[which(df.div$Label %in% c("H-M","C-M")),] %>% 
  filter(type=='All species') %>% 
  wilcox_test(diversity~Label|study_name, data=.)

# data:  diversity by
# Label (C-M, H-M) 
# stratified by study_name
# Z = 4.7068, p-value = 2.516e-06
# alternative hypothesis: true mu is not equal to 0

df.div[which(df.div$Label %in% c("H-F","C-F")),] %>% 
  filter(type=='All species') %>% 
  wilcox_test(diversity~Label|study_name, data=.)

# data:  diversity by
# Label (C-F, H-F) 
# stratified by study_name
# Z = 1.838, p-value = 0.06606
# alternative hypothesis: true mu is not equal to 0


# anova
summary(aov(rank~Label+study_name+, 
            data=df.div %>% 
              filter(type=='All species') %>% 
              mutate(rank=rank(diversity))))
# oberserved:
# Df    Sum Sq Mean Sq F value   Pr(>F)    
#   Label          3   2774392  924797   7.669 4.38e-05 ***
#   study_name    11  67878235 6170749  51.173  < 2e-16 ***
#   Residuals   1409 169905032  120586                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#######################
name2label<-c("H-F"="F-H",
              "H-M"="M-H",
              #"A-F"="F-A",
              #"A-M"="M-A",
              "C-F"="F-C",
              "C-M"="M-C")
df.div$Label<-unname(name2label[df.div$Label])


df.div$country<-meta[which(meta$study_name%in%df.div$study_name),]$country
df.div$study_name<-paste0(df.div$study_name,"(",df.div$country,")")
df.div$Label <- factor(df.div$Label,levels = c("F-H","F-C","M-H","M-C"))

COLOR=c(
  "F-H"="#FFD9C0",
  # "F-A"="#CC9C75",
  "F-C"="#A25B5B",
  "M-H"="#8AB6D6",
  # "M-A"="#84A7E6",
  "M-C"="#23689B"
)

df.div[which(df.div$study_name!="HanniganGD_2017(CAN)"),]%>% 
  ggplot(aes(x=Label, fill=Label, y=diversity)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  #geom_boxplot(position=position_dodge(0.9)) +
  geom_boxplot(lwd=1.2,fatten = 1.2)+
  ylab('Alpha diversity (Observed species)') +
  xlab('') +
  scale_fill_manual(values=COLOR) + 
  geom_signif(comparisons = list(
    c("F-H", "F-C"),
    c("M-H", "M-C"),
    c("F-H", "M-H"),
    c("F-C", "M-C")
  ),
  test = "wilcox.test",
  map_signif_level = T,step_increase = 0.05)+
  facet_wrap(~study_name, scales = 'free_y',nrow =3,
             strip.position = "top") + 
  theme_bw(base_size = 16,base_rect_size=2) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  guides(fill=FALSE)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,hjust = 1,angle = 45))+
  theme(axis.text.y = element_text(size = 20))+
  theme(
    legend.text = element_text(size=24),
    legend.title  = element_text(size=24),
    strip.text = element_text(size=12))



##########
# 整体绘图

re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y
}

#删除含有outliers(NA)的行
df2<-df.div%>%
  group_by(Label)%>%
  mutate(value = re_outlier(diversity))


df2<-df2%>%
  group_by(Label)

df2$Label<-factor(df2$Label,levels = c("F-H", 
                                       "F-C",
                                       "M-H",
                                       "M-C"))


na.omit(df2)%>% 
  ggplot(aes(x=Label, fill=Label, y=diversity)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2)+
  ylab('Alpha diversity (Observed species)') +
  xlab('') +
  scale_fill_manual(values=COLOR) + 
  # facet_wrap(~study_name, scales = 'free_y',nrow =2,
  #            strip.position = "top") + 
  theme_bw(base_size = 16,base_rect_size=2) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  guides(fill=FALSE)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,hjust = 1,angle = 45))+
  theme(axis.text.y = element_text(size = 20))+
  theme(
    legend.text = element_text(size=24),
    legend.title  = element_text(size=24),
    strip.text = element_text(size=18))+
  # stat_compare_means( label ="p.format",
  #                     method = "wilcox.test",
  #                     size=5)+
  geom_signif(comparisons = list(
    # c("C-F", "C-M"),
    # c("H-M", "C-M"),
    # c("H-F", "C-F"),
    # c("H-F", "H-M")
    c("F-H", "F-C"),
    c("M-H", "M-C"),
    c("F-H", "M-H"),
    c("F-C", "M-C")
  ),
  test = "wilcox.test",
  map_signif_level = F,step_increase = 0.05)

#####################
## beta多样性分析
#####################
library("labdsv")  # pco
library("coin")
library("vegan")
library("yaml")
library("ggpubr")
library("cowplot")
library("tidyverse")

###########beta多样性的图###############
load("~/project/meta_analyse/data_11study/data_11study.RData")
load("/home/panxl/Metagenomic_data/result/all_data.RData")
sample_meta<-all_sample_metadata[which(all_sample_metadata$disease %in% c("CRC","healthy","adenoma")),]

ttt<-all_sample[all_sample$sample_id %in% sample_meta$sample_id,]
ttt<-ttt[,c("sample_id","alcohol_numeric","smoker")]

meta<-left_join(sample_meta,ttt)
table(meta$disease)

# Get Data
meta=meta[which(meta$disease!="adenoma"),]
meta <- meta[!(meta$country=="CAN" & meta$study_name=="HanniganGD_2017"),]
# abundce=all_rel.red[,meta$sample_id]
# count=all_ab.red[,meta$sample_id]

t<-as.data.frame(table(meta$study_name)) %>% arrange(Freq)
meta$study_name<-factor(meta$study_name,levels = t$Var1)

# meta<-all_sample_metadata
feat.all=all_ab.red[,meta$sample_id]    ## 1424

# compute PCoA
library(phyloseq)
dist = vegdist(t(feat.all), method = 'bray')
pco.results = pco(dist, k=2)

axis.1.title <- paste('PCo1 [', 
                      round((pco.results$eig[1]/sum(pco.results$eig))*100,1),
                      '%]', sep='')
axis.2.title <- paste('PCo2 [', 
                      round((pco.results$eig[2]/sum(pco.results$eig))*100,1),
                      '%]', sep='')

df.plot <- tibble(Axis1 = -1*pco.results$points[,1],
                  Axis2 = pco.results$points[,2],
                  Sample_ID = rownames(pco.results$points),
                  Group=meta$Label,
                  Study=meta$study_name)

df.plot$Study<-factor(df.plot$Study,levels = t$Var1)
# ##############################################################################
# subplots
table(meta$study_name)

study_color<-c(
  'FengQ_2015'="#8F62C3",
  'GuptaA_2019'="#135AAF",   ## 深蓝色
  "HanniganGD_2017"="#F5C31E",
  'ThomasAM_2018a'="#F73E00",
  "ThomasAM_2018b"="#0FCBC4",   ##青色
  'ThomasAM_2019_c'="#1F7553",
  "VogtmannE_2016"="#E56A54",
  "WirbelJ_2018"="#F73E00",
  "YachidaS_2019"="#D4EEA7",
  'YangJ_2020'="#57D65E",
  'YuJ_2015'="#9F1326",
  "ZellerG_2014"="#F8AC6A"
)

barplot(rep(1,times=12),col=study_color)

group_value=c(
  "F-H"="#FFD9C0",
  # "F-A"="#CC9C75",
  "F-C"="#A25B5B",
  "M-H"="#8AB6D6",
  # "M-A"="#84A7E6",
  "M-C"="#23689B"
)

name2label<-c("H-F"="F-H",
              "H-M"="M-H",
              #"A-F"="F-A",
              #"A-M"="M-A",
              "C-F"="F-C",
              "C-M"="M-C")
df.plot$Group<-unname(name2label[df.plot$Group])



guide = guide_legend(
  #direction = "horizontal",
  title.position = "right",   ## 标题放在左边
  title.hjust = 0.5,  ## 1是向上对齐，0是向下对齐，0.5居中
  # label.vjust = 1,
  title.theme = element_text(angle = 90)
)




# main plot
g.main <- df.plot %>% 
  ggplot(aes(x=Axis1, y=Axis2, shape=Group,col=Study)) +
  geom_point(size=2) + 
  scale_colour_manual(values=study_color) + 
  scale_shape_manual(values=c(19,1,15,0)) + 
  scale_x_continuous(position='top') +
  xlab(axis.1.title) + ylab(axis.2.title) +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        axis.ticks=element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())+
  theme(legend.position = c(0.92,0.4),
        legend.background = element_rect(fill = "transparent",
                                         colour = NA),
        legend.title = element_blank())

# study boxplot axis 1
g.s.1 <- df.plot %>% 
  mutate(Study=factor(Study, levels=names(study_color))) %>% 
  ggplot(aes(y=Axis1, x=Study, fill=Study)) + 
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.8)+
  geom_boxplot(lwd=0.8,fatten = 0.8)+
  scale_fill_manual(values=study_color, guide=FALSE) +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), axis.title.x = element_blank(),
        panel.grid = element_blank()) + 
  coord_flip()

# study boxplot axis 2
g.s.2 <- df.plot %>% 
  mutate(Study=factor(Study, levels=names(study_color))) %>% 
  ggplot(aes(y=Axis2, x=Study, fill=Study)) + 
  # geom_boxplot() + 
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.8)+
  geom_boxplot(lwd=0.8,fatten = 0.8)+
  scale_fill_manual(values=study_color, guide=FALSE) +
  scale_x_discrete(position='top') +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

# group plot axis1
g.g.1 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis1, fill=Group)) +
  # geom_boxplot() +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.8)+
  geom_boxplot(lwd=0.8,fatten = 0.8)+
  scale_fill_manual(values=group_value, guide=FALSE) + 
  ylab(axis.1.title) + 
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank()) + 
  coord_flip()
# group plot axis2
g.g.2 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis2, fill=Group)) +
  #geom_boxplot() +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.8)+
  geom_boxplot(lwd=0.8,fatten = 0.8)+
  scale_fill_manual(values=group_value, guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  ylab(axis.2.title) + 
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())

library(cowplot)
plot_grid(g.main, g.s.2, g.g.2, g.s.1, NULL, NULL, g.g.1, NULL, NULL,
          nrow=3,
          rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))


## 显著性就看下面的这些
# wilcox_test(df.plot[which(df.plot$Group%in%c("F-C","M-C")),]$Axis1~as.factor(df.plot[which(df.plot$Group%in%c("F-C","M-C")),]$Group) )
# wilcox_test(df.plot$Axis1~as.factor(df.plot$Group) )
kruskal.test(df.plot$Axis1~as.factor(df.plot$Group))
# Kruskal-Wallis chi-squared = 12.705, df = 3, p-value = 0.005319
# Asymptotic Wilcoxon-Mann-Whitney Test
# 
# data:  df.plot$Axis1 by as.factor(df.plot$Group) (older, young)
# Z = 17.358, p-value < 2.2e-16
# alternative hypothesis: true mu is not equal to 0

#wilcox_test(df.plot[which(df.plot$Group%in%c("F-C","M-C")),]$Axis2~as.factor(df.plot[which(df.plot$Group%in%c("F-C","M-C")),]$Group) )
kruskal.test(df.plot$Axis2~as.factor(df.plot$Group))
# Kruskal-Wallis chi-squared = 15.422, df = 3, p-value = 0.00149
# p-value = 0.7353
#wilcox_test(df.plot$Axis2~as.factor(df.plot$Group)|as.factor(df.plot$Study))
# Z = -0.32335, p-value = 0.7464

# kruskal.test(df.plot$Axis1~as.factor(df.plot$Group))

kruskal.test(df.plot$Axis1~as.factor(df.plot$Study))
#Kruskal-Wallis chi-squared = 199.36, df = 11, p-value < 2.2e-16
# Kruskal-Wallis rank sum test
# 
# data:  df.plot$Axis1 by as.factor(df.plot$Study)
# Kruskal-Wallis chi-squared = 106.75, df = 6, p-value < 2.2e-16

kruskal.test(df.plot$Axis2~as.factor(df.plot$Study))
#Kruskal-Wallis chi-squared = 552.94, df = 11, p-value < 2.2e-16
# Kruskal-Wallis chi-squared = 17.576, df = 6, p-value = 0.007384



wilcox_test(df.plot[which(df.plot$Group %in% c("F-C","M-C")),]$Axis1 ~ as.factor(df.plot[which(df.plot$Group %in% c("F-C","M-C")),]$Group)|as.factor(df.plot[which(df.plot$Group %in% c("F-C","M-C")),]$Study))
wilcox_test(df.plot[which(df.plot$Group %in% c("F-C","M-C")),]$Axis2 ~ as.factor(df.plot[which(df.plot$Group %in% c("F-C","M-C")),]$Group)|as.factor(df.plot[which(df.plot$Group %in% c("F-C","M-C")),]$Study))



#####################
## 差异分析和热图
#####################

#BiocManager::install("Maaslin2")
library(Maaslin2)

# ##############################################################################
load("~/project/meta_analyse/data_11study/data_11study.RData")
load("/home/panxl/Metagenomic_data/result/all_data.RData")
sample_meta<-all_sample_metadata[which(all_sample_metadata$disease %in% c("CRC","healthy","adenoma")),]

ttt<-all_sample[all_sample$sample_id %in% sample_meta$sample_id,]
ttt<-ttt[,c("sample_id","alcohol_numeric","smoker")]

meta<-left_join(sample_meta,ttt)
table(meta$disease)

# Get Data
meta=meta[which(meta$disease!="adenoma"),]
meta <- meta[!(meta$country=="CAN" & meta$study_name=="HanniganGD_2017"),]
abundce=all_rel.red[,meta$sample_id]
count=all_ab.red[,meta$sample_id]
table(meta$study_condition,meta$gender)

write.table(abundce,file = "../meta_analyse/output/supplementary/abundce.txt",
            sep = "\t",
            quote = T,row.names = T)

write.table(meta,file = "../meta_analyse/output/supplementary/meta.txt",
            sep = "\t",
            quote = T,row.names = T)

###############################################################################


#####C-M  C-F#########
df_input_metadata=meta[which(meta$Label %in% c("C-M","C-F")),]
df_input_metadata$Label<-factor(df_input_metadata$Label)
df_input_metadata<-df_input_metadata[order(df_input_metadata$Label),]
df_input_metadata$study_name<-factor(df_input_metadata$study_name)
df_input_metadata<-cbind(sample_id=df_input_metadata$sample_id,
                         df_input_metadata[,c(1,3:length(df_input_metadata))])

# df_input_data<-count[,df_input_metadata$sample_id]
df_input_data<-abundce[,df_input_metadata$sample_id]

df_input_data<-t(df_input_data)
df_input_data<-data.frame(sample_id=rownames(df_input_data),
                          df_input_data)

## 重写重读一遍
write_tsv(df_input_data,file = "../meta_analyse/output/maasilin/df_input_data.tsv")
df_input_data=read.csv("../meta_analyse/output/maasilin/df_input_data.tsv",sep = "\t")

write_tsv(df_input_metadata,file = "../meta_analyse/output/maasilin/df_input_metadata.tsv")
df_input_metadata=read.csv("../meta_analyse/output/maasilin/df_input_metadata.tsv",sep = "\t")

df_input_metadata$Label<-factor(df_input_metadata$Label)
df_input_metadata$study_name<-factor(df_input_metadata$study_name)

fit_data =Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "../meta_analyse/output/9-29/",
  fixed_effects = c('Label') , # include fixed effects as needed
  random_effects = c("study_name") # include random effects as needed
  #reference = 'Label, H-F'
  ,
  min_prevalence = 0.1,
  min_abundance = 0.0001
)

significant_results=read.csv("/home/zhangy/project/meta_analyse/output/9-29/significant_results.tsv",sep = "\t")
C_mvsf<-significant_results
write_tsv(C_mvsf,file = "../meta_analyse/output/9-29/C_mvsf.tsv")
# identical(C_mvsf,significant_results)






#####H-M  H-F#########
df_input_metadata=meta[which(meta$Label %in% c("H-M","H-F")),]
df_input_metadata$Label<-factor(df_input_metadata$Label)
df_input_metadata<-df_input_metadata[order(df_input_metadata$Label),]
df_input_metadata$study_name<-factor(df_input_metadata$study_name)
df_input_metadata<-cbind(sample_id=df_input_metadata$sample_id,
                         df_input_metadata[,c(1,3:length(df_input_metadata))])

df_input_data<-abundce[,df_input_metadata$sample_id]
df_input_data<-t(df_input_data)
df_input_data<-data.frame(sample_id=rownames(df_input_data),
                          df_input_data)

## 重写重读一遍
write_tsv(df_input_data,file = "../meta_analyse/output/maasilin/df_input_data.tsv")
df_input_data=read.csv("../meta_analyse/output/maasilin/df_input_data.tsv",sep = "\t")

write_tsv(df_input_metadata,file = "../meta_analyse/output/maasilin/df_input_metadata.tsv")
df_input_metadata=read.csv("../meta_analyse/output/maasilin/df_input_metadata.tsv",sep = "\t")

df_input_metadata$Label<-factor(df_input_metadata$Label)
df_input_metadata$study_name<-factor(df_input_metadata$study_name)

fit_data =Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  output = "../meta_analyse/output/9-29-Health/",  
  fixed_effects = c('Label') , # include fixed effects as needed 
  random_effects = c("study_name") # include random effects as needed
  #reference = 'Label, H-F'
  ,
  min_prevalence = 0.1,
  min_abundance = 0.0001    
)

significant_results=read.csv("/home/zhangy/project/meta_analyse/output/9-29-Health/significant_results.tsv",sep = "\t")
H_mvsf<-significant_results

write_tsv(H_mvsf,file = "../meta_analyse/output/9-29-Health/H_mvsf.tsv")


######看交集与差集##########

H_mvsf<-significant_results
C_mvsf<-significant_results

C_mvsf=read.csv("../meta_analyse/output/maasilin_res/C_mvsf.tsv",sep = "\t")
H_mvsf=read.csv("../meta_analyse/output/maasilin_res/H_mvsf.tsv",sep = "\t")

C_mvsf_diff<-C_mvsf[which(C_mvsf$pval<0.05),]
H_mvsf_diff<-H_mvsf[which(H_mvsf$pval<0.05),]

a<-setdiff(H_mvsf_diff$feature,C_mvsf_diff$feature)
b<-setdiff(C_mvsf_diff$feature,H_mvsf_diff$feature)

specific<-C_mvsf_diff[which(C_mvsf_diff$feature%in%b),]
specific$value<-ifelse(specific$coef>0,"C-M","C-F")
View(specific[which(specific$feature%in%b),])



####韦恩图
library(ggvenn)
ggvenn(list("Diff micro(Health)"=H_mvsf_diff$feature,
            "Diff micro(CRC)"=C_mvsf_diff$feature),
       show_percentage = F,
       fill_alpha = 1,
       text_size = 8,
       stroke_color = "white",
       fill_color = c("#97C4B8","#E5CB9F"))

## 热图
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#336699", "white", "#CC3333"))
# mat<-count[b,meta[which(meta$Label %in% c("C-F","C-M")),]$sample_id]
mat<-abundce[b,meta[which(meta$Label %in% c("C-F","C-M")),]$sample_id]

library(ComplexHeatmap)
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#23689B","#A25B5B")),
                       labels =c("C-M","C-F"),      
                       labels_gp = gpar(col = "white", fontsize = 18)))

row_ha = rowAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#23689B","#A25B5B")), 
                       labels =c("C-M","C-F"), 
                       labels_gp = gpar(col = "white", fontsize = 18))
)

rownames(meta)<-meta$sample_id
Group = factor(meta[colnames(mat),]$Label,
               levels = c("C-M","C-F"))



row_Group = factor(specific$value,
                   levels =c("C-M","C-F"))

n <- t(scale(t(log10(mat+0.0001))))
round(sum(n[1,]),2)
range(n)

Heatmap(n,name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        left_annotation = row_ha,
        column_split = Group,
        row_split = row_Group,
        show_heatmap_legend = T,
        border = "black",
        show_column_names = F,
        column_names_gp = gpar(fontsize=14),
        show_row_dend = FALSE,
        show_row_names = T,row_names_gp = gpar(fontsize=14),
        column_title = NULL,cluster_rows = T,cluster_columns = F,
        row_names_max_width = max_text_width(
          rownames(n), 
          gp = gpar(fontsize = 14)
        ))
