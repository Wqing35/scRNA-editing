####计算count与EI的相关系数并可视化
#读取数据
lncRNA_16 <- read.table("~/SPRINT/test_hippo_sets/GW16/featurecounts_out/regular_lncRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(lncRNA_16) <- c("exon_region","count","EI")

lncRNA_18 <- read.table("~/SPRINT/test_hippo_sets/GW18/featurecounts_out/regular_lncRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(lncRNA_18) <- c("exon_region","count","EI")

lncRNA_20 <- read.table("~/SPRINT/test_hippo_sets/GW20_1/featurecounts_out/regular_lncRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(lncRNA_20) <- c("exon_region","count","EI")

lncRNA_22 <- read.table("~/SPRINT/test_hippo_sets/GW22_2/featurecounts_out/regular_lncRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(lncRNA_22) <- c("exon_region","count","EI")

lncRNA_25 <- read.table("~/SPRINT/test_hippo_sets/GW25/featurecounts_out/regular_lncRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(lncRNA_25) <- c("exon_region","count","EI")
#mRNA
mRNA_16 <- read.table("~/SPRINT/test_hippo_sets/GW16/featurecounts_out/regular_mRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(mRNA_16) <- c("exon_region","count","EI")

mRNA_18 <- read.table("~/SPRINT/test_hippo_sets/GW18/featurecounts_out/regular_mRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(mRNA_18) <- c("exon_region","count","EI")

mRNA_20 <- read.table("~/SPRINT/test_hippo_sets/GW20_1/featurecounts_out/regular_mRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(mRNA_20) <- c("exon_region","count","EI")

mRNA_22 <- read.table("~/SPRINT/test_hippo_sets/GW22_2/featurecounts_out/regular_mRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(mRNA_22) <- c("exon_region","count","EI")

mRNA_25 <- read.table("~/SPRINT/test_hippo_sets/GW25/featurecounts_out/regular_mRNA_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(mRNA_25) <- c("exon_region","count","EI")
#all1k
all_16 <- read.table("~/SPRINT/test_hippo_sets/GW16/featurecounts_out/regular_all_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(all_16) <- c("exon_region","count","EI")

all_18 <- read.table("~/SPRINT/test_hippo_sets/GW18/featurecounts_out/regular_all_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(all_18) <- c("exon_region","count","EI")

all_20 <- read.table("~/SPRINT/test_hippo_sets/GW20_1/featurecounts_out/regular_all_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(all_20) <- c("exon_region","count","EI")

all_22 <- read.table("~/SPRINT/test_hippo_sets/GW22_2/featurecounts_out/regular_all_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(all_22) <- c("exon_region","count","EI")

all_25 <- read.table("~/SPRINT/test_hippo_sets/GW25/featurecounts_out/regular_all_EI_inExon_over1.txt",sep = '\t',header = T)
colnames(all_25) <- c("exon_region","count","EI")
# 
# #####要计算每一个exon的count与EI的相关性，用5个样本的同一个exon与对应的EI（5 X 5）？？？？
# #会舍弃掉部分未shared的exon片段
# tmp_lncRNA <- intersect(lncRNA_16$exon_region,lncRNA_18$exon_region)
# tmp_lncRNA <- intersect(tmp_lncRNA,lncRNA_20$exon_region)
# tmp_lncRNA <- intersect(tmp_lncRNA,lncRNA_22$exon_region)
# tmp_lncRNA <- intersect(tmp_lncRNA,lncRNA_25$exon_region)
# 
# all_cors_lncRNA <- c()
# for(i in 1:length(tmp_lncRNA)){
#   character_1 <- c(lncRNA_16[lncRNA_16$exon_region==tmp_lncRNA[i],3],lncRNA_18[lncRNA_18$exon_region==tmp_lncRNA[i],3],lncRNA_20[lncRNA_20$exon_region==tmp_lncRNA[i],3],lncRNA_22[lncRNA_22$exon_region==tmp_lncRNA[i],3],lncRNA_25[lncRNA_25$exon_region==tmp_lncRNA[i],3])
#   character_2 <- c(lncRNA_16[lncRNA_16$exon_region==tmp_lncRNA[i],2],lncRNA_18[lncRNA_18$exon_region==tmp_lncRNA[i],2],lncRNA_20[lncRNA_20$exon_region==tmp_lncRNA[i],2],lncRNA_22[lncRNA_22$exon_region==tmp_lncRNA[i],2],lncRNA_25[lncRNA_25$exon_region==tmp_lncRNA[i],2])
#   each_cor <- cor(character_1,character_2,method = 'pearson')
#   all_cors_lncRNA <- c(all_cors_lncRNA,each_cor)
# }
# 
#####计算exon region的并集以增加数据集大小
#mRNA
tmp_mRNA_union <- union(mRNA_16$exon_region,mRNA_18$exon_region)
tmp_mRNA_union <- union(tmp_mRNA_union,mRNA_20$exon_region)
tmp_mRNA_union <- union(tmp_mRNA_union,mRNA_22$exon_region)
tmp_mRNA_union <- union(tmp_mRNA_union,mRNA_25$exon_region)

#保留在至少两个数据集中找到的region
final_region_mRNA <- c()
for(i in 1:length(tmp_mRNA_union)){
  flag1 <- mRNA_16$exon_region==tmp_mRNA_union[i]
  flag2 <- mRNA_18$exon_region==tmp_mRNA_union[i]
  flag3 <- mRNA_20$exon_region==tmp_mRNA_union[i]
  flag4 <- mRNA_22$exon_region==tmp_mRNA_union[i]
  flag5 <- mRNA_25$exon_region==tmp_mRNA_union[i]
  
  flag <- c(flag1,flag2,flag3,flag4,flag5)
  if(length(grep("TRUE",flag,ignore.case = T)) > 3){
    final_region_mRNA <- c(final_region_mRNA,tmp_mRNA_union[i])
  }
}
all_cors_mRNA_union <- c()
for(i in 1:length(final_region_mRNA)){
  character_1 <- c(mRNA_16[mRNA_16$exon_region==final_region_mRNA[i],3],mRNA_18[mRNA_18$exon_region==final_region_mRNA[i],3],mRNA_20[mRNA_20$exon_region==final_region_mRNA[i],3],mRNA_22[mRNA_22$exon_region==final_region_mRNA[i],3],mRNA_25[mRNA_25$exon_region==final_region_mRNA[i],3])
  character_2 <- c(mRNA_16[mRNA_16$exon_region==final_region_mRNA[i],2],mRNA_18[mRNA_18$exon_region==final_region_mRNA[i],2],mRNA_20[mRNA_20$exon_region==final_region_mRNA[i],2],mRNA_22[mRNA_22$exon_region==final_region_mRNA[i],2],mRNA_25[mRNA_25$exon_region==final_region_mRNA[i],2])
  each_cor <- cor(character_1,character_2,method = 'pearson')
  all_cors_mRNA_union <- c(all_cors_mRNA_union,each_cor)
}
all_cors_mRNA_union[is.na(all_cors_mRNA_union)] <- 0.000001

#lncRNA
tmp_lncRNA_union <- union(lncRNA_16$exon_region,lncRNA_18$exon_region)
tmp_lncRNA_union <- union(tmp_lncRNA_union,lncRNA_20$exon_region)
tmp_lncRNA_union <- union(tmp_lncRNA_union,lncRNA_22$exon_region)
tmp_lncRNA_union <- union(tmp_lncRNA_union,lncRNA_25$exon_region)

#保留在至少两个数据集中找到的region
final_region_lncRNA <- c()
for(i in 1:length(tmp_lncRNA_union)){
  flag1 <- lncRNA_16$exon_region==tmp_lncRNA_union[i]
  flag2 <- lncRNA_18$exon_region==tmp_lncRNA_union[i]
  flag3 <- lncRNA_20$exon_region==tmp_lncRNA_union[i]
  flag4 <- lncRNA_22$exon_region==tmp_lncRNA_union[i]
  flag5 <- lncRNA_25$exon_region==tmp_lncRNA_union[i]
  
  flag <- c(flag1,flag2,flag3,flag4,flag5)
  if(length(grep("TRUE",flag,ignore.case = T)) > 3){
    final_region_lncRNA <- c(final_region_lncRNA,tmp_lncRNA_union[i])
  }
}
all_cors_lncRNA_union <- c()
for(i in 1:length(final_region_lncRNA)){
  character_1 <- c(lncRNA_16[lncRNA_16$exon_region==final_region_lncRNA[i],3],lncRNA_18[lncRNA_18$exon_region==final_region_lncRNA[i],3],lncRNA_20[lncRNA_20$exon_region==final_region_lncRNA[i],3],lncRNA_22[lncRNA_22$exon_region==final_region_lncRNA[i],3],lncRNA_25[lncRNA_25$exon_region==final_region_lncRNA[i],3])
  character_2 <- c(lncRNA_16[lncRNA_16$exon_region==final_region_lncRNA[i],2],lncRNA_18[lncRNA_18$exon_region==final_region_lncRNA[i],2],lncRNA_20[lncRNA_20$exon_region==final_region_lncRNA[i],2],lncRNA_22[lncRNA_22$exon_region==final_region_lncRNA[i],2],lncRNA_25[lncRNA_25$exon_region==final_region_lncRNA[i],2])
  each_cor <- cor(character_1,character_2,method = 'pearson')
  all_cors_lncRNA_union <- c(all_cors_lncRNA_union,each_cor)
}
all_cors_lncRNA_union[is.na(all_cors_lncRNA_union)] <- 0.000001

#all
tmp_all_union <- union(all_16$exon_region,all_18$exon_region)
tmp_all_union <- union(tmp_all_union,all_20$exon_region)
tmp_all_union <- union(tmp_all_union,all_22$exon_region)
tmp_all_union <- union(tmp_all_union,all_25$exon_region)

#保留在至少两个数据集中找到的region
final_region_all <- c()
for(i in 1:length(tmp_all_union)){
  flag1 <- all_16$exon_region==tmp_all_union[i]
  flag2 <- all_18$exon_region==tmp_all_union[i]
  flag3 <- all_20$exon_region==tmp_all_union[i]
  flag4 <- all_22$exon_region==tmp_all_union[i]
  flag5 <- all_25$exon_region==tmp_all_union[i]
  
  flag <- c(flag1,flag2,flag3,flag4,flag5)
  if(length(grep("TRUE",flag,ignore.case = T)) > 3){
    final_region_all <- c(final_region_all,tmp_all_union[i])
  }
}
all_cors_all_union <- c()
for(i in 1:length(final_region_all)){
  character_1 <- c(all_16[all_16$exon_region==final_region_all[i],3],all_18[all_18$exon_region==final_region_all[i],3],all_20[all_20$exon_region==final_region_all[i],3],all_22[all_22$exon_region==final_region_all[i],3],all_25[all_25$exon_region==final_region_all[i],3])
  character_2 <- c(all_16[all_16$exon_region==final_region_all[i],2],all_18[all_18$exon_region==final_region_all[i],2],all_20[all_20$exon_region==final_region_all[i],2],all_22[all_22$exon_region==final_region_all[i],2],all_25[all_25$exon_region==final_region_all[i],2])
  each_cor <- cor(character_1,character_2,method = 'pearson')
  all_cors_all_union <- c(all_cors_all_union,each_cor)
}

all_cors_all_union[is.na(all_cors_all_union)] <- 0.000001
###########random的结果包含的region种类多，平均cor值高于mRNA，此处不对random组做展示
#3个region的pair数量保持一致
# all_cors_random1k_union <- sample(all_cors_all_union,1k)
# all_cors_mRNA1k_union <- sample(all_cors_mRNA_union,1k)
# all_cors_lncRNA1k_union <- sample(all_cors_lncRNA_union,1k)

#取一次随机结果可信吗？？？
#取10次，做平均
mean_random1k <- c()
mean_mRNA1k <- c()
mean_lncRNA1k <- c()
for(time in c(1:10)){
  tmp_all_cors_random1k_union <- sample(all_cors_all_union,500)
  tmp_all_cors_mRNA1k_union <- sample(all_cors_mRNA_union,500)
  tmp_all_cors_lncRNA1k_union <- sample(all_cors_lncRNA_union,500)

  mean_random1k <- c(mean_random1k,mean(tmp_all_cors_random1k_union))
  mean_mRNA1k <- c(mean_mRNA1k,mean(tmp_all_cors_mRNA1k_union))
  mean_lncRNA1k <- c(mean_lncRNA1k,mean(tmp_all_cors_lncRNA1k_union))
}

data_all=data.frame(x=1:length(final_region_all), y=all_cors_all_union, z=rep('all',times=length(final_region_all)))
colnames(data_all) <- c('Num','Correlation_coeffecient','Group')

#mRNA+lncRNA仍然用随机取1000pair的方法
# mean_mRNA1k <- c()
# mean_lncRNA1k <- c()
# for(time in c(1:10)){
#   tmp_all_cors_mRNA1k_union <- sample(all_cors_mRNA_union,1000)
#   tmp_all_cors_lncRNA1k_union <- sample(all_cors_lncRNA_union,1000)
# 
#   mean_mRNA1k <- c(mean_mRNA1k,mean(tmp_all_cors_mRNA1k_union))
#   mean_lncRNA1k <- c(mean_lncRNA1k,mean(tmp_all_cors_lncRNA1k_union))
# }


###############合并mRNA+lncRNA的系数并分组
############生成以上count与EI相关系数的概率密度图
#生成all exon区域的count与EI关系
# ord <- sample(c(1:length($Num)),1000)
# all_set <- data_all_union[ord,]
# all_set$Group <- rep('all',times=1000)

library(ggplot2)
data_all_union=data.frame(x=1:(length(final_region_mRNA)+length(final_region_lncRNA)), y=c(all_cors_mRNA_union,all_cors_lncRNA_union),z=rep(c('mRNA RNA editing','lncRNA RNA editing'),times=c(length(final_region_mRNA),length(final_region_lncRNA))))
colnames(data_all_union) <- c('Num','Correlation_coeffecient','Group')

data_all_union_Wzall <- rbind(data_all_union,data_all)
colnames(data_all_union_Wzall) <- c('Num','Correlation_coeffecient','Group')

pdf("~/SPRINT/analysis/regular/figure/分组相关系数密度图_WizRandom.pdf")
ggplot(data=data_all_union_Wzall, aes(x=Correlation_coeffecient, group=Group,color=Group)) +
  scale_color_manual(
     values=c("lncRNA RNA editing"="orange", "mRNA RNA editing"="purple","all"="grey")) +
  #   labels=c("lncRNA","mRNA"))+
  geom_density(adjust=1.5, alpha=1,lwd=1,linetype = 1) +
  geom_vline(xintercept = mean(all_cors_mRNA_union),color = "purple",lwd = 1,linetype = 2) +
  geom_vline(xintercept = mean(all_cors_lncRNA_union),color = "orange",lwd = 1,linetype = 2) +
  geom_vline(xintercept = mean(all_cors_all_union),color = "grey",lwd = 1,linetype = 2) 
dev.off()

