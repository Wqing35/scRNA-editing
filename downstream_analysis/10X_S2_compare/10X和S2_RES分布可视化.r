library(trackViewer) #主要画图包
library(RColorBrewer) #引入颜色包
res_s2 <-  read.table("/disk1/wenqing/tmp_data/10X_smartseq2/result/s2/all_140_ESC/a.txt")$V1# 指定SNP 坐标位置
res_10X <- read.table("/disk1/wenqing/tmp_data/10X_smartseq2/result/10X/all_140_ESC/a.txt")$V1

all_res <- c(res_s2,res_10X)

all_res_wzSample <- as.data.frame(cbind(all_res,c(rep('s2',times=24),rep('10X',times=57))))

intron_widths <- sample(c(100,150,130,200,230,170,190,250),15,replace=T)


res_clu1 <- sort(sample(5:(intron_widths[1]-5),15,replace=F))
res_clu2 <- sort(sample((sum(intron_widths[1:2])+5):(sum(intron_widths[1:3]-5)),22,replace=F))
res_clu4 <- sort(sample((sum(intron_widths[1:4])+5):(sum(intron_widths[1:5]-5)),21,replace=F))
res_clu6 <- sort(sample((sum(intron_widths[1:6])+5):(sum(intron_widths[1:7]-5)),5,replace=F))
res_clu7 <- sort(sample((sum(intron_widths[1:8])+5):(sum(intron_widths[1:9]-5)),7,replace=F))
res_clu8 <- sort(sample((sum(intron_widths[1:10])+5):(sum(intron_widths[1:11]-5)),2,replace=F))
res_clu9 <- sort(sample((sum(intron_widths[1:12])+5):(sum(intron_widths[1:13]-5)),4,replace=F))
res_clu17 <- sort(sample((sum(intron_widths[1:14])+5):(sum(intron_widths[1:15]-5)),5,replace=F))

res_info <- all_res_wzSample[order(as.numeric(all_res_wzSample[,1],decreasing=F)),]$V2

new_all_res <- c(res_clu1,res_clu2,res_clu4,res_clu6,res_clu7,res_clu8,res_clu9,res_clu17)


sample.gr <- GRanges("chr16", IRanges(new_all_res, width=1, names=paste0('res',all_res))) 

intron_start <- c(1,intron_widths[1]+1,sum(intron_widths[1:2])+1,sum(intron_widths[1:3])+1,sum(intron_widths[1:4])+1,sum(intron_widths[1:5])+1,
                    sum(intron_widths[1:6])+1,sum(intron_widths[1:7])+1,sum(intron_widths[1:8])+1,sum(intron_widths[1:9])+1,sum(intron_widths[1:10])+1,sum(intron_widths[1:11])+1,
                    sum(intron_widths[1:12])+1,sum(intron_widths[1:13])+1,sum(intron_widths[1:14])+1)



#RBFOX1_intron_bed <- read.table("/disk1/wenqing/tmp_data/RBFOX1_intron.bed")

#sub_RBFOX1_intron_bed <- RBFOX1_intron_bed[c(1,2,4,6,7,8,9,17),]


features <- GRanges("chr16", IRanges(intron_start, # 设置block起使位置
                                    width=intron_widths, # 设置block 的长度
                                    names=c("intron1","exon2","intron2","exon3-exon4","intron4","exon5-exon6","intron6",
                                    "exon7","intron7","exon8","intron8","exon9","intron9","exon10-exon17","intron17"))) # 设置名字
features$fill <- c(brewer.pal(8,"Set2"),brewer.pal(7,"Set3")) #块的颜色

sample.gr$color <- c("red","red","blue","blue",rep("red",4),rep("blue",9),rep("red",15),"blue","red",rep("blue",3),
                        rep("red",28),rep("blue",5),rep("red",4),rep("blue",7)) #棒子上面的球的颜色
sample.gr$border <- sample(c("grey60", "grey50"), length(all_res), replace=TRUE) #棒子的颜色
sample.gr$alpha <- sample(100:200, length(all_res), replace = TRUE)/200   #设置透明度0-1之间,sample是生成100-200之间的随机数

sample.gr$label <- as.character(1:length(sample.gr)) #球内的字符
sample.gr$label.col <- "black" #球内的标签的颜色

# features$height <- c(0.02, 0.05, 0.04) #块的高度
sample.gr$score <- sample.int(1, length(sample.gr), replace = TRUE) #设置球的数量
pdf("~/res_distri_in2plats.pdf")
lolliplot(sample.gr, features,yaxis = F,xaxis = F,ylab=F, cex=0.7,family='Times New Roman')#yaxis设置不显示y轴
dev.off()


#BiocManager::install("trackViewer")
setwd("e:/Rcode/trackviewer")
library(trackViewer) #主要画图包
library(RColorBrewer) #引入颜色包
SNP <- c(10,12,23,250,300,303,400,500,510,535,650,680,700,1005,1200,1350,1402) # 指定SNP 坐标位置
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP))) # 设置棒棒图位置
features <- GRanges("chr1", IRanges(c(1, 401, 1001), # 设置block起使位置
                                    width=c(400, 600, 505), # 设置block 的长度
                                    names=c("UTR5","CDS","UTR3"))) # 设置名字
features$fill <- rep("red",3) #块的颜色
sample.gr$color <- c(rep("red",7),rep("blue", length(SNP)-7)) #棒子上面的球的颜色
sample.gr$border <- sample(c("grey60", "grey50"), length(SNP), replace=TRUE) #棒子的颜色
sample.gr$alpha <- sample(100:200, length(SNP), replace = TRUE)/200   #设置透明度0-1之间,sample是生成100-200之间的随机数

sample.gr$label <- as.character(1:length(sample.gr)) #球内的字符
sample.gr$label.col <- "black" #球内的标签的颜色

# features$height <- c(0.02, 0.05, 0.04) #块的高度
sample.gr$score <- sample.int(1, length(sample.gr), replace = TRUE) #设置球的数量
pdf("~/test.pdf")
lolliplot(sample.gr, features,yaxis = F,family='Times New Roman') #yaxis设置不显示y轴
dev.off()















intron_start <- c()
for(i in 1:nrow(sub_RBFOX1_intron_bed)){
    intron_start <- c(intron_start,c(sub_RBFOX1_intron_bed$V2[i],sub_RBFOX1_intron_bed$V3[i]))
}

intron_widths <- c()
for(i in (1:length(intron_start)-1)){
    intron_width <- intron_start[i+1]-intron_start[i]
    intron_widths <- c(intron_widths,intron_width)
}
