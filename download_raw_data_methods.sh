####下载SRR原始数据的一些尝试


#1: prefetch批量下载太缓慢，下载为sra格式后，还需使用fast-dump转换为诶gz格式、再解压
#不推荐！！！！1

#2.尝试使用Aspera批量下载SRR
#配置完成Aspera后，下载速度达到418Mb/s！！！！
#step1:安装Aspera
conda install -c hcc aspera-cli
#step2:获得密钥文件
which ascp
#/disk1/wenqing/anaconda3/envs/wq_py2/bin/ascp
#将bin及之后的代码替换为etc/asperaweb_id_dsa.openssh
#step3:ascp下载命令：
ascp  -vQT -l 500m -P33001 -k 1 -i \
/disk1/wenqing/anaconda3/envs/wq_py2/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR636/004/SRR6367154/SRR6367154_1.fastq.gz  ./
######### 主要使用参数
#-v 详细模式
#-Q 用于自适应流量控制，磁盘限制所需
#-T 设置为无需加密传输
#-l 最大下载速度，一般设为500m
#-P TCP 端口，一般为33001
#-k 断点续传，通常设为 1
#-i 免密下载的密钥文件

#step4:批量下载
#进入EBI数据库使用GSE号搜索后获得fastq的tsv文件
#批量下载代码为：
file='/disk1/wenqing/tmp_data/PFC_s2/all_fastq/undownload_read1.txt'

cat $file |while read id 
do
        ascp -vQT -l 500m -P33001 -k 1 -i \
        /disk1/wenqing/anaconda3/envs/wq_py2/etc/asperaweb_id_dsa.openssh \
        era-fasp@$id  .
done

#3.下载fastq.gz文件后，快速解压gz文件的方式
unpigz -d SRRXXX_*.fastq.gz

#4. 根据barcode将celltype fq提取出来

#5. 所有样本的celltype fq分别输入sprint_main_hisat2.py

#####注意：每个样本的大小在17-112G之间，上述所有处理步骤实现并行处理的方式为：
#####将所有涉及到下载、解压、提取celltype fq、sprint的操作，输入样本都分为两份
#####并尽可能提高每个步骤使用的线程数

