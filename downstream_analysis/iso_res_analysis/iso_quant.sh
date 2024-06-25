##AS analysis

mkdir -p SUPPA2/ref
cd SUPPA2/ref

###########step1：下载所需基因组文件
nohup wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz &
nohup wget  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz &
gunzip gencode.v37.annotation.gtf.gz
gunzip  gencode.v37.transcripts.fa.gz


###########step2：安装所需软件
conda create -n rna
conda activate rna  
conda install -y trim-galore
# python-3.7.9         | 57.3 MB,  openjdk-10.0.2       | 189.2 MB
conda install -y star hisat2 bowtie2 # 没有其它依赖
conda install -y subread  bedtools deeptools
conda install -y salmon=1.4.0
conda install  -y -c bioconda suppa

###########step3: salmon生成索引
nohup salmon index -t  gencode.v37.transcripts.fa  -i   gencode.v37.transcripts.salmon.index &

###########step4: 
# 使用SUPPA generateEvents 命令根据基因组的gtf注释文件生成所有的可变剪切事件，格式保存为ioe格式
# 针对  SE SS MX RI FL 
suppa.py generateEvents \
-i   gencode.v37.annotation.gtf \
-o gencode.v37.events -e SE SS MX RI FL -f ioe
#合并所有的ioe文件
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.ioe > gencode.v37.all.events.ioe

###########step5: fastqc
sample=(GW08 GW12 GW16_1_3 GW16_1_4 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${sample[@]};
do
    dir=/disk1/wenqing/tmp_data/PFC_s2/data/"$one"
    mkdir /disk1/wenqing/tmp_data/PFC_s2/qc/"$one"/
    fastqc  "$dir"/GABAergic_neurons_read*.fastq -o  /disk1/wenqing/tmp_data/PFC_s2/qc/"$one"/
done

###########step6: 剪切reads
sample=(GW08 GW12 GW16_1_3 GW16_1_4 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${sample[@]};
do
    dir=/disk1/wenqing/tmp_data/PFC_s2/data/"$one"/
    mkdir /disk1/wenqing/tmp_data/PFC_s2/clean_data/"$one"/
    trim_galore -j 10 -q 30 --phred33 --length 36 --stringency 3  -o /disk1/wenqing/tmp_data/PFC_s2/clean_data/"$one"/  "$dir"/GABAergic_neurons_read*.fastq
done

#seqtk trimfq -b 30 input.fastq > trimmed.fastq

###########step7: salmon对转录本定量
#所有样本中共35713个转录本的psi不为1/0/Na，太少
index=/disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.salmon.index/
for one in ${sample[@]};
do
    cleanData=/disk1/wenqing/tmp_data/PFC_s2/clean_data/"$one"/
    salmon_output=/disk1/wenqing/tmp_data/PFC_s2/salmon_output/"$one"/
    mkdir "$salmon_output"
    nohup salmon quant -i $index  -l ISF --gcBias \
        -1 "$cleanData"/GABAergic_neurons_read1_trimmed.fq -2 "$cleanData"/GABAergic_neurons_read2_trimmed.fq  -p 10 \
        -o  "$salmon_output" 1>"$salmon_output"/salmon.log 2>&1  &
 done

###########salmon对未剪切的fastq文件定量
index=/disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.salmon.index/
for one in ${sample[@]};
do
    Data=/disk1/wenqing/tmp_data/PFC_s2/data/"$one"/
    salmon_output=/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/"$one"/
    mkdir "$salmon_output"
    nohup salmon quant -i $index  -l ISF --gcBias \
        -1 "$Data"/GABAergic_neurons_read1.fastq -2 "$Data"/GABAergic_neurons_read2.fastq  -p 20 \
        -o  "$salmon_output" 1>"$salmon_output"/salmon.log 2>&1  &
done

###########step8: 整合各样本的转录本丰度并转换转录本ID
multipleFieldSelection.py \
-i  GW*/quant.sf -k 1 -f 4 \
-o iso_tpm_untrimmed_ver.txt

perl -alne '{/(\|.*\|)\t/; ;s/$1//g;s/\|//g;print}' iso_tpm_untrimmed_ver.txt > iso_tpm_formatted_untrimmed_ver.txt


ioe_merge_file=/disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.all.events.ioe
#ls $ioe_merge_file
suppa.py psiPerEvent \
-i $ioe_merge_file -e /disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/iso_tpm_formatted_untrimmed_ver.txt -o project_events_untrimmed_ver 1>psiPerEvent_log.txt 2>&1 &
# 会输出  project_events.psi 文件