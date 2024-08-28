############针对带有readInfo的res的分析流程

###########################各个as事件psi变化（先看RI）和对应isoform上RNA editing level的相关性+变化趋势
#################分离有编辑发生/没有编辑发生的as事件，比较两组编辑水平的高低（非对应分组的统计检验）

####step1: 准备可变剪接转录本和其他转录本的区间（附加比对到转录本上的read）
cat gencode.v37.annotation.gtf| awk '{if ($3=="transcript") print $1"\t"$4"\t"$5"\t"$16"\t"$12"\t"$18}' |sed 's/[";]//g' > ./gencode_hg19_transcript.bed
#共234485个转录本
#提取RI对应的可变事件转录本和对应的其他total转录本
cat gencode.v37.events_RI_strict.ioe|awk '{if($4!="alternative_transcripts") print $4}'|sed 's/,/\n/g' |sort|uniq -c |awk '{print $2}' > RI_as_events.txt
cat gencode.v37.events_RI_strict.ioe|awk '{if($4!="total_transcripts") print $5}'|sed 's/,/\n/g' |sort|uniq -c |awk '{print $2}' > RI_total_events.txt
#all
cat gencode.v37.all.events.ioe|awk '{if($4!="alternative_transcripts") print $4}'|sed 's/,/\n/g' |sort|uniq -c |awk '{print $2}' > All_as_events.txt
cat gencode.v37.all.events.ioe|awk '{if($4!="total_transcripts") print $5}'|sed 's/,/\n/g' |sort|uniq -c |awk '{print $2}' > All_total_events.txt
#其他类型的可变剪接事件
all_as_types=(A3 A5 AF AL MX SE)
for type in ${all_as_types[@]};
do
    cat gencode.v37.events_"$type"_strict.ioe|awk '{if($4!="alternative_transcripts") print $4}'|sed 's/,/\n/g' |sort|uniq -c |awk '{print $2}' > "$type"_as_events.txt
    cat gencode.v37.events_"$type"_strict.ioe|awk '{if($4!="total_transcripts") print $5}'|sed 's/,/\n/g' |sort|uniq -c |awk '{print $2}' > "$type"_total_events.txt
done

#准备对应的bed文件——用R
#先整体看关系，在分别提不同的转录本看趋势


###################step2: 提取比对到每个转录本的read id——用python对每一个转录本进行遍历，累加res type为AG的as数量
####step2.1: 计算每个转录本的res数量
#具体实现参见../scRNA-editing/upstream_analysis/为bulk_res加readInfo.sh
#../scRNA-editing/downstream_analysis/iso_res_analysis/计算每个转录本上res的数量.ipynb

##注意！！！！！！！！！！！！！！！！！！！
##鉴定res的比对结果（每个res所在的read）与salmon的比对位置（每个转录本对应的read）不一致！！！！！！！
#怀疑是将fastq作为input输入salmon，收到salmon quasi-mapping策略的影响，hisat2与其比对的结果不一致
#先拿一个数据量大的样本（GW26_1_1）尝试看是否有当前样本数据量不够的原因

####step2.2: 计算每个转录本的Editing level(计算深度时注意限制在比对到当前转录本的read上)
#带有readInfo的zz文件使用regular_combined.zz

#保留res数量不为0的转录本, 并转换为计算EI所需的格式
directory=/disk1/wenqing/tmp_data/PFC_s2/result/GW26_1_1/GABAergic_neurons/tmp/readInfo_dpFiltered/tmp
cat "$directory"/RI_total_events_wzResNum.txt|awk '{if($7!=0) print $1"\t"$2"\t"$3"\t"$5"\t"$7}' > RI_total_events_wzResNum_Over0.txt
cat "$directory"/RI_as_events_wzResNum.txt|awk '{if($7!=0) print $1"\t"$2"\t"$3"\t"$5"\t"$7}' > RI_as_events_wzResNum_Over0.txt
#通过readInfo限制纳入计算深度的读段时，转录本无法得到有效（一致）的区间
##经all.bam和mapping.sam任意read的测试验证这一点——用‘GAGACCATCCTGGCTAATACGGTGAAACCCCGTCTCTACTAAAAATACAAAAGATTAGCCGGGTGAGGTGGCATGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGTAGGCGGAGCTCGCAGTGAG’
#这个seq为例，在bam文件中该序列比对到chr2上，而在salmon的mapping.sam文件中，该序列比对到ENST00000334574.12转录本上，据注释文件、该转录本位于chr15上


################换用使用hisat2进行比对并定量转录本的工具——featureCounts
#cat gencode.v37.annotation.gtf| awk '{if ($3=="transcript") print $0}' > ./gencode.v37.transcripts.gtf
#bam_file='/disk1/wenqing/tmp_data/PFC_s2/result/GW26_1_1/GABAergic_neurons/tmp/genome/all.srt.0000.bam'
#featureCounts -T 20 -a /disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.gtf -o ./counts.txt -O -p "$bam_file"

################stringtie: 无法生成read匹配转录本的信息，根据每条read
#stringtie -G /disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.annotation.gtf -o ./test.gtf -A test_abund.txt -p 20 --bam_output test_tracking.bam -e "$bam_file"
#stringtie inspect -t ENST00000538951 your_prefix_tracking.bam

################使用hisat2的RNA-seq比对到基因组的bam结果作为input输入到salmon中进行定量
#不可行！因为salmon定量所需的gtf注释文件来源于转录组、而bam文件为基因组比对结果
bam_file='/disk1/wenqing/tmp_data/PFC_s2/result/GW26_1_1/GABAergic_neurons/tmp/genome/all.bam'
salmon_output='/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_bam_ver'
nohup salmon quant -t /disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.fa -l ISF -a "$bam_file" -o "$salmon_output" -p 60 &

################使用hisat2将RNA-seq比对到转录组的bam结果作为input
#先比对到转录组
read1=/disk1/wenqing/tmp_data/PFC_s2/data/GW26_1_1/GABAergic_neurons_read1.fastq
read2=/disk1/wenqing/tmp_data/PFC_s2/data/GW26_1_1/GABAergic_neurons_read2.fastq
refgenome=/disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.fa
hisat2 -q -x "$refgenome" -U "$read1" -S name_read1.sam
hisat2 -q -x "$refgenome" -U "$read2" -S name_read2.sam
#检查与比对到基因组相比，两者的比对信息是否能对应上——抽检100条序列，看两者比对的区间是否一致
seqtk sample -s123 GABAergic_neurons_read2.fastq 100 > GABA_random_100_seq.fastq
#使用R统计这一百条序列在genome和transcript的比对情况

#salmon定量比对到转录组上的bam文件
bam_file='/disk1/wenqing/tmp_data/PFC_s2/transcript_mapping/all.bam'
salmon_output='/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_bam_ver'
nohup salmon quant -t /disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.fa -l ISF -a "$bam_file" -o "$salmon_output" -p 60 --writeMappings "$salmon_output"/mappings.sam &

#生成read和转录本的对应关系
#比对到转录组上的bam文件没有对应的readID，但reads的比对顺序和比对到genome的顺序一致（即read和readID的情况按genome的结果即可）
#对每一个转录本，先从bam文件获得对应的seq，再根据seq获取对应的readID
samtools view all.bam |awk '{if($2!=4) print $0}' |sort -k3,3 > all.sam
#根据quant.sf获取tpm不为0的转录本ID
cat quant.sf |awk '{if($4!=0 && $4!="TPM") print $1"\t"$4}'|sort -k1,1 > iso_wz_tpmOver0.srt.txt
#根据all.sam比对结果和iso_wz_tpmOver0.srt.txt交叉获得比对到转录本ID和seq的对应信息
join -1 1 -2 3 -o 1.1,1.2,2.10 iso_wz_tpmOver0.srt.txt /disk1/wenqing/tmp_data/PFC_s2/transcript_mapping/all.sam > ./iso_wz_seqAndTPMOver0.txt
sort -k3,3 ./iso_wz_seqAndTPMOver0.txt > ./iso_wz_seqAndTPMOver0.srt.txt
#根据readID和seq的对应关系，将上述./iso_wz_seqAndTPMOver0.srt.txt转换为转录本ID和readID的关系
join -1 3 -2 2 -o 1.1,1.2,2.1 ./iso_wz_seqAndTPMOver0.srt.txt /disk1/wenqing/tmp_data/PFC_s2/result/GW26_1_1/GABAergic_neurons/tmp/genome/seq_wz_readID_hasRes.srt.txt |sort|uniq -c |awk '{print $2"\t"$3"\t"$4}'|sort -k1,1 > ./iso_wz_readAndTPMOver0.txt

#随机抽查几个readID，检查对应的基因组范围和转录本坐标是否一致
#从比对到转录组的结果中随机抽取100个比对成功的序列
cat all.srt.sam|awk '{print $10}' |shuf -n 100 |sort -k1,1 > 100_seqs.srt.txt
#根据100个序列提取比对到基因组上的结果
join -1 10 -2 1 -o 1.1,1.2,1.3,1.4,1.10 /disk1/wenqing/tmp_data/PFC_s2/result/GW26_1_1/GABAergic_neurons/tmp/genome/all.srt.sam 100_seqs.srt.txt  
#100个序列全部比对到转录组上，却只有一半比对到基因组上
#查看比对结果（主要是坐标问题）
#结果明显不一致!!!!!!!!!!!!

###########转录组位置和基因组位置对应原理清楚了
#但：salmon使用fastq的比对结果和hisat2比对结果不一致
#还是从bam开始
#新思路：不走参考基因组比对结果，根据转录组的比对结果


########################################################实际运行脚本在这里################################################
################################################################################################################

######step1:将转录组的sam文件转换为zz，方便后续计算EI
#注意转换为zz时，read的区间计算问题(在zz转换时更改还是在寻找对应transcript的reads时更改？)
#!!!!!具体见/disk1/wenqing/scRNA-editing/calculate_EI/transcript_sam2zz.py脚本
index=/disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.salmon.index/
sample=(GW12 GW16_1_3 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${sample[@]};
do
    input_dir=/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/"$one"
    Data=/disk1/wenqing/tmp_data/PFC_s2/data/"$one"/
    salmon_output=/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/"$one"/
    salmon quant -i $index -l ISF --gcBias \
            -1 "$Data"/GABAergic_neurons_read1.fastq -2 "$Data"/GABAergic_neurons_read2.fastq  -p 60 \
            -o "$salmon_output" --writeMappings "$salmon_output"/mappings.sam 1>"$salmon_output"/salmon.log 2>&1 
    cat "$input_dir"/mappings.sam |grep -v @ > "$input_dir"/pure_mappings.sam
    rm "$input_dir"/mappings.sam
done
python2 /disk1/wenqing/scRNA-editing/calculate_EI/transcript_sam2zz.py


######step2:根据bam文件确定seq、readID和转录本的对应关系
join -1 10 -2 2 -o 2.1,1.3 /disk1/wenqing/tmp_data/PFC_s2/transcript_mapping/all.srt.sam /disk1/wenqing/tmp_data/PFC_s2/result/GW26_1_1/GABAergic_neurons/tmp/genome/seq_wz_readID_hasRes.srt.txt |sed 's/|/\t/g' |awk '{print $1"\t"$2}' > readID_wz_iso_hasRes.clean.verFromBam.txt
#从比对到转录组的bam得到的已定量转录本和readID的对应数量远少于从fq出发的结果
#两者的区别是：使用salmon的quasi-mapping还是hisat2
#先看相关性（后比较从fq开始的结果
#！！！！！！！部分iso有res富集，EI却计算为NA——genome mapping和转录组mapping的不一致
#！！！！！！！as & total events的EI和tpm呈负相关，p值不显著

#！！！！！！！因此：改成fq的转录本定量
#fq开始的转录本定量：计算EI时没有NA的出现
#！！！！但结果是相关性不成趋势，p值也不显著


#获取所有类型没有发生可变剪接的转录本
as_types=(A3 A5 AF AL MX SE RI)
for one in ${as_types[@]};
do
    echo "$one"
    #grep -v -f "$one"_as_events_wzResNum_Over0.wzEI.txt "$one"_total_events_wzResNum_Over0.wzEI.txt > "$one"_nonAs_events_wzResNum_Over0.wzEI.txt
    grep -v -f "$one"_as_events.bed "$one"_total_events.bed > "$one"_nonAs_events.bed
done

#####把其他的样本跑了先，再计算这些转录本psi/IRI和EI的关系




###另外构建转录组的参考文件（暂时不管）
bedtools getfasta -fi /disk1/wenqing/tmp_data/hg19/hg19.fa -bed gencode_hg19_transcript.bed -fo hg19_transcript.fa
awk 'BEGIN{FS="\t"; OFS="\t"} NR==FNR{names[$1]=$2; next} $0 ~ />/{print ">" names[substr($0,2)]; next}{print}' names_mapping_info.uniq.txt hg19_transcript.fa > hg19_transcript.renamed.fa
hisat2-build -p 20 ../SUPPA2/ref/hg19_transcript.renamed.fa hisat2_trans_ver2