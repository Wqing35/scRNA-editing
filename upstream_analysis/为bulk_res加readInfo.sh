#!/bin/bash

####################在单个细胞中call形成snv cluster####################
#先拿GW16_1_4 GABA尝试

#########step1: 整理input fastq，双端测序数据的barcode及reads提取,并输入改进后的单个细胞的pipeline中进行res的鉴定
python2 /disk1/wenqing/scRNA-editing/sprint/sprint_main_hisat2_onecell.py -rp /disk1/wenqing/tmp_data/hg19_repeat.txt -c 0 -cr 8 -p 60 -1 /disk1/wenqing/tmp_data/PFC_s2/all_fastq/SRR6367166_1.fastq -2 /disk1/wenqing/tmp_data/PFC_s2/all_fastq/SRR6367166_2.fastq -br /disk1/wenqing/tmp_data/PFC_s2/20240425_use/GW16_1_4/GABAergic_neurons_barcode.txt /disk1/wenqing/tmp_data/hg19/hg19.fa /disk1/wenqing/tmp_data/PFC_s2/result/GW16_1_4_onecell /disk1/wenqing/hisat2-2.2.0/hisat2 /disk1/wenqing/SPRINT/SPRINT_master/samtools_and_bwa/samtools

#########step2：经上述流程得到的res AG-ratio也极低，怀疑大量的snv深度太低，即测序质量不高
#在这里进一步测试，具体为：将标注了read来源的snv去掉read info，合并为bulk的形式，过滤深度小于等于10的snv，用bulk的流程得到新的结果
python2 /disk1/wenqing/scRNA-editing/sprint/sprint_main_hisat2.py -rp /disk1/wenqing/tmp_data/hg19_repeat.txt -c 0 -p 60 -1 /disk1/wenqing/tmp_data/PFC_s2/all_fastq/SRR6367166_1.fastq -2 /disk1/wenqing/tmp_data/PFC_s2/all_fastq/SRR6367166_2.fastq /disk1/wenqing/tmp_data/hg19/hg19.fa /disk1/wenqing/tmp_data/PFC_s2/result/GW16_1_4_onecell/tmp/filter_snv_depth /disk1/wenqing/hisat2-2.2.0/hisat2 /disk1/wenqing/SPRINT/SPRINT_master/samtools_and_bwa/samtools

#结果显示为，过滤snv的深度后，AG-ratio有显著的提升，因此还是按照bulk的方式call res


####################设计为celltype bulk的pipeline加入read info####################
#要将一个isoform的的iri/psi和该isoform的RNA editing level联系起来，需要按照在最开始的阶段就对每一个位点标注read信息

#########step1：先确认对isoform定量时，可以追溯read来源
index=/disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.salmon.index/
sample=(GW08 GW12 GW16_1_3 GW16_1_4 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${sample[@]};
do
    Data=/disk1/wenqing/tmp_data/PFC_s2/data/"$one"/
    salmon_output=/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/"$one"/
    nohup salmon quant -i $index -l ISF --gcBias \
            -1 "$Data"/GABAergic_neurons_read1.fastq -2 "$Data"/GABAergic_neurons_read2.fastq  -p 60 \
            -o "$salmon_output" --writeMappings "$salmon_output"/mappings.sam 1>"$salmon_output"/salmon.log 2>&1 &
done


#########step2：更改bulk流程，加入每一个snv的read info
######注意，加入readinfo的新脚本相比于bulk版本会过滤掉在比对到基因组多个位置但产出相同位置snv的read
#####两个版本生成的res数量（根据深度过滤snv后的）及AG ratio存储在本地表格——bulk与具有readInfo的res统计对比.xlsx
sample=(GW08 GW12 GW16_1_3 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${sample[@]};
do
    #在原有路径下创建子目录并工作
    input_celltype_fq_dir=/disk1/wenqing/tmp_data/PFC_s2/data/"$one"
    original_dir=/disk1/wenqing/tmp_data/PFC_s2/result/"$one"/GABAergic_neurons/tmp/
    output=/disk1/wenqing/tmp_data/PFC_s2/result/"$one"/GABAergic_neurons/tmp/readInfo_dpFiltered/tmp
    mkdir -p "$output"
    cd "$output"
    ln -s "$original_dir"/genome_all.zz.dedup ./
    ln -s "$original_dir"/transcript_all.zz.dedup ./
    ln -s "$original_dir"/all_combined.zz.sorted ./
    ln -s "$original_dir"/baseq.cutoff ./

    echo "$one"

    python2 /disk1/wenqing/scRNA-editing/sprint/sprint_main_hisat2_readInfo.py -rp /disk1/wenqing/tmp_data/hg19_repeat.txt -c 0 -p 60 -1 "$input_celltype_fq_dir"/GABAergic_neurons_read1.fastq -2 "$input_celltype_fq_dir"/GABAergic_neurons_read1.fastq /disk1/wenqing/tmp_data/hg19/hg19.fa "$output" /disk1/wenqing/hisat2-2.2.0/hisat2 /disk1/wenqing/SPRINT/SPRINT_master/samtools_and_bwa/samtools
done
   
#########step3：将res的read id和原始的read shebang行对应起来
#注意，read体量过大，先提取quant.sf中丰度不为0的isoform
sample=(GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3)
for one in ${sample[@]};
do
    cd /disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/"$one"/
    cat quant.sf |awk '{if($5!=0) print $1}'|grep -v Name | sort -k1,1 > iso_quantOver0.srt.txt
    #提取isoform和对应的序列（注意不是每一行的第十列都是序列）
    cat mappings.sam |grep -v @|awk '$10 ~ /^[AGCT]/ {print $3"\t"$10}' | sort -k1,1 > only_iso_seq.srt.txt

    #根据提取的quant.sf中的isoform id提取mapping.sam中对应的序列
    #salmon的mapping文件pcr的重复序列也在其中，所以有定量结果的iso和seq的一一对应关系有部分是重复的
    join -1 1 -2 1 -o 2.1,2.2 iso_quantOver0.srt.txt only_iso_seq.srt.txt | sort -k2,2 | uniq -c | awk '{print $2"\t"$3}' > ./seq_wz_iso_quantOver0.srt.txt
done

#提取每一个样本GABA的read id和对应的seq
sample=(GW08 GW12 GW16_1_3 GW16_1_4 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${sample[@]};
do
    cd /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/GABAergic_neurons/tmp/readInfo_dpFiltered/tmp/
    sort -k7,7 regular.res > regular.readSrt.res
    cd /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/GABAergic_neurons/tmp/genome/
    samtools view all.bam |awk '{print $1"\t"$10}' | sort -k1,1 > ./seq_wz_readID.srt.txt
    #仅保留有编辑位点富集的readID和seq
    join -1 1 -2 7 -o 1.1,1.2 seq_wz_readID.srt.txt ../readInfo_dpFiltered/tmp/regular.readSrt.res |sort -k2,2 > seq_wz_readID_hasRes.srt.txt
done


#生成iso、read id相对应的文件
sample=(GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${sample[@]};
do
    dict1=/disk1/wenqing/tmp_data/PFC_s2/result/"$one"/GABAergic_neurons/tmp/genome/
    dict2=/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/"$one"/
    #结果是多对多，一个read对应多个iso、一个iso有多read比对上去
    join -1 2 -2 2 -o 1.1,2.1 "$dict1"/seq_wz_readID_hasRes.srt.txt "$dict2"/seq_wz_iso_quantOver0.srt.txt > "$dict1"/readID_wz_iso_hasRes.txt
    cat "$dict1"/readID_wz_iso_hasRes.txt |sed 's/|/\t/g'|awk '{print $1"\t"$2}'|sort |uniq -c | awk '{print $2"\t"$3}' > "$dict1"/readID_wz_iso_hasRes.clean.txt
    rm "$dict1"/seq_wz_readID.srt.txt
done

#每个样本的salmon定量过程文件mapping.sam及后续中间文件体量过大，注意一个样本一个样本的分析






