#!/bin/bash

#####鉴定trascripts的类型并计算lncRNA和mRNA的editing rate

##step1: 整理de sites备用
cat all_edited_sites_with_p_values.txt|awk '{if($14 < 0.05) print $1}'|sed 's/_/\t/g'|awk '{print $1"\t"$2-1"\t"$2}' > all_de_edited.res
bedtools sort -faidx ~/tmp_data/names.txt -i all_de_edited.res > all_de_edited.srt.res



##step2: 使用gtf文件，获得每一个lncRNA/mRNA intron对应的区间
#使用gencode和noncode联合注释，（注意每一个gene-id先经ENSMBL转换为gene name和gene type）
bedtools sort -faidx ~/tmp_data/names.txt -i miRNA.withNoname.bed > miRNA.withNoname.srt.bed
bedtools sort -faidx ~/tmp_data/names.txt -i mRNA.withNoname.noMT.bed > mRNA.withNoname.noMT.srt.bed
bedtools sort -faidx ~/tmp_data/names.txt -i lncRNA.withNoname.bed > lncRNA.withNoname.srt.bed

bedtools sort -faidx ~/tmp_data/names.txt -i mRNA.noMT.bed > mRNA.noMT.srt.bed


##step3: 三类区间和差异编辑位点做交集
######de 编辑位点太少（666个），富集不到这些区域(加上没有对应gene name的区域也富集不到编辑位点)
#####Neurons共8322个差异编辑位点，也富集不到这些区域
bedtools intersect -wo -a lncRNA.withNoname.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res |more
bedtools intersect -wo -a miRNA.withNoname.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res |more
bedtools intersect -wo -a mRNA.withNoname.noMT.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res |more

##step4: 不仅限于差异编辑位点，在任意一个样本出现的编辑位点均纳入范围，在lncRNA中能富集到少数，miRNA和mRNA仍然无结果
##？？？？？？？：在形成编辑水平矩阵时，去除了大部份编辑水平是1的位点，
cat all_edited_sites_with_p_values.txt|awk '{print $1}' > all_edited.res
sed 's/_/\t/g' all_edited.res|awk '{if(NF==2) print $0}'|awk '{print $1"\t"$2-1"\t"$2}' > all_edited.processed.res
bedtools sort -faidx ~/tmp_data/names.txt -i all_edited.processed.res > all_edited.processed.srt.res

bedtools intersect -wo -a lncRNA.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res |more
bedtools intersect -wo -a miRNA.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res |more
bedtools intersect -wo -a mRNA.noMT.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res |more

###gencode的注释没有intron区间信息，用gene.bed富集差异编辑位点有结果，大部份编辑位点都落在intron区域
##所以尝试用lncRNA、miRNA、mRNA的intron区间富集res??????
##不对，所用三者bed文件已经是包含intron的gene bed文件
bedtools intersect -wo -a ./refGene_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res |awk '{print $4"_"$6}'|sort|uniq -c|wc -l
bedtools intersect -wo -a ./refGene_exon.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res |awk '{print $4"_"$6}'|sort|uniq -c|wc -l


######在lncRNA和mRNA上的res（无论是整个基因还是intron区）就是比较少？？？？？
##那这些res都落在哪些区域？？？
samples=(GW08 GW12 GW16_1_3 GW16_1_4 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${samples[@]};
do
    echo $one
    cat /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/GABAergic_neurons/regular.res.depth.orig|wc -l
    #bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/refseq/refGene_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/GABAergic_neurons/regular.res.depth.orig|awk '{print $5"_"$7}'|sort|uniq -c|wc -l
    bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/refseq/refGene_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res |awk '{print $5"_"$7}'|sort|uniq -c|wc -l
    #bedtools intersect -wo -a mRNA.noMT.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/GABAergic_neurons/regular.res.depth.orig|awk '{print $5"_"$7}'|sort|uniq -c|wc -l
done

##先根据refseq intron区域的注释得到有res富集的基因，再看这些基因是什么类型
bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/refseq/refGene_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res |awk '{print $4}'|sort|uniq -c|awk '{print $2}' |sed 's/"//g' | sed 's/;//g' > genes_with_de_edited_res.txt
bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/refseq/refGene_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res > detailed_info_genes_with_de_edited_res.txt

#不仅限于差异编辑位点
bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/refseq/refGene_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res > detailed_info_genes_with_all_edited_res.txt



cat hg19_spidex.txt |awk '{if($7 < 0)print "chr"$1"_"$3}' > sites_pro_splice.txt


sort -k1,1 sites_depress_splice.txt > sites_depress_splice.srt.txt
sort -k1,1 sites_pro_splice.txt > sites_pro_splice.srt.txt

#####所有落在剪接位点上的res
sort -k1,1 all_edited.res > all_edited.srt.res

join -1 1 -2 1 ~/SPIDEX/sites_depress_splice.srt.txt all_edited.srt.res > depress_SpliceSites.res
join -1 1 -2 1 ~/SPIDEX/sites_pro_splice.srt.txt all_edited.srt.res > pro_SpliceSites.res

#####落在lncRNA/mRNA分子剪接位点上的res
cat detailed_info_genes_with_all_edited_res_with_geneType.txt|awk '{if($9=="lncRNA") print $6"_"$8}'|sort|uniq -c|awk '{print $2}' > in_lncRNA.res
cat detailed_info_genes_with_all_edited_res_with_geneType.txt|awk '{if($9=="protein_coding") print $6"_"$8}'|sort|uniq -c|awk '{print $2}' > in_mRNA.res

sort -k1,1 in_lncRNA.res > in_lncRNA.srt.res
sort -k1,1 in_mRNA.res > in_mRNA.srt.res
###在各分子上的res，主要和促进剪接的位点做交集
join -1 1 -2 1 ~/SPIDEX/sites_pro_splice.srt.txt in_lncRNA.srt.res > lncRNA_pro_SpliceSites.res
join -1 1 -2 1 ~/SPIDEX/sites_pro_splice.srt.txt in_mRNA.srt.res > mRNA_pro_SpliceSites.res
