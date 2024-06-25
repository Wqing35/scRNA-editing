######扩大lncRNA的注释
######扩大lncRNA的注释
##所有结果数字分别为GABA、Neuron、OPC

#1. LNCipedia——33
bedtools intersect -wo -a LNCipedia_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res > lncRNA_all_edited_res.txt
cat lncRNA_all_edited_res.txt|awk '{print $4"_"$6}'|sort|uniq -c|awk '{print $2}' > all_sites_lncRNA.res
sort -k1 all_sites_lncRNA.res > all_sites_lncRNA.srt.res
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt all_sites_lncRNA.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.res

#2.FANTOM5——91
https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv1_raw/FANTOM_CAT.lv1_raw.gtf.gz
https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv1_raw/FANTOM_CAT.lv1_raw.only_lncRNA.gtf.gz
https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv1_raw/FANTOM_CAT.lv1_raw.only_mRNA.gtf.gz

#all_raw.only.lncRNA.gtf中有这几种lncRNA类型
#  7100 "lncRNA_antisense";
#  14865 "lncRNA_divergent";
#  44791 "lncRNA_intergenic";
#  8031 "lncRNA_sense_intronic";


more FANTOM_CAT.lv1_raw.only_lncRNA.gtf|awk '{if($3=="gene") print $1"\t"$4"\t"$5"\t"$10}' > FANTOM_lncRNA_gene.bed

more FANTOM_CAT.lv1_raw.only_mRNA.gtf |awk '{if($3=="gene") print $1"\t"$4"\t"$5"\t"$10}' > FANTOM_mRNA_gene.bed


bedtools sort -faidx /disk1/wenqing/tmp_data/names.txt -i ./FANTOM_lncRNA_gene.bed > ./FANTOM_lncRNA_gene.srt.bed
bedtools sort -faidx /disk1/wenqing/tmp_data/names.txt -i ./FANTOM_mRNA_gene.bed > ./FANTOM_mRNA_gene.srt.bed

bedtools intersect -wo -a ./FANTOM_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res > lncRNA_all_edited_res.txt
#3087
bedtools intersect -wo -a ./FANTOM_mRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res > mRNA_all_edited_res.txt
#6983

######落在lncRNA splice sites的位点
#1.先提所有在lncRNA上的位点
more lncRNA_all_edited_res.txt|awk '{print $5"_"$7}'|sort|uniq -c |awk '{print $2}' > all_sites_lncRNA.res
sort -k1,1 all_sites_lncRNA.res > all_sites_lncRNA.srt.res

#2.与spice位点做交叉
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt all_sites_lncRNA.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.res



#3.MiTranscriptome——7
bedtools intersect -wo -a ./lncRNA.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res > lncRNA_all_edited_res.txt

#4.BIGTranscriptome——0
bedtools intersect -wo -a BIGTranscriptome_lncRNA_gene.processed.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.processed.srt.res > lncRNA_all_edited_res.txt
cat lncRNA_all_edited_res.txt|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}'|more
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt all_sites_lncRNA.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.res



########合并统计4个database的结果——109，3646
cat LNCipedia/SpliceSites_lncRNA.uniq.res BIGTranscriptome/SpliceSites_lncRNA.uniq.res MiTranscriptome/SpliceSites_lncRNA.uniq.res FANTOM/SpliceSites_lncRNA.uniq.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.in4database.res
cat LNCipedia/all_sites_lncRNA.srt.res BIGTranscriptome/all_sites_lncRNA.srt.res MiTranscriptome/all_sites_lncRNA.srt.res FANTOM/all_sites_lncRNA.srt.res|sort|uniq -c|awk '{print $2}' > all_sites_lncRNA.in4database.res


##注意：这样统计lncRNA上位于非剪接位点的res数量比位于剪接位点的res数量多，
#尝试在差异编辑的res中进行上述分析




####################和差异编辑位点的交叉
#1. LNCipedia——33
cd ../LNCipedia
bedtools intersect -wo -a LNCipedia_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res > lncRNA_de_edited_res.txt
cat lncRNA_de_edited_res.txt|awk '{print $4"_"$6}'|sort|uniq -c|awk '{print $2}' > de_sites_lncRNA.res
sort -k1 de_sites_lncRNA.res > de_sites_lncRNA.srt.res
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt de_sites_lncRNA.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.de.res

#2.FANTOM5——91
cd ../FANTOM
bedtools intersect -wo -a ./FANTOM_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res > lncRNA_de_edited_res.txt
cat lncRNA_de_edited_res.txt|awk '{print $5"_"$7}'|sort|uniq -c |awk '{print $2}' > de_sites_lncRNA.res
sort -k1,1 de_sites_lncRNA.res > de_sites_lncRNA.srt.res
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt de_sites_lncRNA.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.de.res

#3.MiTranscriptome——7
cd ../MiTranscriptome
bedtools intersect -wo -a ./lncRNA.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res > lncRNA_de_edited_res.txt
cat lncRNA_de_edited_res.txt|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}' > de_sites_lncRNA.res
sort -k1,1 de_sites_lncRNA.res > de_sites_lncRNA.srt.res
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt de_sites_lncRNA.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.de.res


#4.BIGTranscriptome——0
cd ../BIGTranscriptome
bedtools intersect -wo -a BIGTranscriptome_lncRNA_gene.processed.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_de_edited.srt.res > lncRNA_de_edited_res.txt
cat lncRNA_de_edited_res.txt|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}'  > de_sites_lncRNA.res
sort -k1,1 de_sites_lncRNA.res > de_sites_lncRNA.srt.res
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt de_sites_lncRNA.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.de.res



########合并统计4个database的结果——109，3646
cat LNCipedia/SpliceSites_lncRNA.uniq.de.res BIGTranscriptome/SpliceSites_lncRNA.uniq.de.res MiTranscriptome/SpliceSites_lncRNA.uniq.de.res FANTOM/SpliceSites_lncRNA.uniq.de.res |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.in4database.de.res
cat LNCipedia/de_sites_lncRNA.srt.res BIGTranscriptome/de_sites_lncRNA.srt.res MiTranscriptome/de_sites_lncRNA.srt.res FANTOM/de_sites_lncRNA.srt.res|sort|uniq -c|awk '{print $2}' > de_sites_lncRNA.in4database.res


###差异编辑位点位于剪接位点？？？？？？？
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/tmp.de.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites.de.res
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/all_edited.srt.res |sort|uniq -c|awk '{print $2}' > SpliceSites.all.res


#####################不能证明lncRNA上的res更倾向于落在剪接位点（plos文章）
#####转向分析三个group中的res与剪接的关系（扮演的角色）
cat LNCipedia/LNCipedia_lncRNA_gene.srt.bed MiTranscriptome/lncRNA.srt.bed BIGTranscriptome/BIGTranscriptome_lncRNA_gene.processed.bed FANTOM/FANTOM_lncRNA_gene.srt.bed |sort|uniq -c|awk '{print $2"\t"$3"\t"$4}' > 4database_lncRNA_gene.srt.bed

##down组
bedtools intersect -wo -a 4database_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res > lncRNA_down.res
cat lncRNA_down.res|awk '{print $4"_"$6}'|sort|uniq -c|awk '{print $2}' > down_res_in_lncRNA.txt
sort -k1 down_res_in_lncRNA.txt > down_res_in_lncRNA.srt.txt
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt down_res_in_lncRNA.srt.txt |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.down.res

##up组
bedtools intersect -wo -a 4database_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res > lncRNA_up.res
cat lncRNA_up.res|awk '{print $4"_"$6}'|sort|uniq -c|awk '{print $2}' > up_res_in_lncRNA.txt
sort -k1 up_res_in_lncRNA.txt > up_res_in_lncRNA.srt.txt
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt up_res_in_lncRNA.srt.txt |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.up.res

##stable组
bedtools intersect -wo -a 4database_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res > lncRNA_stable.res
cat lncRNA_stable.res|awk '{print $4"_"$6}'|sort|uniq -c|awk '{print $2}' > stable_res_in_lncRNA.txt
sort -k1 stable_res_in_lncRNA.txt > stable_res_in_lncRNA.srt.txt
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt stable_res_in_lncRNA.srt.txt |sort|uniq -c|awk '{print $2}' > SpliceSites_lncRNA.uniq.stable.res

###########三组res在mRNA上的数量
#FANTOM注释相对较全，
bedtools intersect -wo -a ./FANTOM/FANTOM_mRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res > mRNA_down.res
cat mRNA_down.res|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}' > down_res_in_mRNA.txt
sort -k1 down_res_in_mRNA.txt > down_res_in_mRNA.srt.txt
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt down_res_in_mRNA.srt.txt |sort|uniq -c|awk '{print $2}' > SpliceSites_mRNA.uniq.down.res

##up组
bedtools intersect -wo -a ./FANTOM/FANTOM_mRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res > mRNA_up.res
cat mRNA_up.res|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}' > up_res_in_mRNA.txt
sort -k1 up_res_in_mRNA.txt > up_res_in_mRNA.srt.txt
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt up_res_in_mRNA.srt.txt |sort|uniq -c|awk '{print $2}' > SpliceSites_mRNA.uniq.up.res

##stable组
bedtools intersect -wo -a ./FANTOM/FANTOM_mRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res > mRNA_stable.res
cat mRNA_stable.res|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}' > stable_res_in_mRNA.txt
sort -k1 stable_res_in_mRNA.txt > stable_res_in_mRNA.srt.txt
join -1 1 -2 1 ~/SPIDEX/sites_all_splice.txt stable_res_in_mRNA.srt.txt |sort|uniq -c|awk '{print $2}' > SpliceSites_mRNA.uniq.stable.res



##########lncRNA和mRNA的注释区域有重复？？？？？？！！！！！
#lncRNA仅使用FANTOM的
bedtools intersect -wo -a FANTOM_lncRNA_gene.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res > lncRNA_down.res
cat lncRNA_down.res|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}' > down_res_in_lncRNA.txt
sort -k1 down_res_in_lncRNA.txt > down_res_in_lncRNA.srt.txt


##########分成mRNA_intron, mRNA_exon, lncRNA_intron, lincRNA_exon四组

#########Define exons
cat FANTOM_CAT.lv1_raw.only_lncRNA.gtf |
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5,$10}' |
bedtools sort |
bedtools merge -i - -c 4 -o collapse > FANTOM_lncRNA_exon.bed

##########Define introns
#取gene.bed在exon.bed上的差集
cat FANTOM_lncRNA_gene.bed|
awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$4}' |
bedtools sort |
bedtools subtract -a stdin -b FANTOM_lncRNA_exon.bed > FANTOM_lncRNA_intron.bed

cat FANTOM_CAT.lv1_raw.only_mRNA.gtf |
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5,$10}' |
bedtools sort |
bedtools merge -i - -c 4 -o collapse > FANTOM_mRNA_exon.bed

##########Define introns
#取gene.bed在exon.bed上的差集
cat FANTOM_mRNA_gene.bed|
awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$4}' |
bedtools sort |
bedtools subtract -a stdin -b FANTOM_mRNA_exon.bed > FANTOM_mRNA_intron.bed


#lncRNA-exon
bedtools intersect -wo -a /disk1/wenqing/tmp_data/lncRNA_anno/FANTOM/FANTOM_lncRNA_exon.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res > lncRNA_exon_down.res
cat lncRNA_exon_down.res|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}' > down_res_in_lncRNA_exon.txt
sort -k1 down_res_in_lncRNA_exon.txt > down_res_in_lncRNA_exon.srt.txt

#lncRNA-intron
bedtools intersect -wo -a /disk1/wenqing/tmp_data/lncRNA_anno/FANTOM/FANTOM_lncRNA_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res > lncRNA_intron_down.res
cat lncRNA_intron_down.res|awk '{print $5"_"$7}'|sort|uniq -c|awk '{print $2}' > down_res_in_lncRNA_intron.txt
sort -k1 down_res_in_lncRNA_intron.txt > down_res_in_lncRNA_intron.srt.txt

