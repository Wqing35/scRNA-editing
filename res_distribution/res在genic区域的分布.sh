#####去掉染色体上的非常规区域，例如chr_Un7328638
awk 'NF == 4'  uniq_all_res_GABA_forAlu_distri.txt > new_GABA_regular.res
awk 'NF == 4'  uniq_all_res_neuron_forAlu_distri.txt > new_neuron_regular.res
awk 'NF == 4'  uniq_all_res_opc_forAlu_distri.txt > new_opc_regular.res


#####在refseq的基础上获得intron和intergenic区域的注释文件
#########Define exons
cat hg19.refGene.gtf |
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5,$10}' |
bedtools sort |
bedtools merge -i - -c 4 -o collapse > refGene_exon.bed

#########Define transcript
cat hg19.refGene.gtf |
awk 'BEGIN{OFS="\t";} $3=="transcript" {print $1,$4-1,$5}' |
bedtools sort |
bedtools merge -i - > ./refGene_transcript.bed

##########Define introns
#取gene.bed在exon.bed上的差集
cat ~/tmp_data/gene.bed|
awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$4}' |
bedtools sort |
bedtools subtract -a stdin -b refGene_exon.bed > refGene_intron.bed

########Define intergenic
#保证hg19.genome与经bedtools排序后的bed顺序一致
cat ~/tmp_data/gene.bed|
awk 'BEGIN{OFS="\t";} {print $1,$2,$3}' |
bedtools sort |
bedtools complement -i stdin -g ~/tmp_data/hg19.srt.genome > refGene_intergenic.bed

###########Define 3utr
cat hg19.refGene.gtf |
awk 'BEGIN{OFS="\t";} $3=="3UTR" {print $1,$4-1,$5}' |
bedtools sort |
bedtools merge -i - > refGene_3utr.bed

###########Define 5utr
cat hg19.refGene.gtf |
awk 'BEGIN{OFS="\t";} $3=="5UTR" {print $1,$4-1,$5}' |
bedtools sort |
bedtools merge -i - > refGene_5utr.bed


#
bedtools intersect -wo -a refGene_intron.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/res_overlap/filterSNVBy10_ver/new_GABA_regular.res |awk '{print $4"_"$6}'|sort|uniq -c|wc -l
bedtools intersect -wo -a refGene_intergenic.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/res_overlap/filterSNVBy10_ver/new_GABA_regular.res |awk '{print $4"_"$6}'|sort|uniq -c|wc -l
bedtools intersect -wo -a refGene_3utr.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/res_overlap/filterSNVBy10_ver/new_GABA_regular.res |awk '{print $4"_"$6}'|sort|uniq -c|wc -l
bedtools intersect -wo -a refGene_exon.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/res_overlap/filterSNVBy10_ver/new_GABA_regular.res |awk '{print $4"_"$6}'|sort|uniq -c|wc -l
bedtools intersect -wo -a refGene_5utr.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/res_overlap/filterSNVBy10_ver/new_GABA_regular.res |awk '{print $4"_"$6}'|sort|uniq -c|wc -l

26973 4472 1538 1886 14 33169 5778 2563 3041 45 5903 805 413 492 7
