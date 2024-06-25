#!/bin/bash

#########1.得到4组psi对应的gene（ENSMBL Ver）
#####经由ENSMBL database转换为gene symble
cat psi_up_events_names.txt|awk '{print $1}' |sort|uniq -c |awk '{print $2}' > psi_up_genes.txt
cat psi_down_events_names.txt|awk '{print $1}' |sort|uniq -c |awk '{print $2}' > psi_down_genes.txt
cat psi_stable_events_names.txt|awk '{print $1}' |sort|uniq -c |awk '{print $2}' > psi_stable_genes.txt


#########2.得到差异编辑位点所在的基因（3簇）
#####经由ENSMBL database转换为gene symble
cat de_edited_sites_wz_3clusters.txt|awk '{if($4==1) print $1}' > de_edited_sites_clsuter1.txt
cat de_edited_sites_wz_3clusters.txt|awk '{if($4==2) print $1}' > de_edited_sites_clsuter2.txt
cat de_edited_sites_wz_3clusters.txt|awk '{if($4==3) print $1}' > de_edited_sites_clsuter3.txt

sed -i 's/_/\t/g' de_edited_sites_clsuter1.txt
sed -i 's/_/\t/g' de_edited_sites_clsuter2.txt
sed -i 's/_/\t/g' de_edited_sites_clsuter3.txt

cat de_edited_sites_clsuter1.txt |awk '{print $1"\t"$2-1"\t"$2}' > de_edited_sites_clsuter1.res
cat de_edited_sites_clsuter2.txt |awk '{print $1"\t"$2-1"\t"$2}' > de_edited_sites_clsuter2.res
cat de_edited_sites_clsuter3.txt |awk '{print $1"\t"$2-1"\t"$2}' > de_edited_sites_clsuter3.res

bedtools sort -i de_edited_sites_clsuter1.res -faidx /disk1/wenqing/tmp_data/names.txt > de_edited_sites_clsuter1.srt.res
bedtools sort -i de_edited_sites_clsuter2.res -faidx /disk1/wenqing/tmp_data/names.txt > de_edited_sites_clsuter2.srt.res
bedtools sort -i de_edited_sites_clsuter3.res -faidx /disk1/wenqing/tmp_data/names.txt > de_edited_sites_clsuter3.srt.res

bedtools intersect -wo -a ~/tmp_data/gene.bed -b de_edited_sites_clsuter1.srt.res|awk '{print $4}'|sort|uniq -c|awk '{print $2}' > genes_de_edited_sites_cluster1.txt
bedtools intersect -wo -a ~/tmp_data/gene.bed -b de_edited_sites_clsuter2.srt.res|awk '{print $4}'|sort|uniq -c|awk '{print $2}' > genes_de_edited_sites_cluster2.txt
bedtools intersect -wo -a ~/tmp_data/gene.bed -b de_edited_sites_clsuter3.srt.res|awk '{print $4}'|sort|uniq -c|awk '{print $2}' > genes_de_edited_sites_cluster3.txt


sed -i 's/"//g' genes_de_edited_sites_cluster1.txt
sed -i 's/;//g' genes_de_edited_sites_cluster1.txt
sed -i 's/;//g' genes_de_edited_sites_cluster2.txt
sed -i 's/"//g' genes_de_edited_sites_cluster2.txt
sed -i 's/;//g' genes_de_edited_sites_cluster3.txt
sed -i 's/"//g' genes_de_edited_sites_cluster3.txt