#####将不同cluster的res分离出来，并富集到基因上
sed -i 's/_/\t/g' de_edited_sites_wz_3clusters.txt
cat de_edited_sites_wz_3clusters.txt |awk '{if($5==1  && NF==5) print $1"\t"$2-1"\t"$2"\t"$3}' > de_edited_in_cluster1.res
cat de_edited_sites_wz_3clusters.txt |awk '{if($5==2  && NF==5) print $1"\t"$2-1"\t"$2"\t"$3}' > de_edited_in_cluster2.res
cat de_edited_sites_wz_3clusters.txt |awk '{if($5==3  && NF==5) print $1"\t"$2-1"\t"$2"\t"$3}' > de_edited_in_cluster3.res

#enrich to gene
bedtools intersect -wo -a ~/tmp_data/gene.bed -b de_edited_in_cluster1.res > de_edited_res_of_cluster1_inGenes.txt
bedtools intersect -wo -a ~/tmp_data/gene.bed -b de_edited_in_cluster2.res > de_edited_res_of_cluster2_inGenes.txt
bedtools intersect -wo -a ~/tmp_data/gene.bed -b de_edited_in_cluster3.res > de_edited_res_of_cluster3_inGenes.txt

#####不分簇，对全体de edited sites富集到基因上
cat de_edited_sites_wz_3clusters.txt |awk '{if($1!="all" && NF==5) print $1"\t"$2-1"\t"$2"\t"$3}' > de_edited_in_all_clusters.res
bedtools intersect -wo -a ~/tmp_data/gene.bed -b de_edited_in_all_clusters.res > de_edited_res_of_all_clusters_inGenes.txt



cut -d'\t' -f4 all_edited_sites_with_p_values.txt
awk '{if (NR==1) {min=$2; max=$2;} else {min=(min<$2)?min:$2; max=(max>$2)?max:$2;}} END {print "Range: ", max-min}' data.txt

