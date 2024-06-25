###对输出的所有样本中intron对应的iri表格做格式处理
awk 'BEGIN{FS=OFS="\t"} {gsub(/-/,"@",$1); gsub(/-/,"@",$3); gsub(/-/,"@",$4); split($2,a,"-"); $2=a[1]; $3=a[2]; gsub(/@/,"-",$1); gsub(/@/,"-",$3); gsub(/@/,"-",$4); print}' all_inron_iri.txt


###拓展每个intron region上下游100个bp
#exd_len=(200 300 400 500 600 700 800 900 1000)
exd_len=(100 300 500 700 900)
for one_len in ${exd_len[@]};
do
    mkdir ./exd"$one_len"
    awk -v shift="$one_len" '{print $1"\t"$2-shift"\t"$3+shift"\t"$4"\t"$5"\t"$6}' /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/intron_retention_analysis_senescence-master/processed_all_intron_iri.srt.txt > ./exd"$one_len"/processed_all_intron_iri.srt.exd"$one_len".txt
    bedtools intersect -wo -a ./exd"$one_len"/processed_all_intron_iri.srt.exd"$one_len".txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.5/stable.srt.res > ./exd"$one_len"/intron_exd${one_len}_region_with_level_stable_res.txt
    bedtools intersect -wo -a ./exd"$one_len"/processed_all_intron_iri.srt.exd"$one_len".txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.5/up.srt.res > ./exd"$one_len"/intron_exd${one_len}_region_with_level_up_res.txt
    bedtools intersect -wo -a ./exd"$one_len"/processed_all_intron_iri.srt.exd"$one_len".txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.5/down.srt.res > ./exd"$one_len"/intron_exd${one_len}_region_with_level_down_res.txt
done


bedtools intersect -wo -a ./processed_all_intron_iri.srt.exd100.txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/de_edited_sites_clsuter1.srt.res > intron_exd100_region_with_de_edited_sites_cluster1.txt
bedtools intersect -wo -a ./processed_all_intron_iri.srt.exd100.txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/de_edited_sites_clsuter2.srt.res > intron_exd100_region_with_de_edited_sites_cluster2.txt
bedtools intersect -wo -a ./processed_all_intron_iri.srt.exd100.txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/de_editing/GABAergic_neurons/de_edited_sites_clsuter3.srt.res > intron_exd100_region_with_de_edited_sites_cluster3.txt

###获得每个intron region中的三组编辑位点
more editing_level_down_events.txt|awk '{print $1}'|sed 's/_/\t/g' |awk '{print $1"\t"$2-1"\t"$2}'|awk '{if($2!="-1") print $0}' > down.res
more editing_level_up_events.txt|awk '{print $1}'|sed 's/_/\t/g' |awk '{print $1"\t"$2-1"\t"$2}'|awk '{if($2!="-1") print $0}' > up.res
more editing_level_stable_events.txt|awk '{print $1}'|sed 's/_/\t/g' |awk '{print $1"\t"$2-1"\t"$2}'|awk '{if($2!="-1") print $0}' > stable.res


bedtools sort -faidx ~/tmp_data/names.txt -i stable.res > stable.srt.res
bedtools sort -faidx ~/tmp_data/names.txt -i up.res > up.srt.res
bedtools sort -faidx ~/tmp_data/names.txt -i down.res > down.srt.res

bedtools intersect -wo -a /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/intron_retention_analysis_senescence-master/processed_all_intron_iri.srt.exd100.txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/down.srt.res > intron_exd100_region_with_level_down_res.txt
bedtools intersect -wo -a /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/intron_retention_analysis_senescence-master/processed_all_intron_iri.srt.exd100.txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/up.srt.res > intron_exd100_region_with_level_up_res.txt
bedtools intersect -wo -a /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/intron_retention_analysis_senescence-master/processed_all_intron_iri.srt.exd100.txt -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/stable.srt.res > intron_exd100_region_with_level_stable_res.txt


###分别获取有三组res富集intron所在的gene
cat stable_iri_fc.txt|awk '{print $1}'|sed 's#/#\t#g' |awk '{print $2}' > stable_intron_genes.txt
cat up_iri_fc.txt|awk '{print $1}'|sed 's#/#\t#g' |awk '{print $2}' > up_intron_genes.txt
cat down_iri_fc.txt|awk '{print $1}'|sed 's#/#\t#g' |awk '{print $2}' > down_intron_genes.txt

##提取有iri值的gene，做GO富集分析
more down_iri_fc.txt|awk '{ if (match($NF, /^[+-]?[0-9]+([.][0-9]*)?$/)) print $1}' |sed 's#/#\t#g' |awk '{print $2}' > down_iri_genes_for_GO.txt
more up_iri_fc.txt|awk '{ if (match($NF, /^[+-]?[0-9]+([.][0-9]*)?$/)) print $1}' |sed 's#/#\t#g' |awk '{print $2}' > up_iri_genes_for_GO.txt
more stable_iri_fc.txt|awk '{ if (match($NF, /^[+-]?[0-9]+([.][0-9]*)?$/)) print $1}' |sed 's#/#\t#g' |awk '{print $2}' > stable_iri_genes_for_GO.txt
##不限制有iri值
