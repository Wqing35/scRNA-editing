#####计算10X/S2的res分布是否有分布上的bias
####knowGene ver
#10X
bedtools intersect -wo -a /disk1/wenqing/tmp_data/transcriptsTabByUTR_3end_sorted.txt -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/10X/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_3end.txt
bedtools intersect -wo -a /disk1/wenqing/tmp_data/transcriptsTabByUTR_5end_sorted.txt -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/10X/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_5end.txt

cat transcriptWithRegular_res_5end.txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_5end.txt
cat transcriptWithRegular_res_3end.txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_3end.txt


#S2
bedtools intersect -wo -a /disk1/wenqing/tmp_data/transcriptsTabByUTR_3end_sorted.txt -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/S2/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_3end.txt
bedtools intersect -wo -a /disk1/wenqing/tmp_data/transcriptsTabByUTR_5end_sorted.txt -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/S2/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_5end.txt

cat transcriptWithRegular_res_5end.txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_5end.txt
cat transcriptWithRegular_res_3end.txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_3end.txt


####refGene Ver
percents=(0.5 0.6 0.7 0.8)
for one in ${percents[@]};
do
    ##10X
    cd /disk1/wenqing/tmp_data/pbmc/distri_bias_test/gencode/10X
    bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/gencode/gencode_5end."$one".bed -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/10X/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_5end."$one".txt
    bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/gencode/gencode_3end."$one".bed -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/10X/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_3end."$one".txt

    cat transcriptWithRegular_res_5end."$one".txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_5end."$one".txt
    cat transcriptWithRegular_res_3end."$one".txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_3end."$one".txt

    ##S2
    cd /disk1/wenqing/tmp_data/pbmc/distri_bias_test/gencode/S2
    bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/gencode/gencode_5end."$one".bed -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/S2/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_5end."$one".txt
    bedtools intersect -wo -a /disk1/wenqing/tmp_data/gtf_files/gencode/gencode_3end."$one".bed -b /disk1/wenqing/tmp_data/pbmc/GEI/over0/S2/regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > transcriptWithRegular_res_3end."$one".txt

    cat transcriptWithRegular_res_5end."$one".txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_5end."$one".txt
    cat transcriptWithRegular_res_3end."$one".txt|awk '{print $4}'|sort|uniq -c > transcriptWzRes_3end."$one".txt
done

#########比较5end和3end的res数量的区别（ad之和）
percents=(0.5 0.6 0.7 0.8)
for one in ${percents[@]};
do
    awk -F'\t' '{ sum[$4] += $7 } END { for (i in sum) print sum[i], i }' ../transcriptWithRegular_res_3end."$one".txt > transcriptWzRes_3end."$one".txt
    awk -F'\t' '{ sum[$4] += $7 } END { for (i in sum) print sum[i], i }' ../transcriptWithRegular_res_5end."$one".txt > transcriptWzRes_5end."$one".txt
done