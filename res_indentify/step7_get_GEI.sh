#!/bin/bash

#计算Gene区域的Editing index，不局限于Alu区域
#具体流程待写，和step6计算AEI的流程相似

#samples=(GW08 GW12 GW16_1_3 GW16_1_9 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_3 GW26_1_2 GW26_1_3 GW26_1_4 GW26_1_5 GW26_1_6 GW26_1_7 GW26_1_8 GW26_1_9 GW26_1_10)
samples=(GW16_1_3 GW16_1_4 GW16_1_9 GW19_1_1 GW19_1_2 GW19_1_3 GW23_1_1 GW23_1_2 GW23_1_3 GW26_1_1)
for one in ${samples[@]};
do
    #1.
    cp /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/regular.res.depth /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/regular.res.depth.orig
    sed -i 's/:/\t/g' /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/regular.res.depth.orig
    bedtools intersect -wo -a ~/tmp_data/gene.bed -b /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/regular.res.depth.orig | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/AEI/OPC/"$one"/geneWithRegular_res.txt
    #2.手动运行/disk1/wenqing/scRNA-editing/calculate_EI/get_bulk_geneWzRes_ad.ipynb
    #3.
    cat /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/genome_all.zz.dedup /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/transcript_all.zz.dedup.genome.zz > /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/regular_combined.zz
    ln -s /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/regular_combined.zz /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/AEI/OPC/"$one"/
    cd /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/AEI/OPC/"$one"/
    #4.
    python2 /disk1/wenqing/scRNA-editing/calculate_EI/get_EI_inGene.py
    echo "$one" finished!
done



#    cp /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/regular.res.depth /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/regular.res.depth.orig
#    sed -i 's/:/\t/g' /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/regular.res.depth.orig
#    cp /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/regular.snv /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/
#    cp /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/regular.snv /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/
#    #cp /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/regular.snv /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/
#    ln -s /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/all_combined.zz.sorted /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/
#    ln -s /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/all_combined.zz.sorted /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/
#    #ln -s /disk1/wenqing/tmp_data/PFC_s2/result/"$one"/OPC/tmp/all_combined.zz.sorted /disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/"$one"/OPC/
#done