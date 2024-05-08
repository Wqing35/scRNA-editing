#!/bin/bash

#####计算celltype res的ALU Editing Index

#将所需输入文件cp到指定计算路径
input_samples=(GW16_1_9 GW19_1_1 GW23_1_2 GW26_1_1)
input_celltypes=(GABAergic_neurons OPC Neurons)
for one_sample in ${input_samples[@]};
do
    echo "$one_sample" running!
    for one_celltype in ${input_celltypes[@]}
    do
        echo "$one_celltype" running!
        data_directory=/disk1/wenqing/tmp_data/PFC_s2/result/"$one_sample"/"$one_celltype"/tmp
        AEI_directory=/disk1/wenqing/tmp_data/PFC_s2/AEI/"$one_celltype"/"$one_sample"

        cat "$data_directory"/genome_all.zz.dedup "$data_directory"/transcript_all.zz.dedup.genome.zz > "$data_directory"/regular_combined.zz
        cp "$data_directory"/regular.res.depth "$AEI_directory"/regular_merged.res
        ln -s "$data_directory"/regular_combined.zz "$AEI_directory"

        cd "$AEI_directory"
        python2 /disk1/wenqing/scRNA-editing/calculate_EI/get_AEI_inAlu.py
    done
done

