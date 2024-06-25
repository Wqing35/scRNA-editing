#!/bin/bash

set -e

FreePSI=/disk1/wenqing/FreePSI-0.3/FreePSI/bin

READS=/disk1/wenqing/tmp_data/PFC_s2/data/GW08/
PSI_output=/disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/exon_psi/GW08

######所有染色体放进freepsi会导致内存崩溃，故分成单条染色体
num=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
for i in ${num[@]};
do 
    #cat /local/wenqing/data/psi_materials/annotation/hg19_refGene_exonBoundary.bed |awk '{if($1=="\"$i\"") print $0}' > /local/wenqing/data/psi_materials/annotation/"$i"_anno.bed
    #mkdir "$i"_genome
    #samtools faidx /local/wenqing/data/psi_materials/genome/hg19.fa "$i" > /local/wenqing/data/psi_materials/"$i"_genome/"$i".fa

    
    GENOME_DIR=/disk1/wenqing/tmp_data/freepsi/psi_materials/"$i"_genome
    BND_FILE=/disk1/wenqing/tmp_data/freepsi/psi_materials/annotation/"$i"_anno.bed

    K=27
    THREAD=20
    set -x

    ${FreePSI}/freePSI build\
        -k $K -p ${THREAD} \
        -g ${GENOME_DIR} \
        -a ${BND_FILE} \
        -r ${READS}/reads.1.fa \
        -o ${PSI_output}/hashtable.json

    ${FreePSI}/freePSI quant\
        -k $K -p ${THREAD} \
        -i ${PSI_output}/hashtable.json \
        -o ${PSI_output}/

    /disk1/wenqing/anaconda3/envs/py3.9_R4/bin/python ${FreePSI}/postProc.py \
        ${PSI_output}/psi_freePSI_raw.json \
        ${PSI_output}/psi_freePSI.json

    /disk1/wenqing/anaconda3/envs/py3.9_R4/bin/python ${FreePSI}/summary.py \
        ${BND_FILE} \
        ${PSI_output}/psi_freePSI.json \
        ${PSI_output}/psi_freePSI.summary_"$i"
done