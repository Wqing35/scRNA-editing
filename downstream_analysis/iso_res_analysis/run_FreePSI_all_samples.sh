#!/bin/bash

###########计算mRNA/lncRNA exon的psi
#step1:安装jellyfish并准备jellyfish所需的文件
#wget -c https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz
#tar -zxvf jellyfish-2.3.1.tar.gz
#mkdir jellyfishlocation
#cd jellyfish-2.3.1
#./configure  --prefix=/jellyfishlocation
#make  -j 4
#make install


#cd Eigen
#mkdir build
#cd build
#mkdir /disk1/wenqing/eigen_lib
#cmake .. -DCMAKE_INSTALL_PREFIX=/disk1/wenqing/eigen_lib
#make install
#export CPATH=${CPATH}:/disk1/wenqing/eigen_lib/include/eigen3

#step2:安装freepsi,准备freepsi所需要的文件GENOME_DIR、BND_FILE、READS
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/disk1/wenqing/anaconda3/lib/

set -e

# Provide the directory containing Jellyfish (usually named as 'bin')  
Jellyfish=/disk1/wenqing/jellyfish-2.3.1/bin
# Provide the directory containing FreePSI
# E.g. 
#FreePSI=../bin
FreePSI=/disk1/wenqing/FreePSI-0.3/FreePSI/bin

K=27
THREAD=20
samples=(GW12  GW16_1_3  GW16_1_4  GW16_1_9  GW19_1_1  GW19_1_2  GW19_1_3  GW23_1_1  GW23_1_2  GW23_1_3  GW26_1_1)
for one in ${samples[@]};
do
    READS=/disk1/wenqing/tmp_data/PFC_s2/data/"$one"/
    PSI_output=/disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/exon_psi/"$one"/


    set -x
    # Count k-mers in RNA-seq reads using jellyfish
    ${Jellyfish}/jellyfish count -m ${K} -s 100M -t ${THREAD} ${READS}/GABAergic_neurons_read1.fastq ${READS}/GABAergic_neurons_read2.fastq -o ${READS}/reads.1.jf


    ${Jellyfish}/jellyfish dump ${READS}/reads.1.jf -o ${READS}/reads.1.fa

    num=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
    for i in ${num[@]};
    do 
        GENOME_DIR=/disk1/wenqing/tmp_data/freepsi/psi_materials/"$i"_genome
        BND_FILE=/disk1/wenqing/tmp_data/freepsi/psi_materials/annotation/"$i"_anno.bed

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

    cat ${PSI_output}/psi_freePSI.summary_* > ${PSI_output}/psi_freePSI.summary_allChr

done
