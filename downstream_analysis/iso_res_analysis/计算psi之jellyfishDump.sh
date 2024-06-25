#!/bin/bash

set -e


Jellyfish=/local/wenqing/jellyfish/bin

one='all'
#sample=(GW16)
sample=(GW16 GW18 GW20_1 GW22_2 GW25)
for phase in ${sample[@]};
do
    READS=/local/wenqing/data/sc/human_brain/hippo/fqs_"$one"_exon_Ver2/"$phase"_"$one"_exon_fqs
    cd /local/wenqing/data/sc/human_brain/hippo/fqs_"$one"_exon_Ver2/"$phase"_"$one"_exon_fqs
    i=$(ls -l |grep reads.1.jf_|wc -l)
    nums=$(($i-1))
    num=$(seq 0 $nums)
    for number in ${num[@]};
    do
        ${Jellyfish}/jellyfish dump ${READS}/reads.1.jf_$number -o ${READS}/reads.1.fa_$number
    done

    cat ${READS}/reads.1.fa_* > ${READS}/reads.1.fa_all
done
