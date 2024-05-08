#!/bin/bash

celltypes=(GABAergic_neurons Neurons Stem_cell OPC)
for one in ${celltypes[@]};
do
        #GW12
	python2 /disk1/wenqing/scRNA-editing/sprint/sprint_main_hisat2.py -rp /disk1/wenqing/tmp_data/hg19_repeat.txt -c 0 -p 30 -1 /disk1/wenqing/tmp_data/PFC_s2/data/GW12/"$one"_read1.fastq -2 /disk1/wenqing/tmp_data/PFC_s2/data/GW12/"$one"_read2.fastq /disk1/wenqing/tmp_data/hg19/hg19.fa /disk1/wenqing/tmp_data/PFC_s2/result/GW12/"$one" /disk1/wenqing/hisat2-2.2.0/hisat2 /disk1/wenqing/SPRINT/SPRINT_master/samtools_and_bwa/samtools
done
