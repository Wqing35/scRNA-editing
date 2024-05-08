#!/bin/bash 

##合并所有bam
sample=(ctr_male_pfc/ctr_32 ctr_male_pfc/ctr_56 asd_male_pfc/asd_17 asd_male_pfc/asd_18 asd_male_pfc/asd_57)
for one in ${sample[@]};
do 

    cd /disk1/wenqing/tmp_data/ASD/"$one"/tmp/genome
    cd ../genome
    samtools sort all.bam all.srt
    cd ../genome_mskAG
    samtools sort all.bam all.srt
    cd ../genome_mskTC
    samtools sort -@ 10 all.bam all.srt
    cd ../transcript
    samtools sort -@ 10 all.bam all.srt
    cd ../transcript_mskTC
    samtools sort -@ 10 all.bam all.srt
    cd ../transcript_mskAG
    samtools sort -@ 10 all.bam all.srt

    cd ../
    samtools merge -@ 20 ./genome.all.bam genome/all.srt.bam genome_mskTC/all.srt.bam genome_mskAG/all.srt.bam 
    samtools merge -@ 20 ./transcript.all.bam transcript/all.srt.bam transcript_mskTC/all.srt.bam transcript_mskAG/all.srt.bam

    samtools sort -@ 10 ./genome.all.bam ./genome.all.srt
    samtools sort -@ 10 ./transcript.all.bam ./transcript.all.srt

    samtools index ./genome.all.srt.bam
    samtools index ./transcript.all.srt.bam
done

