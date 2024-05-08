#!/bin/bash

#stringtie -e -A asd_18_gene_abund.tab /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18/tmp/genome/all.bam -G /disk1/wenqing/tmp_data/gencode.v19.annotation.gtf -o ./asd_18.gtf -p 20 -B
#stringtie -e -A asd_57_gene_abund.tab /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/tmp/genome/all.bam -G /disk1/wenqing/tmp_data/gencode.v19.annotation.gtf -o ./asd_57.gtf -p 20 -B

#stringtie -e -A ctr_20_gene_abund.tab /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20/tmp/genome/all.bam -G /disk1/wenqing/tmp_data/gencode.v19.annotation.gtf -o ./ctr_20.gtf -p 20 -B
#stringtie -e -A ctr_23_gene_abund.tab /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_23/tmp/genome/all.bam -G /disk1/wenqing/tmp_data/gencode.v19.annotation.gtf -o ./ctr_23.gtf -p 20 -B
#stringtie -e -A ctr_56_gene_abund.tab /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56/tmp/genome/all.bam -G /disk1/wenqing/tmp_data/gencode.v19.annotation.gtf -o ./ctr_56.gtf -p 20 -B


#bedtools intersect -wo -a /disk1/wenqing/tmp_data/gene.bed -b resInAluWithADDP_over1|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$9}' > geneWithRegular_res.txt
#sed -i 's/;//g' geneWithRegular_res.txt
#sed -i 's/\"//g' geneWithRegular_res.txt

#从fq开始
#gffread Homo_sapiens.GRCh37.87.gtf -g Homo_sapiens.GRCh37.dna.primary_assembly.fa -w GRCh37.transcript.fa.tmp
#cut -f 1 -d ' ' GRCh37.transcript.fa.tmp > GRCh37.transcript.fa
#建立索引
#salmon index -t /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.fa -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -k 31
#asd样本
#salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/SRR9262917_2.fastq -o /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17_transcripts_quant_gencode -p 10
#salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/SRR9262957_2.fastq -o /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57_transcripts_quant_gencode -p 2
#salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/SRR9262918_2.fastq -o /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18_transcripts_quant_gencode -p 2
#ctr样本
#salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/SRR9262920_2.fastq -o /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20_transcripts_quant_gencode -p 5
salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/SRR9262923_2.fastq -o /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_23_transcripts_quant_gencode -p 5
salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/SRR9262956_2.fastq -o /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56_transcripts_quant_gencode -p 5




