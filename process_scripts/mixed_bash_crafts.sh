bedtools intersect -wo -a ~/tmp_data/gene.bed -b ./SPRINT_identified_all.res |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$11"\t"$12}' > geneWithRes_all.txt

sed -i 's/;//g' geneWithRes_all.txt
sed -i 's/\"//g' geneWithRes_all.txt


samtools view -hb chr1:145507597-145513536 



input_dir='/disk1/wenqing/tmp_data/ASD/celltype_fq'
ref='/disk1/wenqing/tmp_data/hg19/hg19.fa'
output_dir='/disk1/wenqing/tmp_data/ASD/celltype_bam'
phase=(18 57 20 32 56)
for one in ${phase[@]};
do
    echo "$phase"

    bwa aln -t 20 "$ref" "$input_dir"/"$one"_L2_3_ExN_2.fastq > "$output_dir"/"$one"_L2_3_ExN_2.sai

    bwa samse -n4 "$ref" "$output_dir"/"$one"_L2_3_ExN_2.sai "$input_dir"/"$one"_L2_3_ExN_2.fastq > "$output_dir"/"$one"_L2_3_ExN_2.sam

    samtools view -bS "$output_dir"/"$one"_L2_3_ExN_2.sam > "$output_dir"/"$one"_L2_3_ExN_2.bam

    samtools sort "$output_dir"/"$one"_L2_3_ExN_2.bam "$output_dir"/"$one"_L2_3_ExN_2.srt

done


samtools index 17_L2_3_ExN_2.srt.bam.bam
samtools view -hb 17_L2_3_ExN_2.srt.bam.bam chr1:145507597-145513536 > L2_3_RBM8A_17.bam

samtools index 18_L2_3_ExN_2.srt.bam.bam
samtools view -hb 18_L2_3_ExN_2.srt.bam.bam chr1:145507597-145513536 > L2_3_RBM8A_18.bam

samtools index 57_L2_3_ExN_2.srt.bam.bam
samtools view -hb 57_L2_3_ExN_2.srt.bam.bam chr1:145507597-145513536 > L2_3_RBM8A_57.bam

samtools index 20_L2_3_ExN_2.srt.bam.bam
samtools view -hb 20_L2_3_ExN_2.srt.bam.bam chr1:145507597-145513536 > L2_3_RBM8A_20.bam

samtools index 32_L2_3_ExN_2.srt.bam.bam
samtools view -hb 32_L2_3_ExN_2.srt.bam.bam chr1:145507597-145513536 > L2_3_RBM8A_32.bam

samtools index 56_L2_3_ExN_2.srt.bam.bam
samtools view -hb 56_L2_3_ExN_2.srt.bam.bam chr1:145507597-145513536 > L2_3_RBM8A_56.bam

samtools index L2_3_RBM8A_17.bam
samtools index L2_3_RBM8A_18.bam
samtools index L2_3_RBM8A_57.bam
samtools index L2_3_RBM8A_20.bam
samtools index L2_3_RBM8A_32.bam
samtools index L2_3_RBM8A_56.bam


bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/asd_iso_IR_exon_region_exd.txt -b /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_all.res|more

bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/asd_iso_no_IR_exon_region_exd.txt -b /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_all.res|awk '{print $18}' > asd_res_in_iso_no_IR_region_exd.txt


bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/ctr_iso_no_IR_exon_region_exd.txt -b /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_all.res|awk '{print $18}' > ctr_res_in_iso_no_IR_region_exd.txt

bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/ctr_iso_IR_exon_region_exd.txt -b /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_all.res|awk '{print $18}' > ctr_res_in_iso_IR_region_exd.txt



bedtools intersect -wo -a ~/tmp_data/gene.bed -b tmp/all.res |grep GRIA2./ |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$11"\t"$12}' > geneWithRes_all.txt




cat gencode.v19.annotation.gtf|grep FLNA > ./FLNA.gencode.v19.anno.gtf
cat gencode.v19.annotation.gtf|grep RBFOX1 > ./RBFOX1.gencode.v19.anno.gtf
#FLNA exon区域
cat FLNA.gencode.v19.anno.gtf |
awk 'BEGIN{OFS="\t";} {if ($3=="exon" || $3=="UTR") print $1,$4-1,$5}' |
bedtools sort |
bedtools merge -i - | gzip > FLNA_exon.bed.gz


cat FLNA.gencode.v19.anno.gtf |
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
bedtools sort |
bedtools subtract -a stdin -b FLNA_exon.bed.gz |
gzip > FLNA_intron.bed.gz

#RBFOX1
cat RBFOX1.gencode.v19.anno.gtf |
awk 'BEGIN{OFS="\t";} {if ($3=="exon" || $3=="UTR") print $1,$4-1,$5}' |
bedtools sort |
bedtools merge -i - | gzip > RBFOX1_exon.bed.gz


cat RBFOX1.gencode.v19.anno.gtf |
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
bedtools sort |
bedtools subtract -a stdin -b RBFOX1_exon.bed.gz |
gzip > RBFOX1_intron.bed.gz


#s2
bedtools intersect -wo -a ~/tmp_data/gene.bed -b ./tmp/all.res.depth | grep FLNA|awk '{print $7}' > FLNA.txt
#FLNA的res全部分布在exon
bedtools intersect -wo -a ~/tmp_data/FLNA_exon.bed -b ./tmp/all.res.depth > FLNA_exon.res

bedtools intersect -wo -a ~/tmp_data/gene.bed -b ./tmp/all.res.depth | grep RBFOX1|awk '{print $7}' > RBFOX1.txt
#RBFOX1的res全部分布在intron
bedtools intersect -wo -a ~/tmp_data/RBFOX1_intron.bed -b ./tmp/all.res.depth > RBFOX1_intron.txt



#10X
bedtools intersect -wo -a ~/tmp_data/gene.bed -b ./tmp/all.res.depth | grep FLNA|awk '{print $7}' > FLNA.txt
#FLNA的res全部分布在exon
bedtools intersect -wo -a ~/tmp_data/FLNA_exon.bed -b ./tmp/all.res.depth > FLNA_exon.res

bedtools intersect -wo -a ~/tmp_data/gene.bed -b ./tmp/all.res.depth | grep RBFOX1|awk '{print $7}' > RBFOX1.txt
#RBFOX1的res全部分布在intron
bedtools intersect -wo -a ~/tmp_data/RBFOX1_intron.bed -b ./tmp/all.res.depth > RBFOX1_intron.txt


#######计算每个样本3种细胞类型发生de且IR的转录本RNA编辑的水平
bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/deX_isoform_IR/de_isoform_IR_L2_3_ExN.newExd.txt -b /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/L2_3_ExN.res > /disk1/wenqing/tmp_data/ASD/celltype_deX_iso_IR_res/17/L2_3_de_iso_IR.res

bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/deX_isoform_IR/de_isoform_IR_L2_3_ExN.newExd.txt -b /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18/L2_3_ExN.res > /disk1/wenqing/tmp_data/ASD/celltype_deX_iso_IR_res/18/L2_3_de_iso_IR.res

bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/deX_isoform_IR/de_isoform_IR_L2_3_ExN.newExd.txt -b /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/L2_3_ExN.res > /disk1/wenqing/tmp_data/ASD/celltype_deX_iso_IR_res/57/L2_3_de_iso_IR.res

bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/deX_isoform_IR/de_isoform_IR_L2_3_ExN.newExd.txt -b /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20/L2_3_ExN.res > /disk1/wenqing/tmp_data/ASD/celltype_deX_iso_IR_res/20/L2_3_de_iso_IR.res

bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/deX_isoform_IR/de_isoform_IR_L2_3_ExN.newExd.txt -b /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_32/L2_3_ExN.res > /disk1/wenqing/tmp_data/ASD/celltype_deX_iso_IR_res/32/L2_3_de_iso_IR.res

bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/deX_isoform_IR/de_isoform_IR_L2_3_ExN.newExd.txt -b /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56/L2_3_ExN.res > /disk1/wenqing/tmp_data/ASD/celltype_deX_iso_IR_res/56/L2_3_de_iso_IR.res


seqkit grep -n -r -p "ENST00000308275.8_7" isoformSwitchAnalyzeR_isoform_nt.fasta > test.fasta
seqkit seq -n test.fasta > base_count.txt



bedtools intersect -wo -a /disk1/wenqing/tmp_data/ASD/deX_isoform_IR/de_isoform_IR_L2_3_ExN.newExd.txt -b /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/L2_3_ExN.res > /disk1/wenqing/tmp_data/ASD/celltype_deX_iso_IR_res/17/L2_3_de_iso_IR.res


git clone https://github.com/lh3/seqtk.git;
cd seqtk; make


wget https://github.com/lh3/seqtk/archive/v1.2.tar.gz
tar zxvf v1.2.tar.gz
$ cd seqtk-1.2
$ make


python2 /disk1/wenqing/SPRINT/SPRINT_master/sprint/sprint_main.py -rp /disk1/wenqing/tmp_data/hg19_repeat.txt -c 0 -p 30 -1 /disk1/wenqing/tmp_data/pbmc/10X_dat/10X_same_bpNumAsS2_1.fastq -2 /disk1/wenqing/tmp_data/pbmc/10X_dat/10X_same_bpNumAsS2_2.fastq /disk1/wenqing/tmp_data/hg19/hg19.fa /disk1/wenqing/tmp_data/pbmc/result/10X_Ver2 /disk1/wenqing/SPRINT/SPRINT_master/samtools_and_bwa/bwa /disk1/wenqing/SPRINT/SPRINT_master/samtools_and_bwa/samtools






#####从GNB1、CTSS两个基因来比较s2和10X的区别
#step1: regular.res & regular.snv数量的统计
bedtools intersect -wo -a ~/tmp_data/gene.bed -b regular.snv| awk '{if($4=="GNB1") print $0}' > GNB1.regular.snv
bedtools intersect -wo -a ~/tmp_data/gene.bed -b regular.res| awk '{if($4=="GNB1") print $0}' > GNB1.regular.res

bedtools intersect -wo -a ~/tmp_data/gene.bed -b regular.snv| awk '{if($4=="CTSS") print $0}' > CTSS.regular.snv
bedtools intersect -wo -a ~/tmp_data/gene.bed -b regular.res| awk '{if($4=="CTSS") print $0}' > CTSS.regular.res

#step2: 提取特定基因的bam文件计算平均深度
samtools view -hb -@ 8 -b ./all.bam chr1:150702671-150738433 > CTSS.bam
bamdeal statistics Coverage -i /disk1/wenqing/tmp_data/pbmc/result/10X_Ver3/tmp/genome/CTSS.bam -r /disk1/wenqing/tmp_data/hg19/hg19.fa -o ./CTSS_depth_100 -w 100

#step3: 系统的比较两平台regular.res在基因上的情况
bedtools intersect -wo -a ~/tmp_data/gene.bed -b ./SPRINT_identified_regular.res |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$11"\t"$12}' > geneWithRes_regular.txt
more geneWithRes_regular.txt |awk '{print $4}'|sort|uniq -c > 10X_geneWzRes_regular.txt


#######计算10X和S2基因编辑水平的相关性
sed -i 's/:/\t/g' regular.res.depth
awk '{if($8!=1) print $0}' regular.res.depth > regular.res.depth1
rm regular.res.depth
mv regular.res.depth1 regular.res.depth

bedtools intersect -wo -a ~/tmp_data/gene.bed -b ./regular.res.depth |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$8"\t"$11"\t"$12}' > geneWithRegular_res.txt


gtfToGenePred hg19.knownGene.gtf knownGene.genePred
genePredToBed knownGene.genePred knownGene.bed



    ####位于intron区的regular res:83%
    bedtools intersect -wo -a /disk1/wenqing/tmp_data/my_intron.srt.bed -b /disk1/wenqing/tmp_data/pbmc/result/10X_Ver3/tmp/regular.res.depth |awk '{print $4"_"$6}'|sort|uniq -c|wc -l
    #保留在intron区的res
    bedtools intersect -wo -a /disk1/wenqing/tmp_data/my_intron.srt.bed -b /disk1/wenqing/tmp_data/pbmc/result/10X_Ver3/tmp/regular.res.depth |awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > regular.res.depth.intron
    #去重
    awk '!seen[$1,$3,$4]++' regular.res.depth.intron > regular.res.depth.intron.uniq








