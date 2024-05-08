#!/bin/bash

c=0
for line in `cat /disk1/wenqing/tmp_data/gene_region.bed`
do
  one_region=$line
  all_regions[$c]=$one_region
  ((c++))
done

a=0
for line in `cat /disk1/wenqing/tmp_data/gene_name.bed`
do
  one_name=$line
  all_names[$a]=$one_name
  ((a++))
done

for i in `seq 0 55764`
do
    region=${all_regions[$i]}
    gene=${all_names[$i]}
    
    echo $gene
    #bulk
    samtools view -b /disk1/wenqing/tmp_data/human_brain/brain_19/tmp/genome/all.bam "$region" > /disk1/wenqing/tmp_data/coverage_compare/bam_files/bulk/"$gene"_bulk.bam
    bedtools genomecov -ibam /disk1/wenqing/tmp_data/coverage_compare/bam_files/bulk/"$gene"_bulk.bam -bga -split > /disk1/wenqing/tmp_data/coverage_compare/coverage_res/bulk/"$gene"_depth_clean_bulk.bed
    echo 'done'
    #10X
    samtools view -b /disk1/wenqing/tmp_data/coverage_compare/regular_all.bam "$region" > /disk1/wenqing/tmp_data/coverage_compare/bam_files/10X/"$gene"_10X.bam
    bedtools genomecov -ibam /disk1/wenqing/tmp_data/coverage_compare/bam_files/10X/"$gene"_10X.bam -bga -split > /disk1/wenqing/tmp_data/coverage_compare/coverage_res/10X/"$gene"_depth_clean_10X.bed
    #smart-seq
    #samtools view -b /disk1/wenqing/tmp_data/coverage_compare/smart_all.bam "$region" > /disk1/wenqing/tmp_data/coverage_compare/bam_files/smart_seq/"$gene"_smart.bam
    #bedtools genomecov -ibam /disk1/wenqing/tmp_data/coverage_compare/bam_files/smart_seq/"$gene"_smart.bam -bga -split > /disk1/wenqing/tmp_data/coverage_compare/coverage_res/smart_seq/"$gene"_depth_clean_smart.bed
done