#!/usr/bin/python2
# coding=utf-8

import os  
  
def extract_reads(barcode_file, read1_file, read2_file, output_read1, output_read2):  
    barcodes = set(open(barcode_file).read().strip().split("\n"))  
    read2_handle = open(read2_file, "r")  
    read1_handle = open(read1_file, "r")  
    output_read1_handle = open(output_read1, "w")  
    output_read2_handle = open(output_read2, "w")  
  
    while True:  
        read2_header = read2_handle.readline()  
        read2_seq = read2_handle.readline()  
        read2_plus = read2_handle.readline()  
        read2_qual = read2_handle.readline() 

        read1_header = read1_handle.readline()  
        read1_seq = read1_handle.readline()  
        read1_plus = read1_handle.readline()  
        read1_qual = read1_handle.readline() 
  
        if not read2_header:  
            break  # 到达 read2 文件的末尾  
  
        #read2_seq = read2_seq.strip()  
        if read2_seq[:8] in barcodes:  
            # 提取 read1 的对应读段    
            output_read1_handle.write(read1_header+'_'+read2_seq[:8])  
            output_read1_handle.write(read1_seq)  
            output_read1_handle.write(read1_plus)  
            output_read1_handle.write(read1_qual)  
  
            # 写入 read2 的读段  
            output_read2_handle.write(read2_header)  
            output_read2_handle.write(read2_seq[17:])  
            output_read2_handle.write(read2_plus)  
            output_read2_handle.write(read2_qual[17:])  
  
    read2_handle.close()  
    read1_handle.close()  
    output_read1_handle.close()  
    output_read2_handle.close()  
  
# 使用示例
celltypes = ['GABAergic_neurons','Stem_cell','Neurons','OPC']

#GW16_sample=['SRR6367165','SRR6367166','SRR6367171']
GW16_sample=['SRR6367166','SRR6367171']
j = 0
suffix = ['4','9']
for one in GW16_sample:
    read1_file = "/disk1/wenqing/tmp_data/PFC_s2/all_fastq/"+one+"_1.fastq"  
    read2_file = "/disk1/wenqing/tmp_data/PFC_s2/all_fastq/"+one+"_2.fastq"  
    celltypes = ['GABAergic_neurons','Stem_cell','Neurons','OPC']
    i = suffix[j]
    for celltype in celltypes:
        barcode_file = "/disk1/wenqing/tmp_data/PFC_s2/20240425_use/GW16_1_"+str(i)+"/"+celltype+"_barcode.txt"
        output_read1 = "/disk1/wenqing/tmp_data/PFC_s2/data/GW16_1_"+str(i)+"/"+celltype+"_read1.fastq"  
        output_read2 = "/disk1/wenqing/tmp_data/PFC_s2/data/GW16_1_"+str(i)+"/"+celltype+"_read2.fastq"  
        #print(barcode_file+'\n'+read1_file+'\n'+read2_file+'\n'+output_read1+'\n'+output_read2)
        extract_reads(barcode_file, read1_file, read2_file, output_read1, output_read2)
    j+=1