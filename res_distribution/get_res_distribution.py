#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import pandas as pd  
import matplotlib.pyplot as plt 
#mystat是一个外部定义的函数
import mystat
import re
import string
from optparse import OptionParser
import warnings
import collections
import math
from time import strftime
import subprocess
from os.path import basename
import operator
from __future__ import print_function  
import sys  
  

def main():
    
    def printlog (mesg):
        '''print progress into stderr and log file'''
        mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
        LOG=open('log.txt','a')
        print(mesg, file=LOG)
        print(mesg, file=sys.stderr)

    def genebody_percentile(refbed, mRNA_len_cut = 100):
        '''
        return percentile points of gene body
        mRNA length < mRNA_len_cut will be skipped
        '''
        
        g_percentiles = {}
        transcript_count = 0
        for line in open(refbed,'r'):
            try:
                if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5]
                geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])
                    
                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                transcript_count += 1
            except:
                print("[NOTE:input bed must be 12-column] skipped this line")
                continue
            gene_all_base=[]
            mRNA_len =0
            flag=0
            for st,end in zip(exon_starts,exon_ends):
                ###只记录了exon的长度，导致下游计算res的分布时，res的数量会少很多（res主要分布在intron上）
                gene_all_base.extend(list(range(st+1,end+1)))		#1-based coordinates on genome
            if len(gene_all_base) < mRNA_len_cut:
                continue
            g_percentiles[geneID] = (chrom, strand, mystat.percentile_list (gene_all_base))	#get 100 points from each gene's coordinates
        printlog("Total " + str(transcript_count) + ' transcripts loaded')
        return g_percentiles



 
  
    # 假设edit_sites.bed是RNA编辑位点的BED文件  
    edit_sites_file = '/disk1/wenqing/tmp_data/pbmc/result/s2/tmp/regular.res.depth'  
    # 假设gene_percentiles是genebody_percentile脚本输出的字典  
    gene_percentiles = genebody_percentile(refbed = '/disk1/wenqing/tmp_data/transcript_NM_022455.bed', mRNA_len_cut = 100) 
    
    # 读取RNA编辑位点的BED文件  
    edit_sites_df = pd.read_csv(edit_sites_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'type', 'Supporting_reads', 'strand', 'info'])  
    
    # 创建一个字典来存储每个基因区域的编辑位点数量  
    edit_counts_per_percentile={}
    for i in range(1,101):
        edit_counts_per_percentile[str(i)]=0
    
    # 遍历每个编辑位点  
    for index, row in edit_sites_df.iterrows():  
        chrom = row['chrom']  
        end = int(row['end'])  
        strand = row['strand']  
        #mid_point = (start + end) / 2.0  
        #print(chrom+'_'+str(end)+' running:')
        
        # 遍历每个基因区域的划分结果  
        for geneID, (gene_chrom, gene_strand, percentiles) in gene_percentiles.items():  
            if chrom == gene_chrom and strand == gene_strand:  
                gene_length = percentiles[-1]  
                if end < gene_length:
                    for i in range(1,101): 
                        if i==1:
                            if  end < percentiles[i] and end > (percentiles[i]-100):
                                edit_counts_per_percentile[str(i)]+=1
                                break
                        else:
                            if end > percentiles[i-1] and end < percentiles[i]:
                                edit_counts_per_percentile[str(i)]+=1
                                break
                    # 检查编辑位点的中点是否在当前基因区域的百分比位点范围内  
                    """ if abs((mid_point - percentile)/ gene_length) < 0.01:  
                        edit_counts_per_percentile[percentile] += 1  
                        break   """
    
    # 将编辑位点计数转换为列表，以用于绘图  
    int_percentiles = sorted([int(x) for x in edit_counts_per_percentile.keys()])
    str_percentiles = [str(x) for x in int_percentiles]
    counts = [edit_counts_per_percentile[p] for p in str_percentiles]  
    
    # 创建一个数据字典  
    data = {  
        'percentile': int_percentiles,  
        'counts': counts  
    }  

    """ data={
        'res_pos':list(notIn_exon_res),
        'ord':range(1,len(notIn_exon_res)+1)
    } """
    
    # 创建数据框  
    df = pd.DataFrame(data)  
    # 将数据框写入CSV文件  
    df.to_csv('~/outout_s2.csv', index=False)


main()