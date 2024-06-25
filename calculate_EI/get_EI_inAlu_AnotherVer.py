#!/usr/bin/python2
# coding=utf-8

import subprocess,os,sys
import re
import numpy as np
from collections import defaultdict


#####以Alu region为行的计算方式

""" 
./RNAEditingIndex -d /local/wenqing/tmp_data/GW16 -f all_regular.bam -l ~/tmp_data/GW16/ -o ~/tmp_data/GW16/ -os ~/tmp_data/GW16/ --genome hg19

#将多行序列转换为单行
awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\n":$0}' ALUs.fa > Alu_singleLines.fa 
#计数A碱基
awk '{if($0~/^>/)a=$0;else {gsub(/[tcgn]/,"",$0);printf a"\t"length($0)"\n"}}' Alu_singleLines.fa > AsInALU.txt


input='AsInALU.txt'
output='resInAlu.txt'
fi=open(input)
fo=open(output,'w')

i=1
for line in fi:
    seq=line.rstrip().split('\t')
    region=seq[0]
    chr=region.split(':')[0]
    start=region.split(':')[1].split('-')[0]
    end=region.split(':')[1].split('-')[1]
    name='Alu'+str(i)
    i=i+1
    fo.write(chr+'\t'+start+'\t'+end+'\t'+str(seq[1])+'_'+name+'\n')
fo.close()
fi.close()



sed -i 's/>//g' AsInAlu.txt 


cat regular.res.depth.filterLowQualCell.anno|grep Alu |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > regular_merged.res.alu.modified

bedtools intersect -wb -a regular_merged.res.alu.modified -b ~/tmp_data/Alu_seq/AsInAlus_modified.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$11"_"$12}' > resInAluWithADDP

##筛选coverge reads1/2/3
cat regular_merged.res |awk '{print $7}'|awk -F ':' '{print $2}' > dp.txt
paste -d '\t' regular_merged.res dp.txt > regular_merged.res_WizDepth
cat regular_merged.res_WizDepth|awk '{if ($8 > 1) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > regular.res_over1

bedtools intersect -wb -a resInAlu1.txt -b ~/tmp_data/GW18/GW18_depth_clean_regular.bed | awk '{print $6=$7=$8=""; print $0}' > resInAluWithDepth.txt 

############step2############
#合并regular的reads得到regular_all.zz
cat 





"""

def main():

    #def get_AEI(bed_in_dir,bed_anno_dir,bed_out_dir,res_in_dir,AEI_out_dir)
    def get_AEI(element,depth_index):
        res_in_dir='resInAluWithADDP_over'+str(depth_index)

        fi_reads=open('./'+element+'_combined.zz.sorted.modified_InAlus_trimmed_reads_OnlyAG')
        #fi_reads=pd.read_table(bed_in_dir_reads,sep='\t',header=0)

        alu_dict1={}
        for line in fi_reads:
            seq=line.rstrip().split('\t')
            alu_region=seq[6]            
            try:
                alu_dict1[alu_region].append([seq[0],seq[1],seq[2],seq[7]])
            except Exception, e:
                alu_dict1[alu_region]=[[seq[0],seq[1],seq[2],seq[7]]]        
        fi_reads.close()
        
        fi_res=open(res_in_dir)
        alu_dict={}
        for line in fi_res:
            seq=line.rstrip().split('\t')
            alu_region=seq[6]   
            ad=seq[4].split(':')[0]         
            try:
                alu_dict[alu_region].append([seq[0],seq[1],seq[2],seq[3],ad])
            except Exception, e:
                alu_dict[alu_region]=[[seq[0],seq[1],seq[2],seq[3],ad]]        
        fi_res.close()

        #AEIByBarcode={}

        AGByAlu = defaultdict(list)
        AsByAlu = defaultdict(list)

        for key in alu_dict.keys():
            alu_region = key
            lengthOfA = sum(int(record[3]) for record in alu_dict1.get(key, []))
            
            AG_num = sum(int(record[4]) for record in alu_dict.get(key, []) if record[3] == 'AG')
            
            if lengthOfA > 0:
                AGByAlu[alu_region].append([AG_num])
                AsByAlu[alu_region].append([lengthOfA])

        AGByAlu = dict(AGByAlu)
        AsByAlu = dict(AsByAlu)
        
        sumOfAGByAlu={}
        sumOfAsByAlu={}
        for key in AGByAlu.keys():
            AG_num_set=AGByAlu[key]
            As_num_set=AsByAlu[key]
            #AverAEI=AEIByBarcode[key]
            tmp_AG=[]
            tmp_As=[]
            for ele in AG_num_set:
                if isinstance(ele,int):
                    tmp_AG.append(int(ele))
                else:
                    tmp_AG.append(int(ele[0]))
            for ele in As_num_set:
                if isinstance(ele,int):
                    tmp_As.append(int(ele))
                else:
                    tmp_As.append(int(ele[0]))            
            sumOfAGByAlu[key]=np.sum(np.array(tmp_AG))
            sumOfAsByAlu[key]=np.sum(np.array(tmp_As))
        
        all_AGs=sum(sumOfAGByAlu.values())
        non_zero_indices_AG = [index for index, value in enumerate(sumOfAGByAlu.values()) if value != 0]
        all_As=sum([sumOfAsByAlu.values()[i] for i in non_zero_indices_AG])

        AEI=float(all_AGs)/all_As

        print AEI


    print 'start'
    #dp_idx=[0,1,2,3]
    dp_idx=[0]
    for idx in dp_idx:
        get_AEI('regular',idx)
    print 'done'

main() 