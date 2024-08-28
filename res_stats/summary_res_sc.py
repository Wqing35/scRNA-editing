#!/usr/bin/python2
# coding=utf-8

import subprocess,os,sys
#import tmp as sprint
import re

def main():
    #summaryRES检测的是RNAediting的位点、而非具有res的reads数量
    #因此，单细胞的summary数量需要另外修改代码
    def summaryRES(bed_in_path=0, txt_out_path=0):
        TYPE=['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG']
        print TYPE
        NUM=[0,0,0,0,0,0,0,0,0,0,0,0]
        print NUM

        fi=open(bed_in_path)
        fo=open(txt_out_path,'w')
        old=set()
        for line in fi:
            seq=line.rstrip().split('\t')
            chr=seq[0]
            this_type=seq[3]
            site=seq[2]
            if this_type+':'+chr+':'+site not in old:
                NUM[TYPE.index(this_type)]+=1
                old.add(this_type+':'+chr+':'+site)
        
        SUM=sum(NUM)
        header='#type\tnumber\tpercentage\n'
        i=0
        while i <len(TYPE):
            this_type=TYPE[i]
            this_num=str(NUM[i])
            this_ratio=str(round(float(NUM[i])/float(SUM)*100,2))
            fo.write(this_type+'\t'+this_num+'\t'+this_ratio+'\n')
            i=i+1

        #####################
        this_type='AG+TC'
        this_num =str(NUM[TYPE.index('AG')]+NUM[TYPE.index('TC')])
        this_ratio=str(round(float(this_num)/float(SUM)*100 ,2))
        fo.write(this_type+'\t'+this_num+'\t'+this_ratio+'\n')
        ###################
        
        #####################
        this_type='Total'
        this_num =str(SUM)
        this_ratio=str(round(float(this_num)/float(SUM)*100 ,2))
        fo.write(this_type+'\t'+this_num+'\t'+this_ratio+'\n')
        ###################
    
        fo.close()


    summaryRES(bed_in_path="/disk1/wenqing/tmp_data/PFC_s2/result/GW16_1_4_onecell/tmp/regular.res.depth",txt_out_path="./regular.res.txt")  

main()