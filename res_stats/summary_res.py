#!/usr/bin/python2
# coding=utf-8

import subprocess,os,sys
#import tmp as sprint
import re

def main():
    def summaryRES(bed_in_path=0, txt_out_path=0):
        TYPE=['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG']
        print TYPE
        NUM=[0,0,0,0,0,0,0,0,0,0,0,0]
        print NUM

        fi=open(bed_in_path)
        fo=open(txt_out_path,'w')
        for line in fi:
            seq=line.rstrip().split('\t')
            this_type=seq[3]
            NUM[TYPE.index(this_type)]+=1
        
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
    """ summaryRES(bed_in_path="./tmp/hyper_mskAG_alu.res",txt_out_path="./hyperAG_alu_res.txt") 
    summaryRES(bed_in_path="./tmp/hyper_mskTC_alu.res",txt_out_path="./hyperTC_alu_res.txt")   """
    #summaryRES(bed_in_path="./tmp/regular.res",txt_out_path="./regular_res.txt")
    summaryRES(bed_in_path="/local/wenqing/data/sc/test_data_GW16_1/SPRINT_identified_all.res",txt_out_path="./all_res.txt") 



    

main()