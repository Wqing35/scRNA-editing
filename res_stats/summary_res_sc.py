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

    """ summaryRES(bed_in_path="./tmp/hyper_mskAG_alu.res",txt_out_path="./hyperAG_alu_res.txt") 
    summaryRES(bed_in_path="./tmp/hyper_mskTC_alu.res",txt_out_path="./hyperTC_alu_res.txt")   """
    sample=['20','32','56']
    for phase in sample:
        cell_types=["L2_3_ExN","L4","L5_6","L5_6_CC","IN_VIP","Olig"]
        for type in cell_types:
            res_in_dir='/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_'+phase+'/regular_merged.res_over1_'+type
            summaryRES(res_in_dir,res_in_dir+'.txt') 
    """ summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18/SPRINT_identified_all.res",txt_out_path="/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18/all.res.txt") 
    summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20/SPRINT_identified_all.res",txt_out_path="/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20/all.res.txt") 
    summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_23/SPRINT_identified_all.res",txt_out_path="/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_23/all.res.txt") 
    summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56/SPRINT_identified_all.res",txt_out_path="/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56/all.res.txt")  """

    
    """ summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/SPRINT_identified_all.res",txt_out_path="./all.res.txt") 
    summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/SPRINT_identified_regular.res",txt_out_path="./regular.res.txt") 
    summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/SPRINT_identified_all.res",txt_out_path="../asd_57/all.res.txt") 
    summaryRES(bed_in_path="/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/SPRINT_identified_regular.res",txt_out_path="../asd_57/regular.res.txt")  """

main()