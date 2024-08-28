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
    
    #samples=["GW16_1_3","GW16_1_4","GW19_1_1","GW19_1_3","GW23_1_1","GW23_1_2","GW23_1_3","GW26_1_1","GW19_1_2"]
    """ celltypes=['GABA','neuron','opc']
    for ad in [7,8,9,10]:
        for dp in [11,12,13,14,15,16,17,18,19,20]:
            for celltype in celltypes:
                bed_in_path='/disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/res_stats/'+str(ad)+'_'+str(dp)+'_'+celltype+'_all.res.depth'
                txt_out_path='/disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/res_stats/result/'+str(ad)+'_'+str(dp)+'_'+celltype+'_all.res.txt'
                if os.path.exists(bed_in_path)==False:
                    next
                else:
                    summaryRES(bed_in_path,txt_out_path)  """
    """ sample=['GW16_1_4','GW19_1_1','GW23_1_2','GW26_1_1']
    celltypes=['GABAergic_neurons','Neurons','OPC']
    for one_sample in sample:
        for celltype in celltypes:
            bed_in_path='/disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/'+one_sample+'/'+celltype+'/regular.res.depth'
            txt_out_path='/disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/'+one_sample+'/'+celltype+'/regular.res.txt'
            summaryRES(bed_in_path,txt_out_path) """
    
    #for idx in range(1,51):
    #    bed_in_path='/disk1/wenqing/tmp_data/PFC_s2/test_filter_in_snv_neuron/filter_idx'+str(idx)+'/regular.res.depth'
    #    txt_out_path='/disk1/wenqing/tmp_data/PFC_s2/test_filter_in_snv_neuron/filter_idx'+str(idx)+'/regular.res.txt'
    #    summaryRES(bed_in_path,txt_out_path)
    sample=['GW08','GW12','GW16_1_3','GW16_1_9','GW19_1_1','GW19_1_2','GW19_1_3','GW23_1_1','GW23_1_2','GW23_1_3','GW26_1_1']
    for one_sample in sample:
        bed_in_path_readInfo='/disk1/wenqing/tmp_data/PFC_s2/result/'+one_sample+'/GABAergic_neurons/tmp/readInfo_dpFiltered/tmp/regular.res.depth'
        txt_out_path_readInfo='/disk1/wenqing/tmp_data/PFC_s2/result/'+one_sample+'/GABAergic_neurons/tmp/readInfo_dpFiltered/tmp/regular.res.txt'
        summaryRES(bed_in_path_readInfo,txt_out_path_readInfo)
        
        bed_in_path_bulk='/disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/'+one_sample+'/GABAergic_neurons/regular.res.depth'
        txt_out_path_bulk='/disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/'+one_sample+'/GABAergic_neurons/regular.res.txt'        
        summaryRES(bed_in_path_bulk,txt_out_path_bulk)

main()