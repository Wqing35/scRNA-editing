#!/usr/bin/python2
# coding=utf-8

import subprocess,os
import re
import pandas as pd

def main():

    def get_EI(bed_in_dir,EI_out_dir):
        #######step1:获取res的anno以得到alu区域的res
        #转录本区域先不分alu区域，仅计算Editing Index（EI）  
        
            ################step1：将all/regular_all_combined.zz.sorted中比对得到一个区域的多个reads分行放置
            if os.path.exists('./ts_mapping.zz.sorted.modified')==False:
                print 'Prepare reads_file:'

                out=[]
                fi=open('ts_mapping.zz')
                for line in fi:
                    seq=line.split('\t')
                    out.append([seq[0],int(seq[3].split(':')[0]),int(seq[3].split(':')[-1]),line])
                fi.close()
                fo=open('ts_mapping.zz.sorted', 'w')
                out.sort()
                for one in out:
                    fo.write(one[3]) 

                fi.close()
                fo.close()

                fi=open('ts_mapping.zz.sorted')
                fo=open('ts_mapping.zz.sorted.modified','w')

                for line in fi:
                    seq=line.rstrip().split("\t")
                    ts_id=seq[0].split('|')[0]
                    regions=seq[3].split(";")
                    read_info=seq[8]
                    if len(regions) ==1:
                        regions=regions[0]
                        start=regions.split(":")[0]
                        end=regions.split(":")[1]
                        fo.write(ts_id+'\t'+start+'\t'+end+'\t'+seq[4]+'\t'+seq[5]+'\t'+seq[7]+'\t'+read_info+'\n')
                    else:
                        for i in range(0,len(regions)):
                            tmp_region=regions[i]
                            start=tmp_region.split(":")[0]
                            end=tmp_region.split(":")[1]
                            fo.write(ts_id+'\t'+start+'\t'+end+'\t'+seq[4]+'\t'+seq[5]+'\t'+seq[7]+'\t'+read_info+'\n')                
                fi.close()
                fo.close()

            ###################step2:获取transcript区域的reads
            print 'calculating EI:'
            #对每一个转录本进行遍历
            #read_wz_iso_info=pd.read_csv(read_wz_iso,sep='\t',header=None)
            freads=pd.read_csv('ts_mapping.zz.sorted.modified',sep='\t',header=None)
            fi=open(bed_in_dir)
            fo=open(EI_out_dir,'w')
            for line in fi:
                seq=line.rstrip().split('\t')
                ts_start=int(seq[1])
                ts_end=int(seq[2])
                ts_length=ts_end-ts_start
                transcript=seq[4]
                res_ad_in_ts=int(seq[6])
                if res_ad_in_ts > 0:
                    #比较转录本的长度和read的长度，若后者大于前者，剪切read
                    ts_mapped_zz_lines=freads.loc[freads.iloc[:,0]==transcript,:]
                    total_num=0
                    for idx in range(0,len(ts_mapped_zz_lines.index)):
                        read_seq=ts_mapped_zz_lines.iloc[idx,:]
                        read=read_seq[5]
                        read_start=int(read_seq[1])
                        read_end=int(read_seq[2])
                        read_length=read_end-read_start
                        if read_length > ts_length and read_start+ts_start < ts_end:
                            trimmed_read=read[0:(ts_end-(read_start+ts_start))]
                        elif read_length <= ts_length and read_start+ts_start < ts_end:
                            trimmed_read=read
                                                            
                        snv=read_seq[3]
                        AG_num=len(re.findall("AG:",snv))
                        #TC_num=len(re.findall("TC:",snv))
                        AA_num=len(re.findall("A",trimmed_read))
                        #AG mismatch > TC mismatch只算AG
                        total_num=total_num+AG_num+AA_num
                    if total_num!=0:
                        EI_of_this_ts=float(res_ad_in_ts)/total_num
                    else:
                        EI_of_this_ts='NA'
                    fo.write(line.rstrip()+'\t'+str(EI_of_this_ts)+'\n')
            fi.close()
            fo.close()

    print 'start'
    for sample in ['GW08','GW12','GW16_1_3', 'GW16_1_9', 'GW19_1_1', 'GW19_1_2', 'GW19_1_3', 'GW23_1_1', 'GW23_1_2', 'GW23_1_3', 'GW26_1_1']:
        dir='/disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/EI_in_reads/'+sample+'/'
        for as_type in ['A3','A5','AF','AL','MX','SE','RI']:
            #get_AEI('../'+as_type+'_total_events_wzResNum.txt',as_type+'_total_events_wzResNum_Over0.wzEI.txt')
            get_AEI(dir+as_type+'_as_events_wzResNum.txt',dir+as_type+'_as_events_wzResNum_Over0.wzEI.txt')
    print 'done'

main() 


