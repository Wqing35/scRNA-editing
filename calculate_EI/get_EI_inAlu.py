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
        #bed_in_dir: res输入
        #bed_anno_dir：repeat.txt
        #bed_out_dir：res.anno
        #res_in_dir：alu区域的res 
        #element:regular 或all     
        #depth_index：过滤res dp的参数，用于resInAluWithADDP_over1改名


        #######step1:获取res的anno以得到alu区域的res
        if os.path.exists('resInAluWithADDP_WizDepth')==False:
            print 'Get RES in Alu:'

            fi=open(element+'_merged.res') #must by sorted   bed_in_dir
            fa=open('/disk1/wenqing/tmp_data/hg19_repeat.txt') # bed_anno_dir
            fo=open(element+'_merged.res.anno','w') # bed_out_dir
            anno={}
            for line in fa:
                seq=line[0:-1].split('\t')
                try:
                    anno[seq[0]].append([int(seq[1])+1,int(seq[2]),seq[3],seq[4],seq[5]])	
                except Exception as e:
                    anno[seq[0]]=[ [int(seq[1])+1,int(seq[2]),seq[3],seq[4],seq[5]]  ]
            for a in anno:
                anno[a].sort()
            top=0
            point=0
            lastchr=''
            for line in fi:
                
                seq = line.rstrip().split('\t')
                #if len(seq[0]) <=5:
                if 1==1:
                    output=line.replace('\n','')
                    if seq[0]==lastchr:
                        1+1
                    else:
                        top=0
                        lastchr=seq[0]
                    if seq[0] not in anno:
                        anno[seq[0]]=[ [0,0,0,0,0]  ]
                    if top<len(anno[seq[0]]):			
                        while anno[seq[0]][top][1] < int(seq[2]) :
                
                            top=top+1
                            if top >= len(anno[seq[0]]):
                                break
                        point=top
                        if point < len(anno[seq[0]]):
                            while anno[seq[0]][point][0] <= int(seq[2])     and     anno[seq[0]][point][1] >= int(seq[2]) :
                                    if anno[seq[0]][point][2]+'\t'+anno[seq[0]][point][3]+'\t'+anno[seq[0]][point][4] not in output:
                                        output=output+'\t'+anno[seq[0]][point][2]+'\t'+anno[seq[0]][point][3]+'\t'+anno[seq[0]][point][4]
                                    point=point+1
                                    if point >= len(anno[seq[0]]):
                                        break
                
                    fo.write(output+'\n')
            fi.close()
            fa.close()
            fo.close()   

            subprocess.Popen('bash '+'/disk1/wenqing/scRNA-editing/caculate_EI/get_resInAlus_'+element+'.sh',shell=True).wait()


        fi=open('resInAluWithADDP_WizDepth')
        fo=open('resInAluWithADDP_over'+str(depth_index),'w')
        for line in fi:
            seq=line.rstrip().split('\t')
            dp=int(seq[7])
            if dp > depth_index:
                fo.write(seq[0]+'\t'+seq[1]+'\t'+seq[2]+'\t'+seq[3]+'\t'+seq[4]+'\t'+seq[5]+'\t'+seq[6]+'\n')
        fi.close()
        fo.close()

        if os.path.exists('./resInAluWithADDP_WizDepth'):
            os.remove('resInAluWithADDP')
            os.remove('dp.txt')

        #subprocess.Popen('mv resInAluWithADDP'+' ./resInAluWithADDP_over'+str(depth_index),shell=True).wait()   
        
        if os.path.exists('./'+element+'_combined.zz.sorted.modified_InAlus_trimmed_reads_OnlyAG')==False:
            if os.path.exists('./'+element+'_combined.zz.sorted.modified')==False:
            ################step2：将all/regular_all_combined.zz.sorted中比对得到一个区域的多个reads分行放置
                print 'Prepare reads_file:'

                out=[]
                fi=open(element+'_combined.zz')
                for line in fi:
                    seq=line.rstrip().split('\t')
                    out.append([seq[0],int(seq[3].split(':')[0]),int(seq[3].split(':')[-1]),line])
                fi.close()
                fo=open(element+'_combined.zz.sorted', 'w')
                out.sort()
                for one in out:
                    fo.write(one[3]) 

                fi.close()
                fo.close()

                #将匹配到多个region的reads分行放置
                fi=open(element+'_combined.zz.sorted')
                fo=open(element+'_combined.zz.sorted.modified','w')

                for line in fi:
                    seq=line.rstrip().split("\t")
                    regions=seq[3].split(";")
                    if len(regions) ==1:
                        regions=regions[0]
                        start=regions.split(":")[0]
                        end=regions.split(":")[1]
                        fo.write(seq[0]+'\t'+start+'\t'+end+'\t'+seq[4]+'\t'+seq[5]+'\t'+seq[7]+'\n')
                    else:
                        for i in range(0,len(regions)):
                            tmp_region=regions[i]
                            start=tmp_region.split(":")[0]
                            end=tmp_region.split(":")[1]
                            fo.write(seq[0]+'\t'+start+'\t'+end+'\t'+seq[4]+'\t'+seq[5]+'\t'+seq[7]+'\n')                
                fi.close()
                fo.close()
                

            ###################step3:获取alu区域的reads
            print 'Get reads in Alu:'

            #筛选位于Alu区域的reads
            subprocess.Popen('bash '+'/disk1/wenqing/scRNA-editing/caculate_EI/get_readsInAlus_'+element+'.sh',shell=True).wait()
            #bedtools intersect -wb -a regular_all_combined.zz.sorted.modified -b ~/tmp_data/Alu_seq/AsInAlus_modified.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"_"$12}' > regular_all_combined.zz.sorted.modified_InAlus


            ###################step4:根据alu区域范围，剪切reads长度
            print 'Trim reads according to Alu:'

            fi=open(element+'_combined.zz.sorted.modified_InAlus')
            fo=open(element+'_combined.zz.sorted.modified_InAlus_trimmed_reads_OnlyAG','w')

            for line in fi:
                seq=line.rstrip().split("\t")
                read=seq[5]
                reads_start=int(seq[1])
                reads_end=int(seq[2])
                alu_start=int(seq[7])
                alu_end=int(seq[8])

                if alu_start > reads_start and alu_end < reads_end:
                    # 情况1：alu完全包含在read内部
                    new_start = alu_start - reads_start
                    new_end = alu_end - reads_start
                    trimmed_read = read[new_start - 1:new_end]
                elif alu_start <= reads_start and alu_end >= reads_end:
                    # 情况2：read完全被alu覆盖
                    new_start = 0
                    new_end = reads_end - reads_start
                    trimmed_read = read[:new_end]
                elif alu_start == reads_start:
                    # 情况3：alu起始于read开始处，无论alu是否在read内结束
                    new_end = alu_end - reads_start
                    trimmed_read = read[:new_end]
                elif alu_end == reads_end:
                    # 情况4：alu结束于read结束处，无论alu是否从read起始处开始
                    new_start = alu_start - reads_start
                    trimmed_read = read[new_start:]
                else:
                    # 这里可以添加一个默认处理或错误提示，尽管按照你的说明这种情况不应出现
                    print "Unexpected overlap scenario between ALU and read."
                    exit(1)

                snv=seq[3]
                AG_num=len(re.findall("AG:",snv))
                TC_num=len(re.findall("TC:",snv))
                As_num=AG_num+TC_num
                AA_num=len(re.findall("A",trimmed_read))
                #total_num=As_num+AA_num
                #AG mismatch > TC mismatch只算AG
                total_num=AG_num+AA_num
                fo.write(seq[0]+'\t'+seq[1]+'\t'+seq[2]+'\t'+seq[3]+'\t'+seq[4]+'\t'+trimmed_read+'\t'+seq[9]+'\t'+str(total_num)+'\n')

            fi.close()
            fo.close()
        
        ##得到modified.trimmed_AG后，去掉当前路径下的zz文件，以释放存储空间
        if os.path.exists('./'+element+'_combined.zz.sorted.modified_InAlus_trimmed_reads_OnlyAG'):
            os.remove('./'+element+'_combined.zz.sorted')
            os.remove('./'+element+'_combined.zz.sorted.modified')
            os.remove('./'+element+'_combined.zz.sorted.modified_InAlus')


        ###########step5:获取AEI
        print 'Get AEI:'

        AEI_out_dir='AEI_AG_over'+str(depth_index)+'.txt'
        res_in_dir='resInAluWithADDP_over'+str(depth_index)

        fo=open(AEI_out_dir,'w')
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
        
        AverAEIByAlu={}
        for key in sumOfAsByAlu.keys():
            AG_num=int(sumOfAGByAlu[key])
            As_num=int(sumOfAsByAlu[key])
            AverAEIByAlu[key]=float(AG_num)/As_num

        for one in AverAEIByAlu.keys():
            fo.write(one+'\t'+str(AverAEIByAlu[one])+'\n')
        fo.close()

    print 'start'
    #dp_idx=[0,1,2,3]
    dp_idx=[0]
    for idx in dp_idx:
        get_AEI('regular',idx)
    print 'done'

main() 