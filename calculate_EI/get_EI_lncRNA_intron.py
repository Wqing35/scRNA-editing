#!/usr/bin/python2
# coding=utf-8

import subprocess,os
import re
import pandas as pd

def main():

    #def get_AEI(bed_in_dir,bed_anno_dir,bed_out_dir,res_in_dir,AEI_out_dir)
    def get_AEI(element,region,depth_index):
        #bed_in_dir: res输入
        #bed_anno_dir：repeat.txt
        #bed_out_dir：res.anno
        #res_in_dir：alu区域的res 
        #element:regular 或all     
        #depth_index：过滤res dp的参数，用于resInAluWithADDP_over1改名


        #######step1:获取res的anno以得到alu区域的res
        #gene区域先不分alu区域，仅计算Editing Index（EI）  
        
        if os.path.exists('./'+element+'_combined.zz.sorted.modified_lncRNA_intron_trimmed_reads_OnlyAG')==False:
            ################step1：将all/regular_all_combined.zz.sorted中比对得到一个区域的多个reads分行放置
            if os.path.exists('./'+element+'_combined.zz.sorted.modified')==False:
                print 'Prepare reads_file:'

                out=[]
                fi=open(element+'_combined.zz')
                for line in fi:
                    seq=line.split('\t')
                    out.append([seq[0],int(seq[3].split(':')[0]),int(seq[3].split(':')[-1]),line])
                fi.close()
                fo=open(element+'_combined.zz.sorted', 'w')
                out.sort()
                for one in out:
                    fo.write(one[3]) 

                fi.close()
                fo.close()

                fi=open(element+'_combined.zz.sorted')
                fo=open(element+'_combined.zz.sorted.modified','w')

                for line in fi:
                    seq=line.rstrip().split("\t")
                    regions=seq[3].split(";")
                    barcode=seq[8].split('_')[2]
                    if len(regions) ==1:
                        regions=regions[0]
                        start=regions.split(":")[0]
                        end=regions.split(":")[1]
                        fo.write(seq[0]+'\t'+start+'\t'+end+'\t'+seq[4]+'\t'+seq[5]+'\t'+seq[7]+'\t'+barcode+'\n')
                    else:
                        for i in range(0,len(regions)):
                            tmp_region=regions[i]
                            start=tmp_region.split(":")[0]
                            end=tmp_region.split(":")[1]
                            fo.write(seq[0]+'\t'+start+'\t'+end+'\t'+seq[4]+'\t'+seq[5]+'\t'+seq[7]+'\t'+barcode+'\n')                
                fi.close()
                fo.close()

            ###################step2:获取gene区域的reads
            print 'Get reads in lncRNA intron:'

            #筛选位于gene区域的reads
            subprocess.Popen('bash '+'/disk1/wenqing/scRNA-editing/calculate_EI/get_reads_lncRNA_Intron'+'.sh',shell=True).wait()
            #bedtools intersect -wb -a regular_all_combined.zz.sorted.modified -b ~/tmp_data/Alu_seq/AsInAlus_modified.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"_"$12}' > regular_all_combined.zz.sorted.modified_InAlus


            ###################step4:根据gene范围，剪切reads长度
            print 'Trim reads according to lncRNA intron:'

            fi=open(element+'_combined.zz.sorted.modified_lncRNA_intron')
            fo=open(element+'_combined.zz.sorted.modified_lncRNA_intron_trimmed_reads_OnlyAG','w')

            for line in fi:
                seq=line.rstrip().split("\t")
                read=seq[5]
                reads_start=int(seq[1])
                reads_end=int(seq[2])
                gene_start=int(seq[8])
                gene_end=int(seq[9])

                #剪切reads的几种可能的场景
                #1.该条read包含该gene
                if gene_start > reads_start and gene_end < reads_end:
                    new_start=gene_start
                    new_end=gene_end
                    trimmed_read=read[(gene_start-reads_start-1):(gene_end-reads_start)]
                #2.该条read被包含于该gene
                elif gene_start < reads_start and gene_end > reads_end:
                    new_start=reads_start
                    new_end=reads_end
                    trimmed_read=read[0:(reads_end-reads_start)]
                #3.该read的下游区域未覆盖到gene
                elif gene_start == reads_start and gene_end > reads_end:
                    new_start=reads_start
                    new_end=reads_end
                    trimmed_read=read[0:(reads_end-reads_start)] 
                #4.该read的上游区域未覆盖到gene
                elif gene_start < reads_start and gene_end == reads_end:
                    new_start=reads_start
                    new_end=reads_end
                    trimmed_read=read[0:(reads_end-reads_start)] 
                #5.该read下游过覆盖gene
                elif gene_start == reads_start and gene_end < reads_end:
                    new_start=reads_start
                    new_end=gene_end
                    trimmed_read=read[0:(gene_end-gene_start)]   
                #6.该read上游过覆盖gene
                elif gene_start > reads_start and gene_end == reads_end:
                    new_start=gene_start
                    new_end=gene_end
                    trimmed_read=read[(gene_start-reads_start):(reads_end-reads_start)]                 
                #7.该reads上游过覆盖、下游未覆盖
                elif gene_start > reads_start and gene_end > reads_end:
                    new_start=gene_start
                    new_end=reads_end
                    trimmed_read=read[(gene_start-reads_start):(reads_end-reads_start)]                 
                #8.该reads上游未覆盖、下游过覆盖
                elif gene_start < reads_start and gene_end < reads_end:
                    new_start=reads_start
                    new_end=gene_end
                    trimmed_read=read[0:(gene_end-reads_start)]                 

                snv=seq[3]
                AG_num=len(re.findall("AG:",snv))
                TC_num=len(re.findall("TC:",snv))
                As_num=AG_num+TC_num
                AA_num=len(re.findall("A",trimmed_read))
                #total_num=As_num+AA_num
                #AG mismatch > TC mismatch只算AG
                total_num=AG_num+AA_num
                gene_mark=seq[7]+'_'+seq[8]+'_'+seq[9]
                fo.write(seq[0]+'\t'+str(new_start)+'\t'+str(new_end)+'\t'+seq[3]+'\t'+seq[4]+'\t'+trimmed_read+'\t'+seq[10]+'\t'+seq[11]+'\t'+str(total_num)+'\t'+gene_mark+'\n')

            fi.close()
            fo.close()

        ###########step5:获取AEI
        print 'Get GEI:'

        EI_out_dir=element+'_EI_lncRNA_intron_over'+str(depth_index)+'.txt'
        gene_counts_in_dir=element+'_combined.zz.sorted.modified_lncRNA_intron_trimmed_reads_OnlyAG'

        fi=open(gene_counts_in_dir)
        #######待改！！！！！！！！
        gene_dict1={}
        for line in fi:
            seq=line.rstrip().split('\t')
            gene_id=seq[6]  
            res_As_num=int(seq[7])
            reads_As_num=int(seq[8])
            try:
                gene_dict1[gene_id][0]=res_As_num
                gene_dict1[gene_id][1]=gene_dict1[gene_id][1]+reads_As_num
            except Exception, e:
                gene_dict1[gene_id]=[0,0]
                gene_dict1[gene_id][0]=res_As_num
                gene_dict1[gene_id][1]=reads_As_num
        fi.close()

        EI=[]
        for key in gene_dict1.keys():
            if gene_dict1[key][0]==0 or gene_dict1[key][1]==0:
                EI.append(float(0))
            else:
                EI.append(float(gene_dict1[key][0])/int(gene_dict1[key][1]))
        
        data=pd.DataFrame()
        #data["gene_name"]=gene_dict1.keys()
        data["EI"]=EI
        data.index=gene_dict1.keys()

        data.to_csv(EI_out_dir, sep='\t')

    print 'start'
    #dp_idx=[0,1,2,3]
    #dp_idx=[1]
    #for idx in dp_idx:
    idx=0
    get_AEI('regular','lncRNA_intron',idx)
    print 'done'

main() 

