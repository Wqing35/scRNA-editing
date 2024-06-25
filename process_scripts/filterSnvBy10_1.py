#!/usr/bin/python2
# coding=utf-8

import subprocess,os,sys
import re
import time


def main(): 
    #import subprocess,os,sys
    def snv_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
        f1=open(bed_in_dir1)
        f2=open(bed_in_dir2)
        fo=open(bed_out_dir,'w')
        whole={}
        for line in f1:
            seq=line.rstrip().split('\t')
            whole[':'.join(seq[0:4])+':'+seq[5]]=int(seq[4])
        for line in f2:
            seq=line.rstrip().split('\t')
            if ':'.join(seq[0:4])+':'+seq[5] in whole:
                whole[':'.join(seq[0:4])+':'+seq[5]] +=int(seq[4])
            else:
                whole[':'.join(seq[0:4])+':'+seq[5]]=int(seq[4])
        lst=[]
        for one in whole:
            seq=one.split(':')
            lst.append([seq[0],int(seq[1]),int(seq[2]),seq[3],str(whole[one]),seq[4]])
        lst.sort()
        for one in lst:
            out=[]
            for i in one:
                out.append(str(i))
            fo.write('\t'.join(out)+'\n')

        f1.close()
        f2.close()
        fo.close()

    def annotate(bed_in_dir=0,bed_anno_dir=0,bed_out_dir=0):
        if bed_in_dir ==0 or bed_anno_dir ==0 or bed_out_dir ==0:
            print 'Please check the input directory! annotate(bed_in_dir,bed_anno_dir,bed_out_dir)'
            return 'Please check the input directory! annotate(bed_in_dir,bed_anno_dir,bed_out_dir)'
        else:
            fi=open(bed_in_dir) #must by sorted
            fa=open(bed_anno_dir)
            fo=open(bed_out_dir,'w')
            anno={}
            for line in fa:
                seq=line[0:-1].split('\t')
                try:
                    anno[seq[0]].append([int(seq[1])+1,int(seq[2]),seq[3],seq[4],seq[5]])	
                except Exception, e:
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

    def seperate(bed_in_dir,flag_out_dir,rp_out_dir,nonrp_out_dir,flag):
        fi=open(bed_in_dir)
        fo_flag=open(flag_out_dir,'w')
        fo_rp=open(rp_out_dir,'w')
        fo_nonrp=open(nonrp_out_dir,'w')
        for line in fi:
            if 'Simple_repeat' in line or 'Low_complexity' in line:
                next
            elif flag in line:
                fo_flag.write(line)
            elif 'Repeat_region' in line:
                fo_rp.write(line)
            else:
                fo_nonrp.write(line)
        fi.close()
        fo_flag.close()
        fo_rp.close()
        fo_nonrp.close()

    def get_snv_with_ad(snv_in_dir=0,snv_out_dir=0,flag=0):
        fi=open(snv_in_dir)
        fo=open(snv_out_dir,'w')
        
        for line in fi:
            seq=line.split('\t')
            try:
                if int(seq[4])>=int(flag):
                    fo.write(line)
            except Exception, e:
                print seq[0]+'\t'+seq[2]+"\twithout\tAD flag\n"		



        fi.close()
        fo.close()

    def snv_cluster(bed_in_dir=0,bed_out_dir=0,cluster_distance=-1,cluster_size=-1):
        fi=open(bed_in_dir)
        fo=open(bed_out_dir,'w')
        tmp='chr0:0:AA'
        limitdistance=int(cluster_distance)
        limitnum=int(cluster_size)
        lst=[]
        for line in fi:
                seq=line.split('\t')
                tmpseq=tmp.split(':')
                if seq[0]==tmpseq[0] and int(seq[2])-int(tmpseq[1])<=limitdistance and seq[3]==tmpseq[2]:
                    lst.append(line)
                else:
                    if len(lst)>=limitnum:
                            begin=float(lst[0].split('\t')[1])
                            end=float(lst[-1].split('\t')[2])
                            density=len(lst)/(end-begin)
                            for one in lst:
                                fo.write(one[0:-1]+'\t'+str(len(lst))+'\t'+str(density)+'\n')
                    lst=[]
                    lst.append(line)
                tmp=seq[0]+':'+seq[2]+':'+seq[3]
        if len(lst)>=limitnum:
                begin=float(lst[0].split('\t')[1])
                end=float(lst[-1].split('\t')[2])
                density=len(lst)/(end-begin)
                for one in lst:
                    fo.write(one[0:-1]+'\t'+str(len(lst))+'\t'+str(density)+'\n')
        fi.close()
        fo.close()

    def bed_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
        f1=open(bed_in_dir1)
        f2=open(bed_in_dir2)
        fo=open(bed_out_dir,'w')
        whole=[]
        for line in f1:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5]])
        for line in f2:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5]])
        whole.sort()
        old=set()
        for one in whole:
            if one[0]+':'+str(one[1]) not in old:
                fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\t'+one[4]+'\n')
                old.add(one[0]+':'+str(one[1]))
        f1.close()
        f2.close()
        fo.close()

    def combine_res(bed_in_dir1,bed_in_dir2,bed_in_dir3,bed_out_dir):
        f1=open(bed_in_dir1)
        f2=open(bed_in_dir2)
        f3=open(bed_in_dir3)
        fo=open(bed_out_dir,'w')
        whole=[]
        for line in f1:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5]])
        for line in f2:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5]])
        for line in f3:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5]])
        whole.sort()
        for one in whole:
            fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\t'+one[4]+'\n')
        f1.close()
        f2.close()
        f3.close()
        fo.close()

    def res_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
        f1=open(bed_in_dir1)
        f2=open(bed_in_dir2)
        fo=open(bed_out_dir,'w')
        whole={}
        for line in f1:
            seq=line.rstrip().split('\t')
            whole[':'.join(seq[0:4])+':'+seq[5]]=int(seq[4])
        for line in f2:
            seq=line.rstrip().split('\t')
            if ':'.join(seq[0:4])+':'+seq[5] in whole:
                pass
                #whole[':'.join(seq[0:4])+':'+seq[5]] +=int(seq[4])
            else:
                whole[':'.join(seq[0:4])+':'+seq[5]]=int(seq[4])
        lst=[]
        for one in whole:
            seq=one.split(':')
            lst.append([seq[0],int(seq[1]),int(seq[2]),seq[3],str(whole[one]),seq[4]])
        lst.sort()
        for one in lst:
            out=[]
            for i in one:
                out.append(str(i))
            fo.write('\t'.join(out)+'\n')

        f1.close()
        f2.close()
        fo.close()

    def o2b(bed_in,bed_out):
        fi=open(bed_in)
        fo=open(bed_out,'w')
        for line in fi:
            seq=line.rstrip().split('\t')
            #if float(seq[7])<0.2:
            fo.write('\t'.join(seq[0:6])+'\n')

    def get_depth(zz_in_dir=0,bed_in_dir=0,bed_out_dir=0):
        start_time = time.time()  

        fread=open(zz_in_dir)# './zz_folder/all.zz')
        fsnv=open(bed_in_dir) #'../bed_folder/try_new.bed') #   Hap_SRR521447.bed')
        fo=open(bed_out_dir,'w')#'./tmp/readspersite_new.zer','w')

        class Read:
            def __init__(self,read):
                
                self.snv=read.split('\t')[4].split(';')
                self.inter=read.split('\t')[3].split(';')
                self.direct=read.split('\t')[1]
            def locisin(self,loc):
                isin=0
                for inter in self.inter:
                    inter=inter.split(':')
                    if int(loc)<=int(inter[1]) and int(loc)>=int(inter[0]):
                        isin =1
                        break
                if isin ==0:
                    return 0
                elif isin ==1:
                    return 1
            def snvisin(self,snv):
                if snv in self.snv:
                    return 1
                else:
                    return 0
            def getmin(self):
                return int(self.inter[0].split(':')[0])
            def getmax(self):
                return int(self.inter[-1].split(':')[1])


        reads={}
        for line in fread:
            seq=line.split('\t')
            try:
                    reads[seq[0]].append(line[0:-1])
            except Exception,e :
                    #print seq[0]+' begin'
                    reads[seq[0]]=[line[0:-1]]


        top=0
        chrr=''
        for line in fsnv:
            seq=line.rstrip().split('\t')
            deep=0
            altdeep=0
            try:
                snv=seq[3]+':'+seq[2]
                if seq[0] != chrr:
                    top=0
                    chrr=seq[0]
                if seq[0] not in reads:
                    reads[seq[0]]=[]
                if top < len(reads[seq[0]]):
                    while seq[0]==chrr and top < len(reads[seq[0]]) and Read(reads[seq[0]][top]).getmax() < int(seq[2]):
                        top=top+1
                    point=top
                    while seq[0]==chrr and point < len(reads[seq[0]]) and Read(reads[seq[0]][point]).getmin() <= int(seq[2]):
                        if Read(reads[seq[0]][point]).locisin(seq[2]) ==1:
                            deep=deep+1
                        
                        if Read(reads[seq[0]][point]).snvisin(snv)==1:
                            altdeep=altdeep+1
                    
                        point=point+1
                fo.write(line[0:-1]+'\t'+str(altdeep)+':'+str(deep)+'\n')
            except Exception as e:
                print line
        fread.close()
        fsnv.close()
        fo.close()

        end_time = time.time()  # 结束计时
        print("get_depth 函数执行时间: {:.4f} 秒".format(end_time - start_time))  # 打印执行时间
        
    
    repeat='/disk1/wenqing/tmp_data/hg19_repeat.txt'
    output_directory='/disk1/wenqing/tmp_data/PFC_s2/filter_snv_for_new_res/'
    sample=['GW16_1_3','GW19_1_3','GW23_1_1','GW23_1_3']
    celltypes=['OPC']
    idx=10
    for one_sample in sample:
        for celltype in celltypes:
            output_directory_sub=output_directory+'/'+one_sample+'/'+celltype+'/'

            fi_regular_snv=open(output_directory_sub+'regular.snv')
            fo_regular_snv=open(output_directory_sub+'regular.snv.filtered.'+str(idx),'w')
            for line in fi_regular_snv:
                seq=line.rstrip().split('\t')
                coverage_reads=int(seq[4])
                if coverage_reads >= idx:
                    fo_regular_snv.write(line)
            fo_regular_snv.close()
            fi_regular_snv.close()


            """ if os.path.exists(output_directory_sub)==False:
                os.mkdir(output_directory_sub) """

            cluster_distance=200
            cluster_size_alu_ad1 = 3
            cluster_size_alu_ad2 = 2
            cluster_size_nalurp = 5
            cluster_size_nrp = 7
            #根据repeat file对上述流程的snv进行注释，那些区域是重复等
            annotate(output_directory_sub+'regular.snv.filtered.'+str(idx),repeat,output_directory_sub+'regular.snv.anno')    
            #根据注释将snv分成3个部分
            seperate(output_directory_sub+'regular.snv.anno',output_directory_sub+'regular.snv.anno.alu',output_directory_sub+'regular.snv.anno.nalurp',output_directory_sub+'regular.snv.anno.nrp','Alu')
            get_snv_with_ad(output_directory_sub+'regular.snv.anno.alu',output_directory_sub+'regular.snv.anno.alu.ad2',2)
            #筛选Alu snv的2条支线
            snv_cluster(output_directory_sub+'regular.snv.anno.alu',output_directory_sub+'regular_alu.res.ad1', cluster_distance, cluster_size_alu_ad1)
            snv_cluster(output_directory_sub+'regular.snv.anno.alu.ad2',output_directory_sub+'regular_alu.res.ad2', cluster_distance, cluster_size_alu_ad2)
            bed_or(output_directory_sub+'regular_alu.res.ad1',output_directory_sub+'regular_alu.res.ad2',output_directory_sub+'regular_alu.res')
            #筛选non-Alu snv
            snv_cluster(output_directory_sub+'regular.snv.anno.nalurp',output_directory_sub+'regular_nalurp.res', cluster_distance, cluster_size_nalurp)
            #筛选nrp snv
            snv_cluster(output_directory_sub+'regular.snv.anno.nrp',output_directory_sub+'regular_nrp.res', cluster_distance, cluster_size_nrp)
            #合并三类snv(cluster_distance+cluster size三者分开筛选后合并)
            combine_res(output_directory_sub+'regular_alu.res',output_directory_sub+'regular_nalurp.res',output_directory_sub+'regular_nrp.res',output_directory_sub+'regular_split.res')
            cluster_size_regular_max=max([cluster_size_alu_ad1,cluster_size_alu_ad2,cluster_size_nalurp,cluster_size_nrp])
            #print cluster_size_regular_max
            #三类snv共用cluster_distance+cluster_size_regular_max筛选条件
            #alu
            combine_res(output_directory_sub+'regular.snv.anno.alu',output_directory_sub+'regular.snv.anno.nalurp',output_directory_sub+'regular.snv.anno.nrp',output_directory_sub+'regular.snv.anno.rmsrp')
            snv_cluster(output_directory_sub+'regular.snv.anno.rmsrp',output_directory_sub+'regular_overall.res', cluster_distance, cluster_size_regular_max)
            res_or(output_directory_sub+'regular_split.res',output_directory_sub+'regular_overall.res',output_directory_sub+'regular.res')
        

            #对hyper snv采取与regular snv同样的操作
            """ annotate(output_directory_sub+'hyper_mskTC.snv',repeat,output_directory_sub+'hyper_mskTC.snv.anno')    
            seperate(output_directory_sub+'hyper_mskTC.snv.anno',output_directory_sub+'hyper_mskTC.snv.anno.alu',output_directory_sub+'hyper_mskTC.snv.anno.nalurp',output_directory_sub+'hyper_mskTC.snv.anno.nrp','Alu')
            snv_cluster(output_directory_sub+'hyper_mskTC.snv.anno.alu',output_directory_sub+'hyper_mskTC_alu.res', cluster_distance, cluster_size_alu_hp)
            snv_cluster(output_directory_sub+'hyper_mskTC.snv.anno.nalurp',output_directory_sub+'hyper_mskTC_nalurp.res', cluster_distance, cluster_size_nalurp_hp)
            snv_cluster(output_directory_sub+'hyper_mskTC.snv.anno.nrp',output_directory_sub+'hyper_mskTC_nrp.res', cluster_distance, cluster_size_nrp_hp)
            combine_res(output_directory_sub+'hyper_mskTC_alu.res',output_directory_sub+'hyper_mskTC_nalurp.res',output_directory_sub+'hyper_mskTC_nrp.res',output_directory_sub+'hyper_mskTC_split.res')        
            cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
            combine_res(output_directory_sub+'hyper_mskTC.snv.anno.alu',output_directory_sub+'hyper_mskTC.snv.anno.nalurp',output_directory_sub+'hyper_mskTC.snv.anno.nrp',output_directory_sub+'hyper_mskTC.snv.anno.rmsrp')
            snv_cluster(output_directory_sub+'hyper_mskTC.snv.anno.rmsrp',output_directory_sub+'hyper_mskTC_overall.res', cluster_distance, cluster_size_hyper_max)
            res_or(output_directory_sub+'hyper_mskTC_split.res',output_directory_sub+'hyper_mskTC_overall.res',output_directory_sub+'hyper_mskTC.res')

            #为何hyper_mskAG.snv.anno文件大小为0，后续仍能出结果？？？？？？
            annotate(output_directory_sub+'hyper_mskAG.snv',repeat,output_directory_sub+'hyper_mskAG.snv.anno')    
            seperate(output_directory_sub+'hyper_mskAG.snv.anno',output_directory_sub+'hyper_mskAG.snv.anno.alu',output_directory_sub+'hyper_mskAG.snv.anno.nalurp',output_directory_sub+'hyper_mskAG.snv.anno.nrp','Alu')
            snv_cluster(output_directory_sub+'hyper_mskAG.snv.anno.alu',output_directory_sub+'hyper_mskAG_alu.res', cluster_distance, cluster_size_alu_hp)
            snv_cluster(output_directory_sub+'hyper_mskAG.snv.anno.nalurp',output_directory_sub+'hyper_mskAG_nalurp.res', cluster_distance, cluster_size_nalurp_hp)
            snv_cluster(output_directory_sub+'hyper_mskAG.snv.anno.nrp',output_directory_sub+'hyper_mskAG_nrp.res', cluster_distance, cluster_size_nrp_hp)
            combine_res(output_directory_sub+'hyper_mskAG_alu.res',output_directory_sub+'hyper_mskAG_nalurp.res',output_directory_sub+'hyper_mskAG_nrp.res',output_directory_sub+'hyper_mskAG_split.res')        
            cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
            combine_res(output_directory_sub+'hyper_mskAG.snv.anno.alu',output_directory_sub+'hyper_mskAG.snv.anno.nalurp',output_directory_sub+'hyper_mskAG.snv.anno.nrp',output_directory_sub+'hyper_mskAG.snv.anno.rmsrp')
            snv_cluster(output_directory_sub+'hyper_mskAG.snv.anno.rmsrp',output_directory_sub+'hyper_mskAG_overall.res', cluster_distance, cluster_size_hyper_max)
            res_or(output_directory_sub+'hyper_mskAG_split.res',output_directory_sub+'hyper_mskAG_overall.res',output_directory_sub+'hyper_mskAG.res')

            snv_or(output_directory_sub+'hyper_mskTC.res',output_directory_sub+'hyper_mskAG.res',output_directory_sub+'hyper.res') """


            # subprocess.Popen('cp '+output_directory_sub+'/regular.res '+output+'/SPRINT_identified_regular.res',shell=True).wait()
            get_depth( output_directory_sub+'/all_combined.zz.sorted' , output_directory_sub+'/regular.res', output_directory_sub+'/regular.res.depth')
            #get_depth( output_directory_sub+'/all_combined.zz.sorted' , output_directory_sub+'/hyper.res', output_directory_sub+'/hyper.res.depth')

            #snv_or(output_directory_sub+'/regular.res',output_directory_sub+'/hyper.res',output_directory_sub+'/all.res') 
            #get_depth( output_directory_sub+'/all_combined.zz.sorted' , output_directory_sub+'/all.res', output_directory_sub+'/all.res.depth')
            #AD:DP：某位点RNA编辑率
            print one_sample+' '+celltype+' '+'finished !'

    
main()
