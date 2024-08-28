#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import subprocess,os,sys
import re
import numpy as np
import pandas as pd
from datetime import datetime
from collections import Counter


def main(): 
    #import subprocess,os,sys
    
    print ''
    print "##############################################################################################"
    print ''
    print "   SPRINT: SNP-free RNA editing Identification Toolkit"
    print ""
    print "   http://sprint.tianlab.cn/SPRINT/"
    print ""
    print "   Please contact 15110700005@fudan.edu.cn when questions arise."
    print ""
    print "##############################################################################################"
    
    
    def help_doc():
        print ""
        print "   Attention:"
        print ""
        print "      Before using 'sprint main', please use 'sprint prepare' to build mapping index."
        print ""
        print "   Usage:"
        print ""
        print "      sprint main   [options]   reference_genome(.fa)   output_path   hisat2_path   samtools_path"
        print ""
        print "      options:"
        print "         -1       read(.fq)       # Required !!!"
        print "         -r1      barcode+UMI(.fq) # R1 file from cellranger. Required!!!!"
        print "         -rp      repeat_file      # Optional, you can download it from http://sprint.software/SPRINT/dbrep/"
        print "         -ss      INT              # when input is strand-specific sequencing data, please clarify the direction of read1. [0 for antisense; 1 for sense] (default is 0)"
        #print "         -b       INT             # the format of read file [0: fq, 1: bam] (default is 0)"
        print "         -c       INT              # Remove the fist INT bp of each read (default is 0)"
        print "         -cr      INT              # extract barcode sequence from outputs of cellranger"
        print "         -br      barcode.txt      # input barcode output from cellranger to filter original barcode"
        print "         -p       INT              # Mapping CPU (default is 1)"
        print "         -cd      INT              # The distance cutoff of SNV duplets (default is 200)"
        print "         -csad1   INT              # Regular - [-rp is required] cluster size - Alu - AD >=1 (default is 3)"
        print "         -csad2   INT              # Regular - [-rp is required] cluster size - Alu - AD >=2 (default is 2)"
        print "         -csnar   INT              # Regular - [-rp is required] cluster size - nonAlu Repeat - AD >=1 (default is 5)"
        print "         -csnr    INT              # Regular - [-rp is required] cluster size - nonRepeat - AD >=1 (default is 7)"
        print "         -csrg    INT              # Regular - [without -rp] cluster size - AD >=1 (default is 5)"
        print "         -csahp   INT              # Hyper - [-rp is required] cluster size - Alu - AD >=1 (default is 5)"
        print "         -csnarhp INT              # Hyper - [-rp is required] cluster size - nonAlu Repeat - AD >=1 (default is 5)"
        print "         -csnrhp  INT              # Hyper - [-rp is required] cluster size - nonRepeat - AD >=1 (default is 5)"
        print "         -cshp    INT              # Hyper - [without -rp] cluster size - AD >=1 (default is 5)"
        print ""
        print "   Example:"
        print ""
        print "       sprint main -rp hg19_repeat.txt -c 6 -cr 16 -p 6 -1 r2.fq -r1 r1.fq -br ./barcode.txt hg19.fa output ./hisat2-0.7.12/hisat2 ./samtools-1.2/samtools"
        print ""
        print "       Notes: Default protocol of strand-specific RNA-seq is dUTP (read1: '-'; read2: '+')"
        print ""
        #print sys.argv[0]
        
        sys.exit(0)
    
    
    if len(sys.argv)<2:
        #print sys.argv[0]
        help_doc()
    
    
    #read_format:输入文件的类型，此处默认是fq，且不可指定
    read_format=0
    cutbp=0
    cutbp_r1=0
    cluster_distance=200
    overlap_criterion_set=[0,1]
    cluster_size_alu_ad1 = 3
    cluster_size_alu_ad2 = 2
    cluster_size_nalurp = 5
    cluster_size_nrp = 7
    cluster_size_rg = 5
    cluster_size_hp = 5
    cluster_size_alu_hp = 5
    cluster_size_nalurp_hp = 5
    cluster_size_nrp_hp = 5
    strand_specify=0


    mapcpu = 1
    var_limit=20
    poly_limit=10
    rm_multi=0
    
    paired_end=True
    repeat=False
    #options：参数所在位置序号的列表
    options=[]
    read2=''
    read1=''
    read_r1=''
    barcode_ref=''
    #print sys.argv
    #sprint main -rp hg19_repeat.txt -c 6 -p 6 -1 read1.fq -2 read2.fq hg19.fa output ./hisat2-0.7.12/hisat2 ./samtools-1.2/samtools
    i=1
    while i< len(sys.argv):
        #if sys.argv[i]=='-b':
        #    try:
        #        read_format=int(sys.argv[i+1])
        #        options.append(i)
        #        options.append(i+1)
        #    except Exception, e:
        #        print 'options error!'
        #        help_doc()

        #读取read1，无输入则报错、并给出帮助文档
        if sys.argv[i]=='-1':
            try:
                read1=sys.argv[i+1]
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
        #判断是否双端测序
        elif sys.argv[i]=='-2':
            paired_end=True
            try:
                read2=sys.argv[i+1]
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-br':
            try:
                barcode_ref=sys.argv[i+1]
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-rp':
            try:
                repeat=sys.argv[i+1]
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        #是否特异链测序
        elif sys.argv[i]=='-ss':
            try:
                strand_specify=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        #清理前几个碱基
        elif sys.argv[i]=='-c':
            try:
                cutbp=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        #提取barcode
        elif sys.argv[i]=='-cr':
            try:
                cutbp_r1=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        #分配cpu个数
        elif sys.argv[i]=='-p':
            try:
                mapcpu=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-cd':
            try:
                cluster_distance=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csad1':
            try:
                cluster_size_alu_ad1=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csad2':
            try:
                cluster_size_alu_ad2=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csnar':
            try:
                cluster_size_nalurp=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csnr':
            try:
                cluster_size_nrp=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csrg':
            try:
                cluster_size_rg=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-cshp':
            try:
                cluster_size_hp=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csahp':
            try:
                cluster_size_alu_hp=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csnarhp':
            try:
                cluster_size_nalurp_hp=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
        elif sys.argv[i]=='-csnrhp':
            try:
                cluster_size_nrp_hp=int(sys.argv[i+1])
                options.append(i)
                options.append(i+1)
            except Exception, e:
                print 'options error!'
                help_doc()
                exit()
    
    
        i += 1
    
    all_argv=[]
    i=1
    #补全命令中未手动捕捉的参数位置
    while i< len(sys.argv):
        if i not in options:
            all_argv.append(i)
        i=i+1
    
    if len(all_argv)!=4 or read1=='':
        help_doc()
        exit()
    
    #sys.argv：取参数；refgenome：参考序列文件
    refgenome=sys.argv[all_argv[0]]
    output=sys.argv[all_argv[1]]+'/'
    #temp：output目录下结果文件目录
    tmp=output+'/tmp/'

    #hisat2、samtools所在路径
    hisat2=sys.argv[all_argv[2]]
    samtools=sys.argv[all_argv[3]]
    
    #存储路径若不存在，则手动创建
    if os.path.exists(output)==False:
            os.mkdir(output)
    if os.path.exists(tmp)==False:
            os.mkdir(tmp)
    
    frc=open(tmp+'PARAMETER.txt','w')
    frc.write(sys.argv[0])
    for one in sys.argv[1:]:
        frc.write('   '+one)
    frc.write('\n')
    frc.close()

    def extract_reads(barcode_file, read1, read2, output_read1, output_read2, cutbp_r1):  
        print('extracting your celltype reads!')
        barcodes = set(open(barcode_file).read().strip().split("\n"))  
        read2_handle = open(read2, "r")  
        read1_handle = open(read1, "r")  
        output_read1_handle = open(output_read1, "w")  
        output_read2_handle = open(output_read2, "w")  
    
        idd=1
        while 1==1:
            read2_header = read2_handle.readline()  
            read2_seq = read2_handle.readline()  
            read2_plus = read2_handle.readline()  
            read2_qual = read2_handle.readline() 

            read1_header = read1_handle.readline()  
            read1_seq = read1_handle.readline()  
            read1_plus = read1_handle.readline()  
            read1_qual = read1_handle.readline() 
    
            if not read2_header:  
                break  # 到达 read2 文件的末尾  
    
            #read2_seq = read2_seq.strip()  
            if read2_seq[:cutbp_r1] in barcodes: 
                # 提取 read1 的对应读段   
                output_read1_handle.write('@id_'+str(idd)+'_'+read2_seq[:8]+'_read1'+'\n')  
                output_read1_handle.write(read1_seq)  
                output_read1_handle.write(read1_plus)  
                output_read1_handle.write(read1_qual)  
    
                # 写入 read2 的读段  
                output_read2_handle.write('@id_'+str(idd)+'_'+read2_seq[:8]+'_read2'+'\n')  
                output_read2_handle.write(read2_seq[(cutbp_r1+9):])  
                output_read2_handle.write(read2_plus)  
                output_read2_handle.write(read2_qual[(cutbp_r1+9):])  
            idd=idd+1

        read2_handle.close()  
        read1_handle.close()  
        output_read1_handle.close()  
        output_read2_handle.close()  

    def cut(fq_in_dir=0,fq_r1_in_dir=0,barcode_in_dir=0,fq_out_dir=0,cutnum=0,cutnum_r1=0,name='read1'):
        #barcode读取
        final_barcode=pd.read_csv(barcode_in_dir,sep='\n',header=None).replace('-1','',regex=True)
        fi_r1=open(fq_r1_in_dir)
        cutnum_r1=int(cutnum_r1)
        line1_r1=fi_r1.readline()
        line2_r1=fi_r1.readline()
        line3_r1=fi_r1.readline()
        line4_r1=fi_r1.readline()
        #read读取
        fi=open(fq_in_dir)
        fo=open(fq_out_dir,'w')
        cutnum=int(cutnum)
        line1=fi.readline()
        line2=fi.readline()
        line3=fi.readline()
        line4=fi.readline()
        idd=1
        while line1 !='' and line1_r1 !='':
            barcode=line2_r1[:cutnum_r1]
            CELL_TAG=''
            #添加一步，在比对前就根据表达谱的barcode针对reads进行过滤
            if barcode in np.array(final_barcode):
                if "XC:Z:" in line1:
                    seq=line1.split('_')
                    for one in seq:
                        if one[:5]=='XC:Z:':
                            CELL_TAG=one        
                if CELL_TAG !='':
                    fo.write('@id_'+str(idd)+'_'+CELL_TAG+'_'+barcode+'_'+name+'\n')
                else:
                    fo.write('@id_'+str(idd)+'_'+barcode+'_'+name+'\n')
                fo.write(line2[cutnum:])
                fo.write(line3)
                fo.write(line4[cutnum:])
                #read迭代
            line1=fi.readline()
            line2=fi.readline()
            line3=fi.readline()
            line4=fi.readline()
            idd=idd+1
            #barcode迭代
            line1_r1=fi_r1.readline()
            line2_r1=fi_r1.readline()
            line3_r1=fi_r1.readline()
            line4_r1=fi_r1.readline()
        fi.close()
        fi_r1.close()
        fo.close()    


    #此方法用于判断fq文件第四行碱基质量是以33还是以64开始，若以64为准，则后续碱基质量阈值定为89
    def get_baseq_cutoff(fq_in_dir=0,cutoff_out_dir=0):
        fi=open(fq_in_dir)
        fo=open(cutoff_out_dir,'w')
        line1=fi.readline()
        line2=fi.readline()
        line3=fi.readline()
        line4=fi.readline()
        did=0
        while line1 !='':
            if did==1:
                break
            qua=line4[0:-1]
            for i in qua:
                if ord(i) > 76:
                    fo.write('89') #64+25
                    did=1
                    break
                if ord(i) < 60:
                    fo.write('58') #33+25
                    did=1
                    break

            line1=fi.readline()
            line2=fi.readline()
            line3=fi.readline()
            line4=fi.readline()
        fi.close()
        fo.close()
    
    def fq2sam(TAG,paired_end,read1,read2,tmp,refgenome,hisat2,samtools,mapcpu, read_format):
        refgenome='/disk1/wenqing/tmp_data/ref/hisat2_ref/'+TAG
        """ if TAG=='transcript':
            refgenome='/disk1/wenqing/tmp_data/ref/hisat2_ref/'+TAG+'/hg19._trans'
        elif TAG=='genome_mskAG':
            refgenome='/disk1/wenqing/tmp_data/ref/hisat2_ref/'+TAG+'/hg19_mskAG'
        elif TAG=='genome_mskTC':
            refgenome='/disk1/wenqing/tmp_data/ref/hisat2_ref/'+TAG+'/hg19_mskTC'
        elif TAG=='transcript_mskAG':
            refgenome='/disk1/wenqing/tmp_data/ref/hisat2_ref/'+TAG+'/hg19_trans_mskAG'
        elif TAG=='transcript_mskTC':
            refgenome='/disk1/wenqing/tmp_data/ref/hisat2_ref/'+TAG+'/hg19_trans_mskTC' """

        if paired_end==True:
            mapcpu=max([int(int(mapcpu)/2.0),1])
        ori_tmp=tmp
        tmp=tmp+'/'+TAG+'/'
        if os.path.exists(tmp)==False:
            os.mkdir(tmp)
        
        #hisat2比对
        step1_1=subprocess.Popen(hisat2+' -q -x '+refgenome+' -U '+read1+' -S '+tmp+'name_read1.sam',shell=True)
        if paired_end==True:
            print('pair-end!')
            step1_2=subprocess.Popen(hisat2+' -q -x '+refgenome+' -U '+read2+' -S '+tmp+'name_read2.sam',shell=True)
        step1_1.wait()
        if paired_end==True:
            step1_2.wait()


        if os.path.exists(tmp+'name_read1.sam'):
            if os.path.exists(ori_tmp+'cut_read1.fastq'):
                    os.remove(ori_tmp+'cut_read1.fastq')
            if os.path.exists(ori_tmp+'cut_read2.fastq'):
                    os.remove(ori_tmp+'cut_read2.fastq') 

        #sam转bam
        step1_3=subprocess.Popen(samtools+' view -bS '+tmp+'name_read1.sam >'+tmp+'name_read1.bam',shell=True)
        if paired_end==True:
            step1_4=subprocess.Popen(samtools+' view -bS '+tmp+'name_read2.sam >'+tmp+'name_read2.bam',shell=True)
        step1_3.wait()
        if paired_end==True:
            step1_4.wait()
        #排序
        if paired_end==True:
            step1_5=subprocess.Popen(samtools+' sort '+tmp+'name_read1.bam '+tmp+'name_read1_sorted',shell=True)
            step1_6=subprocess.Popen(samtools+' sort '+tmp+'name_read2.bam '+tmp+'name_read2_sorted',shell=True)
            step1_5.wait()
            step1_6.wait()
            step1_7=subprocess.Popen(samtools+' merge -f '+tmp+'all.bam '+tmp+'name_read1.bam '+tmp+'name_read2.bam',shell=True)
            step1_7.wait()
            #all.bam生成后，删除read1与read2对应的sam、bam、sorted.bam
            if os.path.exists(tmp+'all.bam'):
                    if os.path.exists(tmp+'name_read1.sam'):
                            os.remove(tmp+'name_read1.sam')
                    if os.path.exists(tmp+'name_read1.bam'):
                            os.remove(tmp+'name_read1.bam')
                    if os.path.exists(tmp+'name_read1_sorted.bam'):
                            os.remove(tmp+'name_read1_sorted.bam')
                    if os.path.exists(tmp+'name_read2.sam'):
                            os.remove(tmp+'name_read2.sam')
                    if os.path.exists(tmp+'name_read2.bam'):
                            os.remove(tmp+'name_read2.bam')
                    if os.path.exists(tmp+'name_read2_sorted.bam'):
                            os.remove(tmp+'name_read2_sorted.bam')  
        else:
            step1_5=subprocess.Popen(samtools+' sort '+tmp+'name_read1.bam '+tmp+'all',shell=True)
            step1_5.wait()
            if os.path.exists(tmp+'all.bam'):
                    if os.path.exists(tmp+'name_read1.sam'):
                            os.remove(tmp+'name_read1.sam')
                    if os.path.exists(tmp+'name_read1.bam'):
                            os.remove(tmp+'name_read1.bam')
        step2_1=subprocess.Popen(samtools+' view -h -o '+tmp+'all.sam '+tmp+'all.bam',shell=True)
        step2_1.wait()
        subprocess.Popen('cp '+tmp+'./all.sam '+ori_tmp+'/'+TAG+'_all.sam',shell=True).wait()
        if os.path.exists(tmp+'all.sam'):
            os.remove(tmp+'all.sam')

    def umsam2fq(sam_in_dir,fq_out_dir):
        fi=open(sam_in_dir)
        fo=open(fq_out_dir,'w')
        for line in fi:
            seq=line.rstrip().split('\t')
            if line[0] !='@' and len(bin(int(seq[1])))>=5 and bin(int(seq[1]))[-3]=='1':
                if len(bin(int(seq[1])))>=9 and bin(int(seq[1]))[-7]=='1':
                    seq[0]=seq[0][0:-2]+'_1'
                elif len(bin(int(seq[1])))>=10 and bin(int(seq[1]))[-8]=='1':
                    seq[0]=seq[0][0:-2]+'_2'
                elif line[-1]=='1':
                    seq[0]=seq[0][0:-2]+'_1'
                elif line[-1]=='2':
                    seq[0]=seq[0][0:-2]+'_2'
                fo.write('@'+seq[0]+'\n'+seq[9]+'\n+\n'+seq[10]+'\n')
        fo.close()
        fi.close()

    def antisense_reverse(read):
        read=read.upper()
        read_change_base=""
        for one in read:
            if one == 'A':
                read_change_base += 'T'
            elif one == 'C':
                read_change_base += 'G'
            elif one == 'G':
                read_change_base += 'C'
            elif one == 'T':
                read_change_base += 'A'
            else:
                read_change_base += 'N'
        read_reverse=read_change_base[::-1]
        return read_reverse


    def maskfq(fq_in_dir,mask_from,mask_to):
        mask_from=mask_from.upper()
        mask_to=mask_to.upper()
        fi=open(fq_in_dir)
        #重写文件名
        fo=open(fq_in_dir[0:-3]+'_'+mask_from+'_to_'+mask_to+'.fq','w')
        line1=fi.readline().replace('\n','')
        line2=fi.readline().upper().replace('\n','')
        line3=fi.readline().replace('\n','')
        line4=fi.readline().replace('\n','')
        while line1 !='':
            #若是fq文件转换描述行显示该read为read1，则将其转换为反义链；用作下文针对mask_from的碱基计数（单细胞无链特异性的区分，此处会与hyper的结果不同）
            if line1[-1]=='1':
                line2=antisense_reverse(line2)
                line4=line4[::-1]
            record="1"
            line2_new=""
            for one in line2.replace('\n',''):
                if one==mask_from:
                    record=record+'1'
                else:
                    record=record+'0'
        
            #mask操作中，将描述行改写成：‘原有描述行’_|_A_to_G_|_‘被mask序列0101001的int转换’_|_read2
            fo.write(line1+'_|_'+mask_from+'_to_'+mask_to+'_|_'+str(int(record,2))+'_|_read2'+'\n')
            #将read中‘’碱基转换成‘’碱基
            fo.write(line2.replace(mask_from,mask_to)+'\n')            
            fo.write('+\n')
            fo.write(line4+'\n')
            line1=fi.readline().replace('\n','')
            line2=fi.readline().upper().replace('\n','')
            line3=fi.readline().replace('\n','')
            line4=fi.readline().replace('\n','')
        fi.close()
        fo.close()

    def poly_check(seq,poly_limit):
        if 'A'*poly_limit not in seq and 'T'*poly_limit not in seq and 'G'*poly_limit not in seq and 'C'*poly_limit not in seq:
            return True
        else:
            return False
                
    def var_check(change_from, change_to, seq,var_limit):
        ALL=['A','T','C','G']
        tmp=[]
        for one in ALL:
            if one !=change_from.upper() and one !=change_to.upper():
                tmp.append(one)
        flag=1
        for one in tmp:
            if seq.count(one) < var_limit/(float(len(tmp))+2):
                flag=0
        if flag==1:
            return True
        else:
            return False
    
    def reverse_base(base):
        base=base.upper()
        if base=='A':
            return 'T'
        elif base=='C':
            return 'G'
        elif base=='G':
            return 'C'
        elif base=='T':
            return 'A'
        else:
            return 'N'        

    def recover_sam(sam_in_dir,sam_out_dir, var_limit=20,poly_limit=10,rm_multi=0):
        fi=open(sam_in_dir)
        fo=open(sam_out_dir,'w')
        for line in fi:
            seq=line.split('\t')
            if line[0]=='@':
                fo.write(line)
            #去除sam文件转换flag==4的reads
            elif seq[1]=='4' and seq[2]=='*':
                break
            elif seq[1]!='4' and len(seq)>=9:
                seq=line.split('\t')
                seq[9]=seq[9].upper()
                seq[1]=int(seq[1])
                #reads flag值至少大于16
                if len(bin(seq[1]))>=7:
                    if bin(seq[1])[-3]!='1':
                        if bin(seq[1])[-5]=='1':
                            seq[1]='16'
                        else:
                            seq[1]='0'
                seq[1]=str(seq[1])
                #得到对应read的flag值、mask_from、mask_to
                #此步骤相当于maskfq.py操作的相反步骤
                record=bin(int(seq[0].split('_|_')[2]))[3:]
                change_from=seq[0].split('_|_')[1].split('_')[0]
                change_to=seq[0].split('_|_')[1].split('_')[2]
            
            
                if len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=='1':   #seq[1]=='16':
                    change_from=reverse_base(change_from)
                    change_to=reverse_base(change_to)
                    record=record[::-1]
                else:
                    record=record
                changed_read=seq[9]
                i=0
                recovered_read=''
                while i<len(seq[9]):
                    if record[i]=='1' and seq[9][i]==change_to:
                        recovered_read += change_from
                    elif record[i]=='1' and seq[9][i]!=change_to:
                        #print "recover error in "+seq[0]
                        recovered_read += 'N'
                    else:
                        recovered_read += seq[9][i]
                    i=i+1
                seq[9]=recovered_read
                #fo.write(seq[0])
                #if len(record)==len(seq[9]) and 'I' not in seq[5] and 'D' not in seq[5] and len(record)-changed_read.count(change_to) > 25 and poly_check(seq[9],poly_limit):
                #筛选C+T数量大于20、删除碱基N的数量>=10的reads
                if len(record)==len(seq[9]) and len(record)-changed_read.count(change_to) > var_limit and poly_check(seq[9],poly_limit): #and var_check(change_from,change_to,seq[9],var_limit):
                    if (rm_multi==1 and "XA:Z:" not in line) or rm_multi==0:
                        fo.write(seq[0])
                        j=1
                        while j<len(seq):
                            fo.write('\t'+seq[j])
                            j=j+1
            
            #        1+1
        fo.close()
        fi.close()
    
    def sam2zz(sam_in_dir=0,fa_in_dir=0,zz_out_dir=0):
        #-----------------------------------------------------------
        #reading refgenome
        fref=open(fa_in_dir)
        chrom={}
        chrr=''
        line=fref.read()
        line=line.split('>')
        for seq in line:
                if ' ' in seq:
                    chrr=seq[0:seq.find(' ')]
                else:
                    chrr=seq[0:seq.find('\n')]
                chrom[chrr]=seq[seq.find('\n'):].replace('\n','')

        #0 base

        fref.close()
        #------------------------------------------------------------
        #
        #donee:dolist;lst:donumlist
        def doCG(a):
            donee=[]
            lst=re.findall( '(\d+|\+|-|\*|/)', a )
            for i in a:
                if i == 'I' or i == 'D' or i== 'M' or i=='S' or i=='P' or i=='N' :
                    donee.append(i)
            return donee,lst
        #donefunction

        def doneCG(CG,chrr,pos,seq,qseq):#pos is 1 base
            donee,lst=doCG(CG)
            errorsite=''
            intersite=''
            quasite=''
            locsite=''
            pieceloc=''
            refseq=''
            seqseq=''
            refpos=int(pos)-1
            seqpos=0
            step=0
            while step<len(donee):
                if donee[step]=='I':
                    seqpos=seqpos+int(lst[step])
                elif donee[step]=='D':
                    refpos=refpos+int(lst[step])
                elif donee[step]=='N':
                    refpos=refpos+int(lst[step])
                elif donee[step]=='S':
                    seqpos=seqpos+int(lst[step])
                elif donee[step]=='M':
                    
                    refseq=refseq+chrom[chrr][refpos:refpos+int(lst[step])]
                    seqseq=seqseq+seq[seqpos:seqpos+int(lst[step])]
                    j=refpos
                    jj=seqpos
                    while j<refpos+int(lst[step]):
                        try:
                            if chrom[chrr][j].upper() != seq[jj].upper() and chrom[chrr][j].upper() !='N' and seq[jj].upper() != 'N':
                                errorsite=errorsite+chrom[chrr][j].upper()+seq[jj].upper()+':'+str(j+1)+';'
                                quasite=quasite+','+str(ord(qseq[jj]))
                                locsite=locsite+','+str(jj+1)
                                pieceloc=pieceloc+','+str( min( jj+1-seqpos,seqpos+int(lst[step])-jj  ) )
                                        
                        except Exception, e:
                            pass
                            #print "error with",chrr,pos,e
                        j=j+1
                        jj=jj+1
                    intersite=intersite+str(refpos+1)+':'+str(refpos+int(lst[step]))+';'
                    refpos=refpos+int(lst[step])
                    seqpos=seqpos+int(lst[step])
                step=step+1
            refseq=refseq.upper()
            seqseq=seqseq.upper()
            return refseq,seqseq,errorsite,intersite,quasite,locsite,pieceloc
        #------------------------------------------------------------
        ###################additional
        '''
        whole={}
        fi=open(sam_in_dir)
        for line in fi:
            seq=line.split('\t')
            if line[0]!='@'and len(seq)>5:
                if seq[0][0]!='@' and seq[2]!='*' and seq[5]!='*':
                    name=seq[0].split('_|_')[0]
                    
                    #tmp = [ seq[4].count(':') ,seq[0] ]
                    if name in whole:
                    #    if tmp[0] < whole[name][0]:
                    #        whole[name]=tmp
                        
                        whole[name] +=1
                    else:
                        whole[name]=1
                        #whole[name]=tmp
        fi.close()
        '''
        #######################            
        fi=open(sam_in_dir) #sam
        fo=open(zz_out_dir,'w') #zz
        for line in fi:
            seq=line.split('\t')
            if line[0]!='@' and len(seq)>5:
                name=seq[0].split('_|_')[0]
                if seq[0][0]!='@' and seq[2]!='*' and seq[5]!='*' :# and whole[name]<2:
                    refseq,seqseq,errorsite,intersite,quasite,locsite,pieceloc=doneCG(seq[5],seq[2],seq[3],seq[9],seq[10])
                    quasite=quasite[1:]
                    locsite=locsite[1:]
                    pieceloc=pieceloc[1:]
                    
                    if len(intersite[0:-1])==0:
                        intersite='*;'
                    if len(errorsite[0:-1])==0:
                        errorsite='*;'
                    if len(quasite)==0:
                        quasite='*'
                    if len(locsite)==0:
                        locsite='*'
                    if len(pieceloc)==0:
                        pieceloc='*'
                    
                    if len(bin(int(seq[1])))>=9 and bin(int(seq[1]))[-7]=='1':
                        fo.write(seq[2]+'\t'+seq[1]+'\t'+seq[4]+'\t'+intersite[0:-1]+'\t'+errorsite[0:-1]+'\t'+quasite+'\t'+locsite+'\t'+seq[9]+'\t'+seq[0]+'_1'+'\t'+pieceloc+'\n')
                    elif len(bin(int(seq[1])))>=10 and bin(int(seq[1]))[-8]=='1':
                        fo.write(seq[2]+'\t'+seq[1]+'\t'+seq[4]+'\t'+intersite[0:-1]+'\t'+errorsite[0:-1]+'\t'+quasite+'\t'+locsite+'\t'+seq[9]+'\t'+seq[0]+'_2'+'\t'+pieceloc+'\n')
                    else:
                        fo.write(seq[2]+'\t'+seq[1]+'\t'+seq[4]+'\t'+intersite[0:-1]+'\t'+errorsite[0:-1]+'\t'+quasite+'\t'+locsite+'\t'+seq[9]+'\t'+seq[0]+'\t'+pieceloc+'\n')
                        
        fi.close()    
        fo.close()
      
    def dedup(in_dir,out_dir):
        old=set()
        fi=open(in_dir)
        fo=open(out_dir,'w')    
        for line in fi:
            seq=line.split('\t')
            if seq[0]+':'+seq[3]+':'+seq[7] not in old:
                fo.write(line)
                old.add(seq[0]+':'+seq[3]+':'+seq[7])

    def mismatch_num(seqlen):
        if seqlen < 17:
            return 1
        elif seqlen <38:
            return 2
        elif seqlen <64:
            return 3
        elif seqlen <93:
            return 4
        elif seqlen <124:
            return 5
        elif seqlen <157:
            return 6
        elif seqlen <190:
            return 7
        elif seqlen <225:
            return 8
        else:
            return 9


    def mask_zz2snv(zz_in_dir=0,bed_out_dir=0,baseq_cutoff_dir=0):

        fi=open(zz_in_dir)
        fo=open(bed_out_dir,'w')

        fqua=open(baseq_cutoff_dir)
        #reads碱基质量阈值
        limitbasequa=int(fqua.readline().replace('\n',''))
        fqua.close()
        limitad=1
        limitloc=5
        #读段比对的质量
        limitmpqua=0
        allsnv={}
        for line in fi:
            truesnv=[]
            seq=line.rstrip().split('\t')
            mismatch=seq[4].split(';')
            basequa=seq[5].split(',')
            loc=seq[9].split(',') #fragment-loc
            mpqua=int(seq[2])
            ####################################################################
            #change the sam flag 'seq[1]' when you didn't use "hisat2 -aln" as mapper
            seq[1]=int(seq[1])
            if len(bin(seq[1]))>=7:
                if bin(seq[1])[-3]!='1':
                    if bin(seq[1])[-5]=='1':
                        seq[1]='16'
                    else:
                        seq[1]='0'
            elif len(bin(seq[1]))>=5:
                if bin(seq[1])[-3]!='1':
                    seq[1]='0'
                                
            else: 
                seq[1]='0'
            #####################################################################
        
            
            
            if basequa[0]!='*' and mpqua >= limitmpqua and mpqua < 200:
                i=0
                baseqlst=[]
                mistype={}
                while i < len(basequa):
                    baseqlst.append(int(basequa[i]))
                    try:
                        mistype[mismatch[i].split(':')[0]] += 1
                    except Exception , e:
                        mistype[mismatch[i].split(':')[0]] = 1
                    if int(basequa[i]) >= limitbasequa  and  int(loc[i]) > limitloc  :

                        truesnv.append([mismatch[i],seq[1],seq[8]])
                        
                    i=i+1

                #fflag=1
                #masktype_tmp=seq[8].split('_|_')[1].split('_to_')
                #masktype=masktype_tmp[0]+masktype_tmp[1]
                miss=[]
                for mis in mistype:
                    miss.append(mistype[mis])
                miss.sort()
                if len(miss)>=2:
                    missnum=sum(miss[:-1])
                else:
                    missnum=0
                


                #snv平均质量大于SPRINT初始化时得到的碱基质量（58/...）
                if len(baseqlst)>0 and sum(baseqlst)/float(len(baseqlst)) >= limitbasequa: #and missnum <= mismatch_num(len(seq[7])):
                    for snv in truesnv:
                        #snv：mismatch+flag+read_barcode
                        try:
                            #allsnv：以chr+mismatch类型+位点+read为index，其对应的该snv的数量为值
                            allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line][0]=allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line][0]+1
                            if len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=='1':
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line][1]=allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line][1]+1
                            elif len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-len(bin(int(seq[1])))]=='0':
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line][2]=allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line][2]+1
                        except Exception, e:
                            if zz_in_dir==output+"/tmp/genome_all.zz.dedup":
                                add_line=seq[8]
                            else:
                                add_line=seq[8].split("_|_")[0]
                                #read_name_2=read_name_1.split('_')[2]
                                #add_line=read_name_2
                            if len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=='1':
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line]=[1,1,0]
                            elif len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=='0':
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]+'\t'+add_line]=[1,0,1]

        #snv：chr1 3055665 CT TGGCGCAGTATCGCAT
        snv_bed=[]
        for snv in allsnv:
            #print snv
            seq=snv.split('\t')
            #若某snv位点覆盖深度>0
            if allsnv[snv][0]>=limitad:
                    #决定该snv是位于正链或负链
                    if allsnv[snv][1] > allsnv[snv][2]: 
                        snv_bed.append([seq[0],int(seq[2]),seq[1],'+',allsnv[snv][0],seq[3]])
                    elif allsnv[snv][2] > allsnv[snv][1]:
                        snv_bed.append([seq[0],int(seq[2]),seq[1],'-',allsnv[snv][0],seq[3]])
                    else:
                        snv_bed.append([seq[0],int(seq[2]),seq[1],'.',allsnv[snv][0],seq[3]])

        snv_bed.sort()
        for one in snv_bed:
            #one[0]:chr;one[1]:snv位点；one[2]:mismatch类型；one[3]:正义链or反义链；one[4]:snv数量；one[5]：read_barcode
            fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+str(one[4])+'\t'+one[3]+'\t'+one[5]+'\n')
        fi.close()
        fo.close()

    def tzz2gzz(trans_loc_file_in , transcript_zz_in , genome_zz_out):
        fa=open(trans_loc_file_in)
        l1=fa.readline()
        l2=fa.readline()
        TRANS={}
        while l1 !='':
            trans=l1[1:].rstrip()
            seq=l2.split(';')[:-1]
            TRANS[trans]=seq
            l1=fa.readline()
            l2=fa.readline()

        fa.close()

        def loc_t2g(tCHR,tLOC):
            gCHR=tCHR.split('_|_')[1]
            seq=TRANS[tCHR]
            tLOC=int(tLOC)
            flag=1
            tmp=0
            i=0
            while i<len(seq) and flag==1:
                end=int(seq[i].split(',')[1])
                start=int(seq[i].split(',')[0])
                tmp +=  end-start+1
                if tmp >= tLOC :
                    flag=0
                i += 1
            j=i-1
            dis2end  = tmp-tLOC
            gLOC = int(seq[j].split(',')[1]) - dis2end
            return gLOC

        def range_t2g(tCHR,tstart,tend):
            gCHR=tCHR.split('_|_')[1]
            tstart=int(tstart)
            tend=int(tend)
            tloc = [tstart+i for i in range(tend-tstart+1)]
            gloc = []
            for one in tloc:
                gloc.append(loc_t2g(tCHR,one))
            i=1
            record=[gloc[0]]
            out=[]
            while i < len(gloc):
                if abs(gloc[i]-gloc[i-1])> 1:
                    out.append(str(record[0])+':'+str(record[-1]))
                    record=[gloc[i]]
                else:
                    record.append(gloc[i])
                i+=1
            out.append(str(record[0])+':'+str(record[-1]))
            return out

        fi=open(transcript_zz_in)
        fo=open(genome_zz_out,'w')
        for line in fi:
            seq=line.split('\t')
            tCHR = seq[0]
            gCHR=tCHR.split('_|_')[1]
            if seq[4]!='*':
                snv=seq[4].split(';')
                gsnv=[]
                for one in snv:
                    out=loc_t2g(tCHR,one.split(':')[1])
                    gsnv.append(one.split(':')[0]+':'+str(out))
                seq[4]=';'.join(gsnv)

            trange = seq[3].split(';')
            grange=[]
            for one in trange:
                new_one = range_t2g(tCHR, one.split(':')[0],one.split(':')[1])
                grange += new_one
            seq[3] = ';'.join(grange)
            seq[0]=gCHR
            fo.write('\t'.join(seq))

    def sort_zz(zz_in,zz_out):
        out=[]
        fi=open(zz_in)
        for line in fi:
            seq=line.split('\t')
            out.append([seq[0],int(seq[3].split(':')[0]),int(seq[3].split(':')[-1]),line])
        fi.close()
        fo=open(zz_out, 'w')
        out.sort()
        for one in out:
            fo.write(one[3])

    def transcript_locator(fbed_in_dir=0,ftransloc_in_dir=0,fbed_out_dir=0):
        if fbed_in_dir==0 or ftransloc_in_dir==0:
            print 'fbed_in_dir\tftransloc_in_dir\tfbed_out_dir'
            return 0

        ftransloc=open(ftransloc_in_dir)
        fi=open(fbed_in_dir)
        fo=open(fbed_out_dir,'w')


        transcript={}
        whole=ftransloc.read().split('>')[1:]
        for one in whole:
            line=one.split('\n')
            transcript[line[0]]=line[1]
            #print line[0]


        for line in fi:
            seq=line.rstrip().split('\t')
            #print seq[0]
            if seq[0] in transcript:
                transcript_id=seq[0]
                chrr=transcript_id.split('_|_')[1]
                trans_loc=int(seq[2])
                loc=[]
                bande=transcript[seq[0]].split(';')[:-1]
                for one in bande:
                    be=one.split(',')
                    begin=int(be[0])
                    end=int(be[1])
                    tmp = [begin]*(end-begin+1)
                    i=0
                    while i< len(tmp):
                        tmp[i] = tmp[i]+i
                        i +=1
                    loc += tmp
                ref_loc=loc[trans_loc-1]
                fo.write(chrr+'\t'+str(ref_loc-1)+'\t'+str(ref_loc)+'\t'+'\t'.join(seq[3:])+'\n')#'\t'+chrom[chrr][ref_loc-1]+':'+chrom_t[transcript_id][trans_loc-1]+'\n')

    def transcript_sort(fbed_in_dir=0,fbed_out_dir=0):
        fi=open(fbed_in_dir)
        fo=open(fbed_out_dir,'w')

        whole={}
        for line in fi:
            seq=line.rstrip().split('\t')
            try:
                whole[seq[0]+':'+seq[2]+':'+seq[3]+':'+seq[5]+':'+seq[6]] += int(seq[4])
            except Exception,e:
                whole[seq[0]+':'+seq[2]+':'+seq[3]+':'+seq[5]+':'+seq[6]] = int(seq[4])

        tmp=[]
        for one in whole:
            seq=one.split(':')
            dep=str(whole[one])
            tmp.append([seq[0],int(seq[1]),seq[2],dep,seq[3],seq[4]])
        tmp.sort()
        for one in tmp:
            print one
            fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\t'+one[4]+'\t'+one[5]+'\n')
    
    def snv_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
        f1=open(bed_in_dir1)
        f2=open(bed_in_dir2)
        fo=open(bed_out_dir,'w')
        whole={}
        for line in f1:
            seq=line.rstrip().split('\t')
            whole[':'.join(seq[0:4])+':'+seq[5]+':'+seq[6]]=int(seq[4])
        for line in f2:
            seq=line.rstrip().split('\t')
            if ':'.join(seq[0:4])+':'+seq[5]+':'+seq[6] in whole:
                whole[':'.join(seq[0:4])+':'+seq[5]+':'+seq[6]] +=int(seq[4])
            else:
                whole[':'.join(seq[0:4])+':'+seq[5]+':'+seq[6]]=int(seq[4])
        lst=[]
        for one in whole:
            seq=one.split(':')
            lst.append([seq[0],int(seq[1]),int(seq[2]),seq[3],str(whole[one]),seq[4],seq[5]])
        lst.sort()
        for one in lst:
            out=[]
            for i in one:
                out.append(str(i))
            fo.write('\t'.join(out)+'\n')

        f1.close()
        f2.close()
        fo.close()

    #use .bed to annotate .bed
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

    def snv_cluster(bed_in_dir=0, bed_out_dir=0, cluster_distance=-1, cluster_size=-1):
        fi = open(bed_in_dir)
        tmp = 'chr0:0:AA'
        limitdistance = int(cluster_distance)
        limitnum = int(cluster_size)
        lst = []
        barcode = None
        fo = open(bed_out_dir, 'w')
        
        for line in fi:
            seq = line.split('\t')
            tmpseq = tmp.split(':')
            
            # 解析第7列读段信息中的barcode
            current_barcode = seq[6].split('_')[2]  # 假设barcode格式为id_xxxxxx_xxx_xxx形式
            
            if barcode is None:
                barcode = current_barcode
            elif barcode != current_barcode:
                sites_in_cluster = []
                for tmp_site in lst:
                    sites_in_cluster.append(tmp_site.split('\t')[2])
                if len(np.unique(sites_in_cluster)) >= limitnum:
                    begin = float(lst[0].split('\t')[1])
                    end = float(lst[-1].split('\t')[2])
                    if int(end - begin) != 0:
                        density = len(np.unique(sites_in_cluster)) / (end - begin)
                        for one in lst:
                            fo.write(one[:-1] + '\t' + str(len(np.unique(sites_in_cluster))) + '\t' + str(density) + '\n')
                lst = [line]
                barcode = current_barcode
                continue
            
            # 检查seq是否满足条件
            if seq[0] == tmpseq[0] and int(seq[2]) - int(tmpseq[1]) <= limitdistance and seq[3] == tmpseq[2]:
                lst.append(line)
            else:
                sites_in_cluster = []
                for tmp_site in lst:
                    sites_in_cluster.append(tmp_site.split('\t')[2])
                if len(np.unique(sites_in_cluster)) >= limitnum:
                    begin = float(lst[0].split('\t')[1])
                    end = float(lst[-1].split('\t')[2])
                    if int(end - begin) != 0:
                        density = len(np.unique(sites_in_cluster)) / (end - begin)
                        for one in lst:
                            fo.write(one[:-1] + '\t' + str(len(np.unique(sites_in_cluster))) + '\t' + str(density) + '\n')
                lst = [line]
                tmp = seq[0] + ':' + seq[2] + ':' + seq[3]
        
        # 处理最后一组满足条件的序列
        sites_in_cluster = []
        for tmp_site in lst:
            sites_in_cluster.append(tmp_site.split('\t')[2])
        if len(np.unique(sites_in_cluster)) >= limitnum:
            begin = float(lst[0].split('\t')[1])
            end = float(lst[-1].split('\t')[2])
            if int(end - begin) != 0:
                density = len(np.unique(sites_in_cluster)) / (end - begin)
                for one in lst:
                    fo.write(one[:-1] + '\t' + str(len(np.unique(sites_in_cluster))) + '\t' + str(density) + '\n')
        
        fi.close()
        fo.close()


        """ fo=open(bed_out_dir,'w')
        #final_dict=[]
        for i in range(0,len(tmp_dict)): 
            cluster_barcodes=[]  
            cluster_sites=[]         
            for j in range(0,len(tmp_dict[i])):
                each_line=tmp_dict[i][j]
                each_seq=each_line.rstrip().split('\t')
                each_chr_site=each_seq[0]+'_'+each_seq[2]
                each_barcode=each_seq[6]
                cluster_barcodes.append(each_barcode)
                cluster_sites.append(each_chr_site)
            #最后一列的数字为该簇所含编辑位点的数量
            if len(np.unique(cluster_barcodes))==1:
                for one in tmp_dict[i]:
                    fo.write(one[0:-1]+'\t'+str(len(tmp_dict[i]))+'\n')    
            else:
                #所有cluster严格控制在一个细胞里！！！！！！！！！
                for one_barcode in cluster_barcodes:
                    len_cluster_sites=cluster_barcodes.count(one_barcode)
                    index_one=[]
                    # 获取第一个barcode的下标
                    loc=cluster_barcodes.index(one_barcode)
                    index_one.append(loc)
                    while loc < len(cluster_barcodes)-1:
                        # 从第一个barcode的下一个位置开始查找, 所以加1
                        try:
                            index_one.append(cluster_barcodes.index(one_barcode, loc + 1))
                            loc=cluster_barcodes.index(one_barcode, loc + 1)+1
                        except ValueError as e:
                            break
                    if len_cluster_sites > limitnum: 
                        for k in np.unique(index_one):                       
                            one=tmp_dict[i][k]
                            fo.write(one[0:-1]+'\t'+str(len_cluster_sites)+'\n')                 
        fo.close() """
    
    def bed_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
        f1=open(bed_in_dir1)
        f2=open(bed_in_dir2)
        fo=open(bed_out_dir,'w')
        whole=[]
        for line in f1:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5],seq[6]])
        for line in f2:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5],seq[6]])
        whole.sort()
        old=set()
        for one in whole:
            if one[0]+':'+str(one[1])+':'+one[2]+':'+one[5] not in old:
                fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\t'+one[4]+'\t'+one[5]+'\n')
                old.add(one[0]+':'+str(one[1])+':'+one[2]+':'+one[5])
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
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5],seq[6]])
        for line in f2:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5],seq[6]])
        for line in f3:
            seq=line.replace('\n','').split('\t')
            whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5],seq[6]])
        whole.sort()
        for one in whole:
            fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\t'+one[4]+'\t'+one[5]+'\n')
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
            whole[':'.join(seq[0:4])+':'+seq[5]+':'+seq[6]]=int(seq[4])
        for line in f2:
            seq=line.rstrip().split('\t')
            if ':'.join(seq[0:4])+':'+seq[5]+':'+seq[6] in whole:
                pass
                #whole[':'.join(seq[0:4])+':'+seq[5]] +=int(seq[4])
            else:
                whole[':'.join(seq[0:4])+':'+seq[5]+':'+seq[6]]=int(seq[4])
        lst=[]
        for one in whole:
            seq=one.split(':')
            lst.append([seq[0],int(seq[1]),int(seq[2]),seq[3],str(whole[one]),seq[4],seq[5]])
        lst.sort()
        for one in lst:
            out=[]
            for i in one:
                out.append(str(i))
            fo.write('\t'.join(out)+'\n')

        f1.close()
        f2.close()
        fo.close()

    def get_res_in_1cell(bed_in_dir, bed_out_dir):
        # 读取输入文件并按barcode分组SNV，同时将结果存储到一个列表中
        results = []
        barcode_lst = {}
        with open(bed_in_dir) as fi:
            for line in fi:
                seq = line.rstrip().split('\t')
                barcode = seq[6].split('_')[2]
                barcode_lst.setdefault(barcode, []).append('_'.join(seq[0:4])+'_'+seq[5])

        # 处理每个barcode的SNV列表，计算计数并添加到results列表中
        for barcode, snv_list in barcode_lst.items():
            counter = Counter(snv_list)
            for snv, count in counter.items():
                chr_site=snv.split('_')[0]
                site_0=snv.split('_')[1]
                site_1=snv.split('_')[2]
                type=snv.split('_')[3]
                strand = snv.split('_')[4]
                if strand in ['+', '-', '.']:
                    results.append("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chr_site, site_0, site_1, type, count, strand, barcode))

        # 对results列表进行排序
        sorted_results = sorted(results, key=lambda x: (x.split()[0], int(x.split()[1]), int(x.split()[2])))

        # 将排序后的结果写入输出文件
        with open(bed_out_dir, 'w') as fo:
            fo.writelines([result + '\n' for result in sorted_results])

    def o2b(bed_in,bed_out):
        fi=open(bed_in)
        fo=open(bed_out,'w')
        for line in fi:
            seq=line.rstrip().split('\t')
            #if float(seq[7])<0.2:
            fo.write('\t'.join(seq[0:6])+'\n')

    def get_depth(zz_in_dir=0,bed_in_dir=0,bed_out_dir=0):
        fread=open(zz_in_dir)# './zz_folder/all.zz')
        fsnv=open(bed_in_dir) #'../bed_folder/try_new.bed') #   Hap_SRR521447.bed')
        fo=open(bed_out_dir,'w')#'./tmp/readspersite_new.zer','w')

        class Read:
            def __init__(self,read):
                #snv初始化：例如——'AG:3000287 GT:3000325'
                self.snv=read.split('\t')[4].split(';')
                #snv区间（位点）初始化：例如——'3388959:3389056'
                self.inter=read.split('\t')[3].split(';')
                #该line的flag  
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
            #getmin与getmax取该line区间的最小值与最大值
            def getmin(self):
                return int(self.inter[0].split(':')[0])
            def getmax(self):
                return int(self.inter[-1].split(':')[1])

        #以染色体为index，某染色体（例如chr1）zz的所有line为element
        reads={}
        for line in fread:
            seq=line.split('\t')
            if '_|_' in seq[8]:
                barcode=seq[8].split("_|_")[0].split('_')[2]
            else:
                barcode=seq[8].split("_")[2]
            try:
                    reads[seq[0]+'\t'+barcode].append(line[0:-1])
            except Exception,e :
                    #print seq[0]+' begin'
                    reads[seq[0]+'\t'+barcode]=[line[0:-1]]


        top=0
        chrr=''
        barcode=''
        for line in fsnv:
            seq=line.rstrip().split('\t')
            deep=0
            altdeep=0
            try:
                #snv：mismatch+site   
                snv=seq[3]+':'+seq[2]
                if seq[0] != chrr or seq[6]!=barcode:
                    top=0
                    chrr=seq[0]
                    barcode=seq[6]
                if seq[0]+'\t'+barcode not in reads:
                    reads[seq[0]+'\t'+barcode]=[]
                #若top的值小于以chr+barcode为index的element数量，则进行遍历
                #top的作用：迭代某染色体前，对某一条zz中的line做初始化
                if top < len(reads[seq[0]+'\t'+barcode]):
                    while seq[0]==chrr and seq[6]==barcode and top < len(reads[seq[0]+'\t'+barcode]) and Read(reads[seq[0]+'\t'+barcode][top]).getmax() < int(seq[2]):
                        top=top+1
                    point=top
                    while seq[0]==chrr and seq[6]==barcode and point < len(reads[seq[0]+'\t'+barcode]) and Read(reads[seq[0]+'\t'+barcode][point]).getmin() <= int(seq[2]):
                        if Read(reads[seq[0]+'\t'+barcode][point]).locisin(seq[2]) ==1:
                            deep=deep+1
                        
                        if Read(reads[seq[0]+'\t'+barcode][point]).snvisin(snv)==1:
                            altdeep=altdeep+1
                            
                    
                        point=point+1
                fo.write(line[0:-1]+'\t'+str(altdeep)+':'+str(deep)+'\n')
            except Exception, e:
                print line
        fread.close()
        fsnv.close()
        fo.close()


    #try:
    if 1==1: 

        #预处理阶段
        """ print 'preprocessing...'
        # 此处为针对双端测序reads的整理方式
        #cut(read1,read_r1,barcode_ref,tmp+'cut_read1.fastq',cutbp,cutbp_r1,'read1')
        #extract_reads(barcode_ref, read1, read2, tmp+'cut_read1.fastq', tmp+'cut_read2.fastq', cutbp_r1)
        
        get_baseq_cutoff(read1,tmp+'baseq.cutoff')
                    
        print 'mapping...'
    
        #step1_1：将raw fastq文件经处理后转换成bam文件
        TAG='genome'
        fq2sam(TAG,paired_end,tmp+'cut_read1.fastq',tmp+'cut_read2.fastq',tmp,refgenome,hisat2,samtools,mapcpu,read_format)

        #step1_2
        #将未匹配到基因组上的reads从bam文件中提取出来、并转换为sam格式(包含了read1和read2总共的unmapped reads)
        TAG='genome'
        subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
        #将上述sam文件转换为fq文件
        umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/'+TAG+'_unmapped.fq') 
     
        #若存在转录参考基因组，则将unmapped reads和转录组模版再比对一次
        if os.path.exists(refgenome+'.trans.fa'):
            TAG='transcript'
            fq2sam(TAG,False,tmp+'/genome_unmapped.fq',read2,tmp,refgenome+'.trans.fa',hisat2,samtools,mapcpu,read_format)
    
            subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
            umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/regular_unmapped.fq')    
            maskfq(tmp+'/regular_unmapped.fq','A','G')
        else:
            umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/regular_unmapped.fq')
            maskfq(tmp+'/regular_unmapped.fq','A','G') 
    
        #step2_1:
        #将step1中未匹配到基因组上的、经A-G转换的read加上read2一起，再在经A-G转换的参考基因组上比对一次         
        TAG='genome_mskAG'
        fq2sam(TAG,False,tmp+'/regular_unmapped_A_to_G.fq',read2,tmp,refgenome+'.mskAG.fa',hisat2,samtools,mapcpu,read_format)
    
        subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
        umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/'+TAG+'_unmapped.fq')    
    
        #step2_2:
        # 将step1中未匹配到基因组上的、经A-G转换的read加上read2一起，再在经T-C转换的参考基因组上比对一次
        TAG='genome_mskTC'
        fq2sam(TAG,False,tmp+'/regular_unmapped_A_to_G.fq',read2,tmp,refgenome+'.mskTC.fa',hisat2,samtools,mapcpu,read_format)
    
        subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
        umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/'+TAG+'_unmapped.fq')    
        
        #在有转录参考基因组样本的情况下，经A-G转换的reads与经A-G、T-C转换的转录参考基因组再匹配一次
        if os.path.exists(refgenome+'.trans.fa'):
            TAG='transcript_mskAG'
            fq2sam(TAG,False,tmp+'/genome_mskAG_unmapped.fq',read2,tmp,refgenome+'.trans.fa.mskAG.fa',hisat2,samtools,mapcpu,read_format)
    
            TAG='transcript_mskTC'
            fq2sam(TAG,False,tmp+'/genome_mskTC_unmapped.fq',read2,tmp,refgenome+'.trans.fa.mskTC.fa',hisat2,samtools,mapcpu,read_format)  
    
        if os.path.exists(tmp+'genome_mskAG_unmapped.sam'):
                    if os.path.exists(tmp+'cut_read1.fastq'):
                            os.remove(tmp+'cut_read1.fastq')
                    if os.path.exists(tmp+'cut_read2.fastq'):
                            os.remove(tmp+'cut_read2.fastq')
                    if os.path.exists(tmp+'genome_mskAG_unmapped.fq'):
                            os.remove(tmp+'genome_mskAG_unmapped.fq')
                    if os.path.exists(tmp+'genome_mskAG_unmapped.sam'):
                            os.remove(tmp+'genome_mskAG_unmapped.sam')
                    if os.path.exists(tmp+'genome_mskTC_unmapped.fq'):
                            os.remove(tmp+'genome_mskTC_unmapped.fq')
                    if os.path.exists(tmp+'genome_mskTC_unmapped.sam'):
                            os.remove(tmp+'genome_mskTC_unmapped.sam')
                    if os.path.exists(tmp+'genome_unmapped.fq'):
                            os.remove(tmp+'genome_unmapped.fq')
                    if os.path.exists(tmp+'genome_unmapped.sam'):
                            os.remove(tmp+'genome_unmapped.sam')
                    if os.path.exists(tmp+'transcript_unmapped_A_to_G.fq'):
                            os.remove(tmp+'transcript_unmapped_A_to_G.fq')
                    if os.path.exists(tmp+'transcript_unmapped.fq'):
                            os.remove(tmp+'transcript_unmapped.fq')
                    if os.path.exists(tmp+'transcript_unmapped.sam'):
                            os.remove(tmp+'transcript_unmapped.sam')
                    if os.path.exists(tmp+'regular_unmapped.fq'):
                            os.remove(tmp+'regular_unmapped.fq')
                    if os.path.exists(tmp+'regular_unmapped_A_to_G.fq'):
                            os.remove(tmp+'regular_unmapped_A_to_G.fq') 
                            
        #step3：
        #将A-G/T-C转换的sam文件还原，其中根据poly尾、CT含量对mapped reads进行过滤。且将其转换成zz格式(snv结果文件)
        recover_sam(tmp+'genome_mskAG_all.sam',tmp+'genome_mskAG_all.sam.rcv', var_limit, poly_limit, rm_multi) 
        sam2zz(tmp+'genome_mskAG_all.sam.rcv',refgenome,tmp+'genome_mskAG_all.zz')
        recover_sam(tmp+'genome_mskTC_all.sam',tmp+'genome_mskTC_all.sam.rcv', var_limit, poly_limit, rm_multi)
        sam2zz(tmp+'genome_mskTC_all.sam.rcv',refgenome,tmp+'genome_mskTC_all.zz') 
        sam2zz(tmp+'genome_all.sam',refgenome,tmp+'genome_all.zz')

        if os.path.exists(tmp+'genome_mskAG_all.sam.rcv'):
                os.remove(tmp+'genome_mskAG_all.sam.rcv')
        if os.path.exists(tmp+'genome_mskAG_all.sam'):
                os.remove(tmp+'genome_mskAG_all.sam')
        if os.path.exists(tmp+'genome_mskTC_all.sam.rcv'):
                os.remove(tmp+'genome_mskTC_all.sam.rcv')
        if os.path.exists(tmp+'genome_mskTC_all.sam'):
                os.remove(tmp+'genome_mskTC_all.sam')
        if os.path.exists(tmp+'genome_all.sam'):
                os.remove(tmp+'genome_all.sam')

        #同样考虑trans.fa存在的情况
        if os.path.exists(refgenome+'.trans.fa'):
            recover_sam(tmp+'transcript_mskAG_all.sam',tmp+'transcript_mskAG_all.sam.rcv', var_limit, poly_limit, rm_multi)
            sam2zz(tmp+'transcript_mskAG_all.sam.rcv',refgenome+'.trans.fa',tmp+'transcript_mskAG_all.zz')
            recover_sam(tmp+'transcript_mskTC_all.sam',tmp+'transcript_mskTC_all.sam.rcv', var_limit, poly_limit, rm_multi)
            sam2zz(tmp+'transcript_mskTC_all.sam.rcv',refgenome+'.trans.fa',tmp+'transcript_mskTC_all.zz')
            sam2zz(tmp+'transcript_all.sam',refgenome+'.trans.fa',tmp+'transcript_all.zz')
            
            if os.path.exists(tmp+'transcript_mskAG_all.sam.rcv'):
                os.remove(tmp+'transcript_mskAG_all.sam.rcv')
            if os.path.exists(tmp+'transcript_mskAG_all.sam'):
                os.remove(tmp+'transcript_mskAG_all.sam')
            if os.path.exists(tmp+'transcript_mskTC_all.sam.rcv'):
                os.remove(tmp+'transcript_mskTC_all.sam.rcv')
            if os.path.exists(tmp+'transcript_mskTC_all.sam'):
                os.remove(tmp+'transcript_mskTC_all.sam')
            if os.path.exists(tmp+'transcript_all.sam'):
                os.remove(tmp+'transcript_all.sam')

        
        #step4：       
        #zz文件去重(此处去重为PCR重复的去除)
        if os.path.exists(refgenome+'.trans.fa'):
            dedup(tmp+'transcript_mskAG_all.zz',tmp+'transcript_mskAG_all.zz.dedup')
            dedup(tmp+'transcript_mskTC_all.zz',tmp+'transcript_mskTC_all.zz.dedup') 
            dedup(tmp+'transcript_all.zz',tmp+'transcript_all.zz.dedup') 
    
        #基因组进行同样的操作，且最终仅留下去重的结果
        dedup(tmp+'genome_mskAG_all.zz',tmp+'genome_mskAG_all.zz.dedup') 
        dedup(tmp+'genome_mskTC_all.zz',tmp+'genome_mskTC_all.zz.dedup')
        dedup(tmp+'genome_all.zz',tmp+'genome_all.zz.dedup')
         
        if os.path.exists(tmp+'transcript_mskAG_all.zz'):
                os.remove(tmp+'transcript_mskAG_all.zz')
        if os.path.exists(tmp+'transcript_mskTC_all.zz'):
                os.remove(tmp+'transcript_mskTC_all.zz')
        if os.path.exists(tmp+'transcript_all.zz'):
                os.remove(tmp+'transcript_all.zz')
        if os.path.exists(tmp+'genome_mskAG_all.zz'):
                os.remove(tmp+'genome_mskAG_all.zz')
        if os.path.exists(tmp+'genome_mskTC_all.zz'):
                os.remove(tmp+'genome_mskTC_all.zz')
        if os.path.exists(tmp+'genome_all.zz'):
                os.remove(tmp+'genome_all.zz')

        print 'identifying SNVs...'
   
        if os.path.exists(refgenome+'.trans.fa'):
            #筛选符合条件的snv
            mask_zz2snv(tmp+'transcript_mskAG_all.zz.dedup',tmp+'transcript_mskAG_all.zz.dedup.snv',tmp+'baseq.cutoff') 
            mask_zz2snv(tmp+'transcript_mskTC_all.zz.dedup',tmp+'transcript_mskTC_all.zz.dedup.snv',tmp+'baseq.cutoff') 
            mask_zz2snv(tmp+'transcript_all.zz.dedup',tmp+'transcript_all.zz.dedup.snv',tmp+'baseq.cutoff')

            #将候选snv在转录本模板的定位转换成基因组上的位置
            tzz2gzz(refgenome+'.trans.fa.loc', tmp+'transcript_mskAG_all.zz.dedup', tmp+'transcript_mskAG_all.zz.dedup.genome.zz')
            tzz2gzz(refgenome+'.trans.fa.loc', tmp+'transcript_mskTC_all.zz.dedup', tmp+'transcript_mskTC_all.zz.dedup.genome.zz')
            tzz2gzz(refgenome+'.trans.fa.loc', tmp+'transcript_all.zz.dedup', tmp+'transcript_all.zz.dedup.genome.zz')
             
        #根据前述决定的碱基质量阈值对snv进行筛选并去重
        mask_zz2snv(tmp+'genome_mskAG_all.zz.dedup',tmp+'genome_mskAG_all.zz.dedup.snv',tmp+'baseq.cutoff') 
        mask_zz2snv(tmp+'genome_mskTC_all.zz.dedup',tmp+'genome_mskTC_all.zz.dedup.snv',tmp+'baseq.cutoff')
        mask_zz2snv(tmp+'genome_all.zz.dedup',tmp+'genome_all.zz.dedup.snv',tmp+'baseq.cutoff')
        
        
        #将mapped reads与unmapped reads两条线找到的snv候选合并并排序——all_combined.zz.sorted
        #all_combined.zz.sorted用于后续计算找到的所有RES的深度
        if os.path.exists(refgenome+'.trans.fa'):
            subprocess.Popen('cat '+tmp+'/genome_mskAG_all.zz.dedup '+tmp+'/genome_mskTC_all.zz.dedup '+tmp+'/genome_all.zz.dedup '+tmp+'/transcript_mskAG_all.zz.dedup.genome.zz '+tmp+'/transcript_mskTC_all.zz.dedup.genome.zz '+tmp+'/transcript_all.zz.dedup.genome.zz '+' > '+tmp+'/all_combined.zz',shell=True).wait()
            sort_zz(tmp+'/all_combined.zz', tmp+'/all_combined.zz.sorted')
        else: 
            subprocess.Popen('cat '+tmp+'/genome_mskAG_all.zz.dedup '+tmp+'/genome_mskTC_all.zz.dedup '+tmp+'/genome_all.zz.dedup '+' > '+tmp+'/all_combined.zz',shell=True).wait()
            sort_zz(tmp+'/all_combined.zz', tmp+'/all_combined.zz.sorted') 
        
     
        if os.path.exists(refgenome+'.trans.fa'):
            #将snv在转录本模板的定位转换成基因组上的位置
            transcript_locator(tmp+'transcript_mskAG_all.zz.dedup.snv',refgenome+'.trans.fa.loc', tmp+'transcript_mskAG_all.zz.dedup.snv.genome.snv')
            transcript_locator(tmp+'transcript_mskTC_all.zz.dedup.snv',refgenome+'.trans.fa.loc', tmp+'transcript_mskTC_all.zz.dedup.snv.genome.snv')
            transcript_locator(tmp+'transcript_all.zz.dedup.snv',refgenome+'.trans.fa.loc', tmp+'transcript_all.zz.dedup.snv.genome.snv')

            transcript_sort(tmp+'transcript_all.zz.dedup.snv.genome.snv',tmp+'transcript_all.zz.dedup.snv.genome.snv.sort')
            transcript_sort(tmp+'transcript_mskTC_all.zz.dedup.snv.genome.snv',tmp+'transcript_mskTC_all.zz.dedup.snv.genome.snv.sort')
            transcript_sort(tmp+'transcript_mskAG_all.zz.dedup.snv.genome.snv',tmp+'transcript_mskAG_all.zz.dedup.snv.genome.snv.sort')
            
            #将mapped reads与unmapped reads两条线找到的snv合并
            snv_or(tmp+'transcript_all.zz.dedup.snv.genome.snv.sort',tmp+'genome_all.zz.dedup.snv',tmp+'regular.snv')
            snv_or(tmp+'transcript_mskTC_all.zz.dedup.snv.genome.snv.sort',tmp+'genome_mskTC_all.zz.dedup.snv', tmp+'hyper_mskTC.snv')
            snv_or(tmp+'transcript_mskAG_all.zz.dedup.snv.genome.snv.sort',tmp+'genome_mskAG_all.zz.dedup.snv', tmp+'hyper_mskAG.snv')


        else:
            subprocess.Popen('cp '+tmp+'/genome_all.zz.dedup.snv '+tmp+'/regular.snv',shell=True).wait()
            subprocess.Popen('cp '+tmp+'/genome_mskTC_all.zz.dedup.snv '+tmp+'/hyper_mskTC.snv',shell=True).wait()
            subprocess.Popen('cp '+tmp+'/genome_mskAG_all.zz.dedup.snv '+tmp+'/hyper_mskAG.snv',shell=True).wait()    """
        

        print 'identifying RESs...'
        if repeat !=False:
            annotate(tmp+'regular.snv',repeat,tmp+'regular.snv.anno')    
            #根据注释将snv分成3个部分
            seperate(tmp+'regular.snv.anno',tmp+'regular.snv.anno.alu',tmp+'regular.snv.anno.nalurp',tmp+'regular.snv.anno.nrp','Alu')
            get_snv_with_ad(tmp+'regular.snv.anno.alu',tmp+'regular.snv.anno.alu.ad2',2)
            annotate(tmp+'hyper_mskTC.snv',repeat,tmp+'hyper_mskTC.snv.anno')    
            seperate(tmp+'hyper_mskTC.snv.anno',tmp+'hyper_mskTC.snv.anno.alu',tmp+'hyper_mskTC.snv.anno.nalurp',tmp+'hyper_mskTC.snv.anno.nrp','Alu')
            annotate(tmp+'hyper_mskAG.snv',repeat,tmp+'hyper_mskAG.snv.anno')    
            seperate(tmp+'hyper_mskAG.snv.anno',tmp+'hyper_mskAG.snv.anno.alu',tmp+'hyper_mskAG.snv.anno.nalurp',tmp+'hyper_mskAG.snv.anno.nrp','Alu')

            #snv cluster定义为了在一个barcode中鉴定res，并且每一个结果的res按read来源作为单独的一行
            snv_cluster(tmp+'regular.snv.anno.alu',tmp+'regular_alu.res.ad1', cluster_distance, cluster_size_alu_ad1)
            snv_cluster(tmp+'regular.snv.anno.alu.ad2',tmp+'regular_alu.res.ad2', cluster_distance, cluster_size_alu_ad2)
            ##############修改至此处！！！！！！

            #bed or按照res位置、类型、read信息为唯一标识，合并ad1与ad2（有read信息限制，其实也无法进行加和
            bed_or(tmp+'regular_alu.res.ad1',tmp+'regular_alu.res.ad2',tmp+'regular_alu.res')
            #nalurp
            snv_cluster(tmp+'regular.snv.anno.nalurp',tmp+'regular_nalurp.res', cluster_distance, cluster_size_nalurp)
            #nrp 
            snv_cluster(tmp+'regular.snv.anno.nrp',tmp+'regular_nrp.res', cluster_distance, cluster_size_nrp)
            combine_res(tmp+'regular_alu.res',tmp+'regular_nalurp.res',tmp+'regular_nrp.res',tmp+'regular_split.res')
            cluster_size_regular_max=max([cluster_size_alu_ad1,cluster_size_alu_ad2,cluster_size_nalurp,cluster_size_nrp])
            print cluster_size_regular_max
            #三类snv共用cluster_distance+cluster_size_regular_max筛选条件
            combine_res(tmp+'regular.snv.anno.alu',tmp+'regular.snv.anno.nalurp',tmp+'regular.snv.anno.nrp',tmp+'regular.snv.anno.rmsrp')
            snv_cluster(tmp+'regular.snv.anno.rmsrp',tmp+'regular_overall.res', cluster_distance, cluster_size_regular_max)
            res_or(tmp+'regular_split.res',tmp+'regular_overall.res',tmp+'regular.res')
            #将res按细胞合并，利于下游计算深度
            get_res_in_1cell(tmp+'regular.res_readInfo', tmp+'regular.res_2getDepth')
            
            #maskTC
            snv_cluster(tmp+'hyper_mskTC.snv.anno.alu',tmp+'hyper_mskTC_alu.res', cluster_distance, cluster_size_alu_hp)
            snv_cluster(tmp+'hyper_mskTC.snv.anno.nalurp',tmp+'hyper_mskTC_nalurp.res', cluster_distance, cluster_size_nalurp_hp)
            snv_cluster(tmp+'hyper_mskTC.snv.anno.nrp',tmp+'hyper_mskTC_nrp.res', cluster_distance, cluster_size_nrp_hp)
                
            combine_res(tmp+'hyper_mskTC_alu.res',tmp+'hyper_mskTC_nalurp.res',tmp+'hyper_mskTC_nrp.res',tmp+'hyper_mskTC_split.res')            
            cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
            combine_res(tmp+'hyper_mskTC.snv.anno.alu',tmp+'hyper_mskTC.snv.anno.nalurp',tmp+'hyper_mskTC.snv.anno.nrp',tmp+'hyper_mskTC.snv.anno.rmsrp')
                            
            snv_cluster(tmp+'hyper_mskTC.snv.anno.rmsrp',tmp+'hyper_mskTC_overall.res', cluster_distance, cluster_size_hyper_max)          
            res_or(tmp+'hyper_mskTC_split.res',tmp+'hyper_mskTC_overall.res',tmp+'hyper_mskTC.res')

            #maskAG
            snv_cluster(tmp+'hyper_mskAG.snv.anno.alu',tmp+'hyper_mskAG_alu.res', cluster_distance, cluster_size_alu_hp)
            snv_cluster(tmp+'hyper_mskAG.snv.anno.nalurp',tmp+'hyper_mskAG_nalurp.res', cluster_distance, cluster_size_nalurp_hp)
            snv_cluster(tmp+'hyper_mskAG.snv.anno.nrp',tmp+'hyper_mskAG_nrp.res', cluster_distance, cluster_size_nrp_hp)

            combine_res(tmp+'hyper_mskAG_alu.res',tmp+'hyper_mskAG_nalurp.res',tmp+'hyper_mskAG_nrp.res',tmp+'hyper_mskAG_split.res')            
            cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
            combine_res(tmp+'hyper_mskAG.snv.anno.alu',tmp+'hyper_mskAG.snv.anno.nalurp',tmp+'hyper_mskAG.snv.anno.nrp',tmp+'hyper_mskAG.snv.anno.rmsrp')
            snv_cluster(tmp+'hyper_mskAG.snv.anno.rmsrp',tmp+'hyper_mskAG_overall.res', cluster_distance, cluster_size_hyper_max)
            res_or(tmp+'hyper_mskAG_split.res',tmp+'hyper_mskAG_overall.res',tmp+'hyper_mskAG.res') 

            snv_or(tmp+'hyper_mskTC.res',tmp+'hyper_mskAG.res',tmp+'hyper.res')

            snv_or(tmp+'/regular.res',tmp+'/hyper.res',tmp+'/all.res')


     
        else:
            snv_cluster(tmp+'regular.snv',tmp+'regular.res_tmp',cluster_distance,cluster_size_rg) 
            o2b(tmp+'regular.res_tmp',tmp+'regular.res') 

            snv_cluster(tmp+'hyper_mskTC.snv',tmp+'hyper_mskTC.res',cluster_distance,cluster_size_hp)
            snv_cluster(tmp+'hyper_mskAG.snv',tmp+'hyper_mskAG.res',cluster_distance,cluster_size_hp)
            snv_or(tmp+'hyper_mskTC.res',tmp+'hyper_mskAG.res',tmp+'hyper.res') 

        """ try:
            subprocess.Popen('rm -rf '+tmp+'/*.anno.*',shell=True).wait()
        except Exception, e:
            pass """  

        '''
        if repeat !=False:
            get_res(tmp+'regular_alu.res',tmp+'regular_nalurp.res',tmp+'regular_nrp.res',  tmp+'hyper.res',  tmp+'SPRINT_identified')
            bed_sort(tmp+'SPRINT_identified_A_to_I_regular.res',tmp+'SPRINT_identified_A_to_I_regular.res_sort')
            bed_sort(tmp+'SPRINT_identified_A_to_I_hyper.res',tmp+'SPRINT_identified_A_to_I_hyper.res_sort')
            bed_sort(tmp+'SPRINT_identified_C_to_U.res',tmp+'SPRINT_identified_C_to_U.res_sort')
        
            o2b(tmp+'SPRINT_identified_A_to_I_regular.res_sort',output+'SPRINT_identified_A_to_I_regular.res')
            o2b(tmp+'SPRINT_identified_A_to_I_hyper.res_sort',output+'SPRINT_identified_A_to_I_hyper.res')
            o2b(tmp+'SPRINT_identified_C_to_U.res_sort',output+'SPRINT_identified_C_to_U.res')
        
            snv_or(tmp+'SPRINT_identified_A_to_I_hyper.res',tmp+'SPRINT_identified_A_to_I_regular.res',output+'SPRINT_identified_A_to_I_all.res') 
        '''

        get_depth( tmp+'/all_combined.zz.sorted' , tmp+'/regular.res', tmp+'/regular.res.depth')
        #subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tbarcode\tAD:DP\n" | cat - '+tmp +'regular.res.depth   > '+output+'/SPRINT_identified_regular.res_onecell',shell=True).wait()
        get_depth( tmp+'/all_combined.zz.sorted' , tmp+'/hyper.res', tmp+'/hyper.res.depth')
        #subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tbarcode\tAD:DP\n" | cat - '+tmp +'/hyper.res.depth   > '+output+'/SPRINT_identified_hyper.res_onecell',shell=True).wait()
        #subprocess.Popen('cp '+tmp+'/PARAMETER.txt '+output+'/PARAMETER.txt',shell=True).wait()


        snv_or(tmp+'/regular.res',tmp+'/hyper.res',tmp+'/all.res') 
        get_depth( tmp+'/all_combined.zz.sorted' , tmp+'/all.res', tmp+'/all.res.depth')
        #subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tbarcode\tAD:DP\n" | cat - '+tmp +'/all.res.depth   > '+output+'/SPRINT_identified_all.res_onecell',shell=True).wait()
        print 'finished !'
        
        
        sys.exit(0) 
    try:
        pass
    except Exception,e:
        print ''
        print 'ERROR!'
        print ''
        print e
        print ''
        help_doc() 
    
main()
    
    
    
    
    
    
 
#if __name__=='__main__':   
#    main()
