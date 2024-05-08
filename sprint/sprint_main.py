#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import subprocess,os,sys
import re

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
        print "      sprint main   [options]   reference_genome(.fa)   output_path   bwa_path   samtools_path"
        print ""
        print "      options:"
        print "         -1       read1(.fq)       # Required !!!"
        print "         -2       read2(.fq)       # Optional"
        print "         -rp      repeat_file      # Optional, you can download it from http://sprint.software/SPRINT/dbrep/"
        print "         -ss      INT              # when input is strand-specific sequencing data, please clarify the direction of read1. [0 for antisense; 1 for sense] (default is 0)"
        #print "         -b       INT             # the format of read file [0: fq, 1: bam] (default is 0)"
        print "         -c       INT              # Remove the fist INT bp of each read (default is 0)"
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
        print "       sprint main -rp hg19_repeat.txt -c 6 -p 6 -1 read1.fq -2 read2.fq hg19.fa output ./bwa-0.7.12/bwa ./samtools-1.2/samtools"
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
    cluster_distance=200
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
    
    paired_end=False
    repeat=False
    #options：参数所在位置序号的列表
    options=[]
    read2=''
    read1=''
    #print sys.argv
    #sprint main -rp hg19_repeat.txt -c 6 -p 6 -1 read1.fq -2 read2.fq hg19.fa output ./bwa-0.7.12/bwa ./samtools-1.2/samtools
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

    #bwa、samtools所在路径
    bwa=sys.argv[all_argv[2]]
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
    
    def cut(fq_in_dir=0,fq_out_dir=0,cutnum=0,name='read1'):
        fi=open(fq_in_dir)
        fo=open(fq_out_dir,'w')
        cutnum=int(cutnum)
        line1=fi.readline()
        line2=fi.readline()
        line3=fi.readline()
        line4=fi.readline()
        idd=1
        while line1 !='':
            CELL_TAG=''
            if "XC:Z:" in line1:
                seq=line1.split('_')
                for one in seq:
                    if one[:5]=='XC:Z:':
                        CELL_TAG=one        
            if CELL_TAG !='':
                fo.write('@id_'+str(idd)+'_'+CELL_TAG+'_'+name+'\n')
            else:
                fo.write('@id_'+str(idd)+'_'+name+'\n')
            fo.write(line2[cutnum:])
            fo.write(line3)
            fo.write(line4[cutnum:])
            line1=fi.readline()
            line2=fi.readline()
            line3=fi.readline()
            line4=fi.readline()
            idd=idd+1
        fi.close()
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
    
    def fq2sam(TAG,paired_end,read1,read2,tmp,refgenome,bwa,samtools,mapcpu, read_format):
        if paired_end==True:
            mapcpu=max([int(int(mapcpu)/2.0),1])
        ori_tmp=tmp
        tmp=tmp+'/'+TAG+'/'
        if os.path.exists(tmp)==False:
            os.mkdir(tmp)
        
        #bwa分别构建样本基因组在参考基因组上的索引
        step1_1=subprocess.Popen(bwa+' aln -t '+str(mapcpu)+' '+refgenome+' '+read1+' > '+tmp+'read1.sai',shell=True)
        if paired_end==True:
            step1_2=subprocess.Popen(bwa+' aln  -t '+str(mapcpu)+' '+refgenome+' '+read2+' > '+tmp+'read2.sai',shell=True)

        step1_1.wait()
        if paired_end==True:
            step1_2.wait()
        #bwa samse生成sam文件
        step1_3=subprocess.Popen(bwa+' samse -n4 '+refgenome+' '+tmp+'read1.sai '+read1+' > '+tmp+'name_read1.sam',shell=True)
        if paired_end==True:
            step1_4=subprocess.Popen(bwa+' samse -n4 '+refgenome+' '+tmp+'read2.sai '+read2+' > '+tmp+'name_read2.sam',shell=True)
        step1_3.wait()
        if paired_end==True:
            step1_4.wait()
        #生成sam文件完成后，删除对应的索引文件、碱基剪切后的fq文件
        if os.path.exists(tmp+'name_read1.sam'):
            if os.path.exists(tmp+'read1.sai'):
                    os.remove(tmp+'read1.sai')
            if os.path.exists(ori_tmp+'cut_read1.fastq'):
                    os.remove(ori_tmp+'cut_read1.fastq')
        if os.path.exists(tmp+'name_read2.sam'):
            if os.path.exists(tmp+'read2.sai'):
                    os.remove(tmp+'read2.sai')
            if os.path.exists(ori_tmp+'cut_read2.fastq'):
                    os.remove(ori_tmp+'cut_read2.fastq')

        #sam转bam
        step1_7=subprocess.Popen(samtools+' view -bS '+tmp+'name_read1.sam >'+tmp+'name_read1.bam',shell=True)
        if paired_end==True:
            step1_8=subprocess.Popen(samtools+' view -bS '+tmp+'name_read2.sam >'+tmp+'name_read2.bam',shell=True)
        step1_7.wait()
        if paired_end==True:
            step1_8.wait()
        #排序
        if paired_end==True:
            step1_9=subprocess.Popen(samtools+' sort '+tmp+'name_read1.bam '+tmp+'name_read1_sorted',shell=True)
            step1_10=subprocess.Popen(samtools+' sort '+tmp+'name_read2.bam '+tmp+'name_read2_sorted',shell=True)
            step1_9.wait()
            step1_10.wait()
            #read1与read2合并成all.bam(结果文件名称固定)
            step1_11=subprocess.Popen(samtools+' merge -f '+tmp+'all.bam '+tmp+'name_read1_sorted.bam '+tmp+'name_read2_sorted.bam',shell=True)
            step1_11.wait()
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
            step1_9=subprocess.Popen(samtools+' sort '+tmp+'name_read1.bam '+tmp+'all',shell=True)
            step1_9.wait()
            if os.path.exists(tmp+'all.bam'):
                    if os.path.exists(tmp+'name_read1.sam'):
                            os.remove(tmp+'name_read1.sam')
                    if os.path.exists(tmp+'name_read1.bam'):
                            os.remove(tmp+'name_read1.bam')
        step2_2=subprocess.Popen(samtools+' view -h -o '+tmp+'all.sam '+tmp+'all.bam',shell=True)
        step2_2.wait()
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
            #若是fq文件转换描述行显示该read为read1，则将其转换为反义链；用作下文针对mask_from的碱基计数
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
        
            #mask操作中，将描述行改写成：‘原有描述行’_|_A_to_G_|_‘对应read的flag值’_|_read2
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
            #change the sam flag 'seq[1]' when you didn't use "bwa -aln" as mapper
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
                


                
                if len(baseqlst)>0 and sum(baseqlst)/float(len(baseqlst)) >= limitbasequa: #and missnum <= mismatch_num(len(seq[7])):
                    for snv in truesnv:
                        try:
                            allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]][0]=allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]][0]+1
                            if (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=='1' and snv[2][-1] == '1' ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=='0') and snv[2][-1] == '2' ):
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]][1]=allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]][1]+1
                            elif (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=='1' and snv[2][-1] == '2' ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=='0') and snv[2][-1] == '1' ):
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]][2]=allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]][2]+1
                        except Exception, e:
                            if (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=='1' and snv[2][-1] == '1' ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=='0') and snv[2][-1] == '2' ):
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]]=[1,1,0]
                            elif (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=='1' and snv[2][-1] == '2' ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=='0') and snv[2][-1] == '1' ):
                                allsnv[seq[0]+'\t'+snv[0].split(':')[0]+'\t'+snv[0].split(':')[1]]=[1,0,1]

        snv_bed=[]
        for snv in allsnv:
            seq=snv.split('\t')
            if allsnv[snv][0]>=limitad:
                    if allsnv[snv][1] > allsnv[snv][2]: 
                        snv_bed.append([seq[0],int(seq[2]),seq[1],'+',allsnv[snv][0]])
                    elif allsnv[snv][2] > allsnv[snv][1]:
                        snv_bed.append([seq[0],int(seq[2]),seq[1],'-',allsnv[snv][0]])
                    else:
                        snv_bed.append([seq[0],int(seq[2]),seq[1],'.',allsnv[snv][0]])

        snv_bed.sort()
        for one in snv_bed:
            fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+str(one[4])+'\t'+one[3]+'\n')
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
                whole[seq[0]+':'+seq[2]+':'+seq[3]+':'+seq[5]] += int(seq[4])
            except Exception,e:
                whole[seq[0]+':'+seq[2]+':'+seq[3]+':'+seq[5]] = int(seq[4])

        tmp=[]
        for one in whole:
            seq=one.split(':')
            dep=str(whole[one])
            tmp.append([seq[0],int(seq[1]),seq[2],dep,seq[3]])
        tmp.sort()
        for one in tmp:
            fo.write(one[0]+'\t'+str(one[1]-1)+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\t'+one[4]+'\n')
    
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
        
    
    #try:
    if 1==1: 

        #预处理阶段
        print 'preprocessing...'
        if read_format !=0:
            subprocess.Popen(samtools+' view -h -o '+tmp+'read1.sam '+read1,shell=True).wait()
            sprint.sam2fq(tmp+'read1.sam', tmp+'read1.fq')
            read1=tmp+'read1.fq'
            sprint.cut(read1,tmp+'cut_read1.fastq',cutbp,'read1')
            if paired_end==True:
                subprocess.Popen(samtools+' view -h -o '+tmp+'read2.sam '+read2,shell=True).wait()
                sprint.sam2fq(tmp+'read2.sam', tmp+'read2.fq')
                read2=tmp+'read2.fq'
                sprint.cut(read2,tmp+'cut_read2.fastq',cutbp,'read2')             
        else:
            #一般是从此开始，即read的碱基剪切开始
            if strand_specify==0:
                cut(read1,tmp+'cut_read1.fastq',cutbp,'read1')
                if paired_end==True:
                    cut(read2,tmp+'cut_read2.fastq',cutbp,'read2')
            else:
                sprint.cut(read1,tmp+'cut_read1.fastq',cutbp,'read2')
                if paired_end==True:
                    sprint.cut(read2,tmp+'cut_read2.fastq',cutbp,'read1')

        #决定后续碱基质量的阈值——测试数据为89
        get_baseq_cutoff(read1,tmp+'baseq.cutoff')
                    
        print 'mapping...'
    
        #step1_1：将raw fastq文件经处理后转换成bam文件
        #处理好的read1与read2合并的bam文件防于output/tmp/genome目录下
        #构建基因组索引时，事例命令行分配3个进程数
        TAG='genome'
        fq2sam(TAG,paired_end,tmp+'cut_read1.fastq',tmp+'cut_read2.fastq',tmp,refgenome,bwa,samtools,mapcpu,read_format) 
     
        #step1_2
        #将未匹配到基因组上的reads从bam文件中提取出来、并转换为sam格式(包含了read1和read2总共的unmapped reads)
        subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
        #将上述sam文件转换为fq文件
        umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/'+TAG+'_unmapped.fq')
     
        #若存在转录参考基因组，则将unmapped reads和转录组模版再比对一次（有什么必要？？？？？）
        #工具问题的选用，回帖的意义
        if os.path.exists(refgenome+'.trans.fa'):
            print "yes!!!!!"
            TAG='transcript'
            fq2sam(TAG,False,tmp+'/genome_unmapped.fq',read2,tmp,refgenome+'.trans.fa',bwa,samtools,mapcpu,read_format)
    
            subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
            umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/regular_unmapped.fq')    
            maskfq(tmp+'/regular_unmapped.fq','A','G')
        #若当前目录不存在转录参考基因组样本的情况下，则将tmp/genome/下的unmapped.sam文件转换为tmp/路径下的fq文件
        #且将unmapped reads上的A换成G
        else:
            umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/regular_unmapped.fq')
            maskfq(tmp+'/regular_unmapped.fq','A','G')
    
        #step2_1:
        #将step1中未匹配到基因组上的、经A-G转换的read加上read2一起，再在经A-G转换的参考基因组上比对一次         
        TAG='genome_mskAG'
        fq2sam(TAG,False,tmp+'/regular_unmapped_A_to_G.fq',read2,tmp,refgenome+'.mskAG.fa',bwa,samtools,mapcpu,read_format)
    
        subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
        umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/'+TAG+'_unmapped.fq')    
    
        #step2_2:
        # 将step1中未匹配到基因组上的、经A-G转换的read加上read2一起，再在经T-C转换的参考基因组上比对一次
        TAG='genome_mskTC'
        fq2sam(TAG,False,tmp+'/regular_unmapped_A_to_G.fq',read2,tmp,refgenome+'.mskTC.fa',bwa,samtools,mapcpu,read_format)
    
        subprocess.Popen(samtools+' view -f4 '+tmp+'/'+TAG+'/all.bam > '+tmp+'/'+TAG+'_unmapped.sam',shell=True).wait()
        umsam2fq(tmp+'/'+TAG+'_unmapped.sam',tmp+'/'+TAG+'_unmapped.fq')    
        
        #在有转录参考基因组样本的情况下，经A-G转换的reads与经A-G、T-C转换的转录参考基因组再匹配一次
        if os.path.exists(refgenome+'.trans.fa'):
            TAG='transcript_mskAG'
            fq2sam(TAG,False,tmp+'/genome_mskAG_unmapped.fq',read2,tmp,refgenome+'.trans.fa.mskAG.fa',bwa,samtools,mapcpu,read_format)
    
            TAG='transcript_mskTC'
            fq2sam(TAG,False,tmp+'/genome_mskTC_unmapped.fq',read2,tmp,refgenome+'.trans.fa.mskTC.fa',bwa,samtools,mapcpu,read_format)  
    
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

            #将候选snv在转录本模板的定位转换成基因组上的位置（用处？？？？？？？？）
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
            subprocess.Popen('cp '+tmp+'/genome_mskAG_all.zz.dedup.snv '+tmp+'/hyper_mskAG.snv',shell=True).wait()
        
        


        print 'identifying RESs...'

        #repeat：repeat file
        if repeat !=False:

            #根据repeat file对上述流程的snv进行注释，那些区域是重复等
            annotate(tmp+'regular.snv',repeat,tmp+'regular.snv.anno')    
            #根据注释将snv分成3个部分
            seperate(tmp+'regular.snv.anno',tmp+'regular.snv.anno.alu',tmp+'regular.snv.anno.nalurp',tmp+'regular.snv.anno.nrp','Alu')
            get_snv_with_ad(tmp+'regular.snv.anno.alu',tmp+'regular.snv.anno.alu.ad2',2)
            #筛选Alu snv的2条支线
            snv_cluster(tmp+'regular.snv.anno.alu',tmp+'regular_alu.res.ad1', cluster_distance, cluster_size_alu_ad1)
            snv_cluster(tmp+'regular.snv.anno.alu.ad2',tmp+'regular_alu.res.ad2', cluster_distance, cluster_size_alu_ad2)
            bed_or(tmp+'regular_alu.res.ad1',tmp+'regular_alu.res.ad2',tmp+'regular_alu.res')
            #筛选non-Alu snv
            snv_cluster(tmp+'regular.snv.anno.nalurp',tmp+'regular_nalurp.res', cluster_distance, cluster_size_nalurp)
            #筛选nrp snv
            snv_cluster(tmp+'regular.snv.anno.nrp',tmp+'regular_nrp.res', cluster_distance, cluster_size_nrp)
            #合并三类snv(cluster_distance+cluster size三者分开筛选后合并)
            combine_res(tmp+'regular_alu.res',tmp+'regular_nalurp.res',tmp+'regular_nrp.res',tmp+'regular_split.res')
            cluster_size_regular_max=max([cluster_size_alu_ad1,cluster_size_alu_ad2,cluster_size_nalurp,cluster_size_nrp])
            print cluster_size_regular_max
            #三类snv共用cluster_distance+cluster_size_regular_max筛选条件
            #alu
            combine_res(tmp+'regular.snv.anno.alu',tmp+'regular.snv.anno.nalurp',tmp+'regular.snv.anno.nrp',tmp+'regular.snv.anno.rmsrp')
            snv_cluster(tmp+'regular.snv.anno.rmsrp',tmp+'regular_overall.res', cluster_distance, cluster_size_regular_max)
            res_or(tmp+'regular_split.res',tmp+'regular_overall.res',tmp+'regular.res')
    

            #对hyper snv采取与regular snv同样的操作
            annotate(tmp+'hyper_mskTC.snv',repeat,tmp+'hyper_mskTC.snv.anno')    
            seperate(tmp+'hyper_mskTC.snv.anno',tmp+'hyper_mskTC.snv.anno.alu',tmp+'hyper_mskTC.snv.anno.nalurp',tmp+'hyper_mskTC.snv.anno.nrp','Alu')
            snv_cluster(tmp+'hyper_mskTC.snv.anno.alu',tmp+'hyper_mskTC_alu.res', cluster_distance, cluster_size_alu_hp)
            snv_cluster(tmp+'hyper_mskTC.snv.anno.nalurp',tmp+'hyper_mskTC_nalurp.res', cluster_distance, cluster_size_nalurp_hp)
            snv_cluster(tmp+'hyper_mskTC.snv.anno.nrp',tmp+'hyper_mskTC_nrp.res', cluster_distance, cluster_size_nrp_hp)
            combine_res(tmp+'hyper_mskTC_alu.res',tmp+'hyper_mskTC_nalurp.res',tmp+'hyper_mskTC_nrp.res',tmp+'hyper_mskTC_split.res')            
            cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
            combine_res(tmp+'hyper_mskTC.snv.anno.alu',tmp+'hyper_mskTC.snv.anno.nalurp',tmp+'hyper_mskTC.snv.anno.nrp',tmp+'hyper_mskTC.snv.anno.rmsrp')
            snv_cluster(tmp+'hyper_mskTC.snv.anno.rmsrp',tmp+'hyper_mskTC_overall.res', cluster_distance, cluster_size_hyper_max)
            res_or(tmp+'hyper_mskTC_split.res',tmp+'hyper_mskTC_overall.res',tmp+'hyper_mskTC.res')

            #为何hyper_mskAG.snv.anno文件大小为0，后续仍能出结果？？？？？？
            annotate(tmp+'hyper_mskAG.snv',repeat,tmp+'hyper_mskAG.snv.anno')    
            seperate(tmp+'hyper_mskAG.snv.anno',tmp+'hyper_mskAG.snv.anno.alu',tmp+'hyper_mskAG.snv.anno.nalurp',tmp+'hyper_mskAG.snv.anno.nrp','Alu')
            snv_cluster(tmp+'hyper_mskAG.snv.anno.alu',tmp+'hyper_mskAG_alu.res', cluster_distance, cluster_size_alu_hp)
            snv_cluster(tmp+'hyper_mskAG.snv.anno.nalurp',tmp+'hyper_mskAG_nalurp.res', cluster_distance, cluster_size_nalurp_hp)
            snv_cluster(tmp+'hyper_mskAG.snv.anno.nrp',tmp+'hyper_mskAG_nrp.res', cluster_distance, cluster_size_nrp_hp)
            combine_res(tmp+'hyper_mskAG_alu.res',tmp+'hyper_mskAG_nalurp.res',tmp+'hyper_mskAG_nrp.res',tmp+'hyper_mskAG_split.res')            
            cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
            combine_res(tmp+'hyper_mskAG.snv.anno.alu',tmp+'hyper_mskAG.snv.anno.nalurp',tmp+'hyper_mskAG.snv.anno.nrp',tmp+'hyper_mskAG.snv.anno.rmsrp')
            snv_cluster(tmp+'hyper_mskAG.snv.anno.rmsrp',tmp+'hyper_mskAG_overall.res', cluster_distance, cluster_size_hyper_max)
            res_or(tmp+'hyper_mskAG_split.res',tmp+'hyper_mskAG_overall.res',tmp+'hyper_mskAG.res')

            snv_or(tmp+'hyper_mskTC.res',tmp+'hyper_mskAG.res',tmp+'hyper.res') 


            
    
     
        else:
            snv_cluster(tmp+'regular.snv',tmp+'regular.res_tmp',cluster_distance,cluster_size_rg) 
            o2b(tmp+'regular.res_tmp',tmp+'regular.res') 

            snv_cluster(tmp+'hyper_mskTC.snv',tmp+'hyper_mskTC.res',cluster_distance,cluster_size_hp)
            snv_cluster(tmp+'hyper_mskAG.snv',tmp+'hyper_mskAG.res',cluster_distance,cluster_size_hp)
            snv_or(tmp+'hyper_mskTC.res',tmp+'hyper_mskAG.res',tmp+'hyper.res')

        try:
            subprocess.Popen('rm -rf '+tmp+'/*.anno.*',shell=True).wait()
        except Exception, e:
            pass

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

        # subprocess.Popen('cp '+tmp+'/regular.res '+output+'/SPRINT_identified_regular.res',shell=True).wait()
        get_depth( tmp+'/all_combined.zz.sorted' , tmp+'/regular.res', tmp+'/regular.res.depth')
        subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tAD:DP" | cat - '+tmp +'/regular.res.depth   > '+output+'/SPRINT_identified_regular.res',shell=True).wait()
       # subprocess.Popen('cp '+tmp+'/hyper.res '+output+'/SPRINT_identified_hyper.res',shell=True).wait()
        get_depth( tmp+'/all_combined.zz.sorted' , tmp+'/hyper.res', tmp+'/hyper.res.depth')
        subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tAD:DP" | cat - '+tmp +'/hyper.res.depth   > '+output+'/SPRINT_identified_hyper.res',shell=True).wait()
        subprocess.Popen('cp '+tmp+'/PARAMETER.txt '+output+'/PARAMETER.txt',shell=True).wait()

        #subprocess.Popen('grep "AG" '+tmp+'/hyper.res | grep "+"  > '+tmp+'/hyper_AG+.res',shell=True).wait()
        #subprocess.Popen('grep "TC" '+tmp+'/hyper.res | grep "-"  > '+tmp+'/hyper_TC-.res',shell=True).wait()
        #snv_or(tmp+'/hyper_AG+.res',tmp+'/hyper_TC-.res',tmp+'/hyper_AG.res') 
        #snv_or(tmp+'/regular.res',tmp+'/hyper_AG.res',tmp+'/all.res') 
        snv_or(tmp+'/regular.res',tmp+'/hyper.res',tmp+'/all.res') 
        get_depth( tmp+'/all_combined.zz.sorted' , tmp+'/all.res', tmp+'/all.res.depth')
        #AD:DP：某位点RNA编辑率
        subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tAD:DP" | cat - '+tmp +'/all.res.depth   > '+output+'/SPRINT_identified_all.res',shell=True).wait()
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

