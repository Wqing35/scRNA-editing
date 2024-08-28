#!/usr/bin/python2
# coding=utf-8

import pandas as pd
import re

def main():

    #所有read比对的位置都在转录本的3’端，因此所有start位置都是read的坐标
    #无需更改坐标
    """ def modify_sam_loc(sam_in_dir=0,gtf_in_dir=0,sam_out_dir=0):
        annotation=pd.read_csv(gtf_in_dir,seq='\t')
        fi=open(sam_in_dir)
        for line in fi:
            seq=line.rstrip().split('\t')
            mapping_loc=int(seq[3])
            trans_id=seq[2].split('|')[0]
            #read_wz_iso_info[read_wz_iso_info.iloc[:, 1] == transcript].iloc[:, 0]
            genome_loc_start_for_trans=annotation[annotation.ilocp[:,11] == trans_id].iloc[:,3]
            new_ts_start=int(genome_loc_start_for_trans)+mapping_loc """
            


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

    refgenome='/disk1/wenqing/tmp_data/SUPPA2/ref/gencode.v37.transcripts.fa'
    for sample in ['GW12','GW16_1_3', 'GW16_1_9', 'GW19_1_1', 'GW19_1_2', 'GW19_1_3', 'GW23_1_1', 'GW23_1_2', 'GW23_1_3', 'GW26_1_1']:
        ori_sam_file='/disk1/wenqing/tmp_data/PFC_s2/salmon_output_untrimmed_ver/'+sample+'/pure_mappings.sam'
        out_zz_dir='/disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/EI_in_reads/'+sample+'/ts_mapping.zz'
        sam2zz(ori_sam_file,refgenome,out_zz_dir)

main()