#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import numpy as np
import pandas as pd

#对每个细胞类型抽取10个细胞，每个sample共获得5*10个细胞，以证明RNA编辑与IR的相关性
def main():

    def get_specific_cell(fq_in_dir=0,fq_r1_in_dir=0,barcode_in_dir=0,fq_out_dir=0,cutnum=0,cutnum_r1=0):
        #barcode读取
        final_barcode=barcode_in_dir
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
        while line1 !='' and line1_r1 !='':
            barcode=line2_r1[:cutnum_r1]
            if barcode == final_barcode:
                fo.write(line1)
                fo.write(line2[cutnum:])
                fo.write(line3)
                fo.write(line4[cutnum:])
                #read迭代
            line1=fi.readline()
            line2=fi.readline()
            line3=fi.readline()
            line4=fi.readline()
            #barcode迭代
            line1_r1=fi_r1.readline()
            line2_r1=fi_r1.readline()
            line3_r1=fi_r1.readline()
            line4_r1=fi_r1.readline()
        fi.close()
        fi_r1.close()
        fo.close() 

    cutnum=0
    cutnum_r1=16
    fq_in_dir_17='/disk1/wenqing/tmp_data/ASD/SRR9262917_2.fastq'
    fq_r1_in_dir_17='/disk1/wenqing/tmp_data/ASD/SRR9262917_1.fastq'
    barcode_in_dir_17='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/random_sample_25_barcodes.txt'
    barcode_in_17=pd.read_csv(barcode_in_dir_17,sep='\n',header=None).replace('-1','',regex=True)
    for tmp_barcode in np.array(barcode_in_17):
        fq_out_dir='/disk1/wenqing/tmp_data/ASD/17_'+tmp_barcode[0]+'_2.fastq'
        get_specific_cell(fq_in_dir_17,fq_r1_in_dir_17,tmp_barcode[0],fq_out_dir,cutnum,cutnum_r1)

    fq_in_dir_18='/disk1/wenqing/tmp_data/ASD/SRR9262918_2.fastq'
    fq_r1_in_dir_18='/disk1/wenqing/tmp_data/ASD/SRR9262918_1.fastq'
    barcode_in_dir_18='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18/random_sample_25_barcodes.txt'
    barcode_in_18=pd.read_csv(barcode_in_dir_18,sep='\n',header=None).replace('-1','',regex=True)
    for tmp_barcode in np.array(barcode_in_18):
        fq_out_dir='/disk1/wenqing/tmp_data/ASD/18_'+tmp_barcode[0]+'_2.fastq'
        get_specific_cell(fq_in_dir_18,fq_r1_in_dir_18,tmp_barcode[0],fq_out_dir,cutnum,cutnum_r1)

    fq_in_dir_57='/disk1/wenqing/tmp_data/ASD/SRR9262957_2.fastq'
    fq_r1_in_dir_57='/disk1/wenqing/tmp_data/ASD/SRR9262957_1.fastq'
    barcode_in_dir_57='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/random_sample_25_barcodes.txt'
    barcode_in_57=pd.read_csv(barcode_in_dir_57,sep='\n',header=None).replace('-1','',regex=True)
    for tmp_barcode in np.array(barcode_in_57):
        fq_out_dir='/disk1/wenqing/tmp_data/ASD/57_'+tmp_barcode[0]+'_2.fastq'
        get_specific_cell(fq_in_dir_57,fq_r1_in_dir_57,tmp_barcode[0],fq_out_dir,cutnum,cutnum_r1)

main()
    
