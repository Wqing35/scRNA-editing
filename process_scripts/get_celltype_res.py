#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import pandas as pd
import numpy as np

####提取每种细胞类型的res

def main():

    def get_celltype_res(barcode_in_dir,res_in_dir,res_out_dir):
        final_barcode=pd.read_csv(barcode_in_dir,sep='\n',header=None).replace('-1','',regex=True)

        fi=open(res_in_dir)
        fo=open(res_out_dir,'w')
        for line in fi:
            seq=line.rstrip().split('\t')
            barcode=seq[6]
            if barcode in np.array(final_barcode):
                fo.write(line.rstrip()+'\n')
        fi.close()
        fo.close()

    input_path='/disk1/wenqing/tmp_data/ASD/'
    samples=["asd_male_pfc/asd_17/","asd_male_pfc/asd_18/","asd_male_pfc/asd_57/","ctr_male_pfc/ctr_20/","ctr_male_pfc/ctr_32/","ctr_male_pfc/ctr_56/"]
    #samples=['asd_male_pfc/asd_17/']
    #celltype=["L2_3_ExN","IN_VIP"]
    """ celltype=["Astrocyte","L5_6","L5_6_CC","Olig"]
    for one in celltype:
        for tmp in samples:
            print(tmp+' And '+one+' running:')
            input=input_path+tmp+'regular_merged.res_over1'
            output=input_path+tmp+'regular_merged.res_over1_'+one
            barcode_in_dir=input_path+tmp+one+'_barcodes.txt'
            get_celltype_res(barcode_in_dir,input,output) """


    #celltype=["L4","IN_VIP","L2_3"]
    celltype=["L2_3_ExN"]
    for one in celltype:
        for tmp in samples:
            print(tmp+' And '+one+' running:')
            input=input_path+tmp+'SPRINT_identified_all.res'
            output=input_path+tmp+one+'.res'
            barcode_in_dir=input_path+tmp+one+'_barcodes.txt'
            get_celltype_res(barcode_in_dir,input,output)
    print('done')

main()