#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

####计算每个样本的SNR
#"ExN","Interneuron","L5/6-CC","Olig"
import pandas as pd
import numpy as np

sample=['17','18','57']
for phase in sample:
    cell_types=["ExN","Interneuron","L5_6_CC","Olig"]
    for type in cell_types:
        barcode_in_dir='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/'+phase+'_'+type+'_barcodes.txt'
        barcode_file=pd.read_csv(barcode_in_dir,sep='\n',header=None)

        res_in_dir='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_'+phase+'/regular_merged.res_over1'
        res_out_dir='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_'+phase+'/regular_merged.res_over1_'+type

        fi=open(res_in_dir)
        fo=open(res_out_dir,'w')
        for line in fi:
            seq=line.rstrip().split('\t')
            barcode=seq[6]
            if barcode in np.array(barcode_file[0]):
                fo.write(line)
        fi.close()
        fo.close()
