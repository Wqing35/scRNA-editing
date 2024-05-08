#!/usr/bin/python2
# coding=utf-8

import numpy as np


###对regular_combined.zz.sorted.modified根据geneWzRes的barcode进行筛选
def main():

    def get_celltype_zz(zz_in_dir,geneWzRes_ad_in_dir,zz_out_dir):
        fi_geneWzRes=open(geneWzRes_ad_in_dir)

        barcodes=[]
        for line in fi_geneWzRes:
            seq=line.rstrip().split('\t')
            barcode=seq[6]
            barcodes.append(barcode)
        fi_geneWzRes.close()
        
        uniq_barcodes=np.unique(np.array(barcodes))

        fi_zz=open(zz_in_dir)
        fo=open(zz_out_dir,'w')
        for line in fi_zz:
            seq=line.rstrip().split('\t')
            read_id=seq[8]
            barcode=read_id.split('_')[2]
            if barcode in uniq_barcodes:
                fo.write(line)
        fi_zz.close()
        fo.close()       


    #samples='20'    
    samples=('17','18','57')
    celltypes=('L2_3_ExN','L4','L5_6','L5_6_CC','Olig','IN_VIP')
    #celltypes=('L4','L5_6','L5_6_CC','Olig','IN_VIP')
    for one_sample in samples:
        for one_celltype in celltypes:
            print(one_sample+' and '+one_celltype+' running!')
            zz_in_dir='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AEI/'+one_sample+'/regular_combined.zz'
            geneWzRes_ad_in_dir='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_'+one_sample+'/'+one_celltype+'_geneWithRegular_res_over1.txt'
            zz_out_dir='/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AEI/'+one_sample+'/'+'celltype_AEI/'+one_celltype+'/regular_combined.zz'
            
            get_celltype_zz(zz_in_dir,geneWzRes_ad_in_dir,zz_out_dir)
    print('Done')
main()

    
            