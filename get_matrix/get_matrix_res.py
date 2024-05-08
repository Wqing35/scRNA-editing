#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

def main(): 

    def get_matrix(bed_in_dir,barcode_in_dir,ad_out_dir,dp_out_dir):
        fi=open(bed_in_dir)
        #创建一个h5文件，文件指针是fo
        line=fi.readline()
        rowname=[]
        colname=[]
        while line !='':
            seq=line.rstrip().split('\t')
            chr_end=seq[0]+'_'+seq[2]
            #print chr_end
            barcode_id=seq[6]
            rowname.append(chr_end)
            colname.append(barcode_id)
            line=fi.readline()
        fi.close()

        #对barcode_id及RES去重
        srt_uniq_rowname=np.sort(np.unique(np.array(rowname)))
        srt_uniq_colname=np.sort(np.unique(np.array(colname)))

        ###barcode数过大，先导进过滤的barcode做筛选
        #####读取barcode文件
        barcode_file=pd.read_csv(barcode_in_dir,sep='\n',header=None).replace('-1','',regex=True)
        filtered_srt_uniq_colname=[val for val in np.array(barcode_file[0]) if val in srt_uniq_colname]
        #创建全0数组    
        matrix_ad=np.zeros((len(srt_uniq_rowname),len(filtered_srt_uniq_colname)),dtype='int8')
        matrix_dp=np.zeros((len(srt_uniq_rowname),len(filtered_srt_uniq_colname)),dtype='int8')
        df_ad=pd.DataFrame(matrix_ad, columns=filtered_srt_uniq_colname, index=srt_uniq_rowname) 
        df_dp=pd.DataFrame(matrix_dp, columns=filtered_srt_uniq_colname, index=srt_uniq_rowname) 
        
        #赋值
        ###改写成字典提取元素(缓)
        fi=open(bed_in_dir)    
        line=fi.readline()
        while line !='':
            seq=line.rstrip().split('\t')
            chr_end=seq[0]+'_'+seq[2]
            barcode_id=seq[6] 
            if(barcode_id in filtered_srt_uniq_colname):
                ad_dp=seq[7]
                ad=int(ad_dp.split(":")[0])
                dp=int(ad_dp.split(":")[1])
                df_ad.at[chr_end,barcode_id]=int(df_ad.at[chr_end,barcode_id]+ad)
                df_dp.at[chr_end,barcode_id]=int(df_dp.at[chr_end,barcode_id]+dp)
            line=fi.readline()
        fi.close()

        #转成稀疏矩阵
        df_ad_s=sparse.csr_matrix(df_ad) 
        df_dp_s=sparse.csr_matrix(df_dp) 

        #ad matrix
        with h5py.File(ad_out_dir,'w') as f_ad:
            df_ad_s.sort_indices()  # Compatibility with R
            ad_expr = f_ad.create_group("exprs")
            ad_expr.create_dataset("data", data=df_ad_s.data)
            ad_expr.create_dataset("indices", data=df_ad_s.indices)
            ad_expr.create_dataset("indptr", data=df_ad_s.indptr)
            ad_expr.create_dataset("shape", data=df_ad_s.shape)
            
            ad_barcode = f_ad.create_group('barcode')
            ad_barcode.create_dataset('barcode',data=filtered_srt_uniq_colname)
            ad_res = f_ad.create_group('site')
            ad_res.create_dataset('site',data=df_ad.index.values.astype(str))

        #dp matrix
        with h5py.File(dp_out_dir,'w') as f_dp:
            df_dp_s.sort_indices()  # Compatibility with R
            dp_expr = f_dp.create_group("exprs")
            dp_expr.create_dataset("data", data=df_dp_s.data)
            dp_expr.create_dataset("indices", data=df_dp_s.indices,)
            dp_expr.create_dataset("indptr", data=df_dp_s.indptr)
            dp_expr.create_dataset("shape", data=df_dp_s.shape)
            
            dp_barcode = f_dp.create_group('barcode')
            dp_barcode.create_dataset('barcode',data=filtered_srt_uniq_colname)
            dp_res = f_dp.create_group('site')
            dp_res.create_dataset('site',data=df_ad.index.values.astype(str))

    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/SPRINT_identified_regular.res','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/17_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/17/matrix_ad_over0.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/17/matrix_dp_over0.h5')
    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/SPRINT_identified_regular.res','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/57_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/57/matrix_ad_over0.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/57/matrix_dp_over0.h5')
    """ get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/20_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/all/20/matrix_ad.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/all/20/matrix_dp.h5')
    get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/17_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/17/matrix_ad.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/17/matrix_dp.h5')

    get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_32/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/32_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/all/32/matrix_ad.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/all/32/matrix_dp.h5')
    get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/18_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/18/matrix_ad.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/18/matrix_dp.h5')
    get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/56_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/all/56/matrix_ad.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/all/56/matrix_dp.h5')
    get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/57_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/57/matrix_ad.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/57/matrix_dp.h5') """


    get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_32/regular_merged.res_over1','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/32_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/32/matrix_ad_over1.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/32/matrix_dp_over1.h5')

    """ get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20/regular_merged.res_over0','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/20_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/20/matrix_ad_over0.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/20/matrix_dp_over0.h5')
    get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_20/regular_merged.res_over1','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/20_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/20/matrix_ad_over1.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/20/matrix_dp_over1.h5')

    get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56/regular_merged.res_over0','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/56_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/56/matrix_ad_over0.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/56/matrix_dp_over0.h5')
    get_matrix('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/ctr_56/regular_merged.res_over1','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/56_barcodes.csv','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/56/matrix_ad_over1.h5','/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AdDp_mat/regular/56/matrix_dp_over1.h5') """


    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/regular_merged.res_over1','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/17_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/17/matrix_ad_over1.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/17/matrix_dp_over1.h5')
    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/regular_merged.res_over1','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/57_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/57/matrix_ad_over1.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/57/matrix_dp_over1.h5')

    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/regular_merged.res_over2','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/17_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/17/matrix_ad_over2.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/17/matrix_dp_over2.h5')
    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/regular_merged.res_over2','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/57_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/57/matrix_ad_over2.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/regular/57/matrix_dp_over2.h5')


    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/17_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/17/matrix_ad_over0.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/17/matrix_dp_over0.h5')
    #get_matrix('/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/SPRINT_identified_all.res','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/57_barcodes.csv','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/57/matrix_ad_over0.h5','/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AdDp_mat/all/57/matrix_dp_over0.h5')


main()

