#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


""" awk '{print $1"_"$3}' resInAluWithADDP_over1 > tmp_site.txt
paste -d '\t' resInAluWithADDP_over1 tmp_site.txt > resInAluWithADDP_over1.withSite
bedtools intersect -wb -a /local/wenqing/data/geneElements/gene.bed -b resInAluWithADDP_over1.withSite | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$12"\t"$9}' > geneWithRegular_res.txt
sed -i 's/;//g' geneWithRegular_res.txt
sed -i 's/\"//g' geneWithRegular_res.txt

rm tmp_site.txt
rm resInAluWithADDP_over1.withSite """

def main():

    def get_gene_matrix(bed_in_dir,barcode_in_dir,geneWithRes_in_dir,geneMat_ad_out_dir,geneMat_dp_out_dir,geneMat_out_dir):
        print "creating null matrix"
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

        print "calculating value for created matrix"
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


        #与细胞表达矩阵取交集后，过滤行
        #df_ad.apply(lambda x: x.sum(), axis=0)      #df_ada按照行与列分别求和无小于0的，因此不做过滤

        print "filtering res and cell"
        #过滤含有不超过3个res的细胞

        #apply效率太低，弃用
        #resInCell_loc = df_ad.apply(lambda x: np.sum(x > 0) > 3, axis=0)
        #仅出现在一个细胞中的res（过滤了80%），所以过滤条件改为大于0即可
        #res_loc = sub_ad_mat.apply(lambda x: np.sum(x > 0) > 0, axis=1)
  
        # filtering on the basis of columns
        sub_ad_mat = df_ad.loc[:, df_ad.sum(axis=0) > 3]
        # filtering on the basis of rows
        sub_ad_mat = sub_ad_mat[sub_ad_mat.sum(axis=1) > 0]
        sub_dp_mat = df_dp.loc[sub_ad_mat.index.values,sub_ad_mat.columns.values]


        print "filtering res accoring to genes"
        #过滤不在基因上的res
        geneWithRes=pd.read_csv(geneWithRes_in_dir,sep='\t',header=None)
        ResInGene=np.unique(geneWithRes.iloc[:,4].values)
        overlap_res=set(ResInGene).intersection(set(sub_ad_mat.index))
        geneWithRes_filtered=geneWithRes[geneWithRes.iloc[:,4].isin(overlap_res)]

        #有部分res落在多个基因上，此处仅采取落在一个基因上的res形成gene mat
        final_res=[]
        res_dict={}
        for res in overlap_res:
            res_dict[res]=np.unique(np.array(geneWithRes_filtered.iloc[:,3][geneWithRes_filtered.iloc[:,4]==one_res]))
            if len(res_dict[res])==1:
                final_res.append(res)
        final_res=set(final_res)

        resOngenes_sub_ad_mat = sub_ad_mat.loc[final_res,:] #筛选得到res 1171518
        resOngenes_sub_dp_mat = sub_dp_mat.loc[final_res,:] 
        geneWithRes_filtered=geneWithRes_filtered[geneWithRes_filtered.iloc[:,4].isin(final_res)]


        resOngenes_sub_mat=resOngenes_sub_ad_mat
        #下面的语句将resOngenes_sub_ad_mat改变，故重新对其赋值一次
        resOngenes_sub_mat[resOngenes_sub_mat!=0]=1
        resOngenes_sub_ad_mat = sub_ad_mat.loc[final_res,:] 

        ##生成基因矩阵
        _,idx=np.unique(geneWithRes_filtered.iloc[:,3],return_index=True)
        gene_name=geneWithRes_filtered.iloc[np.sort(idx),3].values

        print "finding res in genes"
        geneHasRes_dict_ad={}
        geneHasRes_dict_dp={}
        geneHasRes_dict={}
        for gene in gene_name:
            res=np.unique(geneWithRes_filtered.iloc[:,4][geneWithRes_filtered.iloc[:,3]==gene].values)
            geneHasRes_dict_ad[gene]=resOngenes_sub_ad_mat.loc[res,:].apply(lambda x: x.sum(), axis=0)#axis=0按行求和
            geneHasRes_dict_dp[gene]=resOngenes_sub_dp_mat.loc[res,:].apply(lambda x: x.sum(), axis=0)#axis=0按行求和
            geneHasRes_dict[gene]=resOngenes_sub_mat.loc[res,:].apply(lambda x: x.sum(), axis=0)

        gene_mat_ad=pd.DataFrame(geneHasRes_dict_ad)
        gene_mat_ad=gene_mat_ad.T
        gene_mat_dp=pd.DataFrame(geneHasRes_dict_dp)
        gene_mat_dp=gene_mat_dp.T
        gene_mat=pd.DataFrame(geneHasRes_dict)
        gene_mat=gene_mat.T

        sparse_gene_mat_ad=sparse.csr_matrix(gene_mat_ad) 
        sparse_gene_mat_dp=sparse.csr_matrix(gene_mat_dp) 
        sparse_gene_mat=sparse.csr_matrix(gene_mat) 

        print "writing to output"
        with h5py.File(geneMat_out_dir,'w') as f:
            sparse_gene_mat.sort_indices()  # Compatibility with R
            expr = f.create_group("exprs")
            expr.create_dataset("data", data=sparse_gene_mat.data)
            expr.create_dataset("indices", data=sparse_gene_mat.indices)
            expr.create_dataset("indptr", data=sparse_gene_mat.indptr)
            expr.create_dataset("shape", data=sparse_gene_mat.shape)
            
            genes = f.create_group('gene')
            genes.create_dataset('gene',data=gene_name.astype(str))
            barcode = f.create_group('barcode')
            barcode.create_dataset('barcode',data=gene_mat.columns.values.astype(str))

        with h5py.File(geneMat_ad_out_dir,'w') as f:
            sparse_gene_mat_ad.sort_indices()  # Compatibility with R
            expr = f.create_group("exprs")
            expr.create_dataset("data", data=sparse_gene_mat_ad.data)
            expr.create_dataset("indices", data=sparse_gene_mat_ad.indices)
            expr.create_dataset("indptr", data=sparse_gene_mat_ad.indptr)
            expr.create_dataset("shape", data=sparse_gene_mat_ad.shape)
            
            genes = f.create_group('gene')
            genes.create_dataset('gene',data=gene_name.astype(str))
            barcode = f.create_group('barcode')
            barcode.create_dataset('barcode',data=gene_mat.columns.values.astype(str))

        with h5py.File(geneMat_dp_out_dir,'w') as f:
            sparse_gene_mat_dp.sort_indices()  # Compatibility with R
            expr = f.create_group("exprs")
            expr.create_dataset("data", data=sparse_gene_mat_dp.data)
            expr.create_dataset("indices", data=sparse_gene_mat_dp.indices)
            expr.create_dataset("indptr", data=sparse_gene_mat_dp.indptr)
            expr.create_dataset("shape", data=sparse_gene_mat_dp.shape)
            
            genes = f.create_group('gene')
            genes.create_dataset('gene',data=gene_name.astype(str))
            barcode = f.create_group('barcode')
            barcode.create_dataset('barcode',data=gene_mat.columns.values.astype(str))              

    print('Start') 
    #asd
    #resInAluWithADDP_over1:位于ALU区域的res，计算gene的eiditing index不再考虑其alu部分？，（ALU区域的res站总res的绝大部份，限制这个条件与否可能不影响后续的分析？？？）
    #get_gene_matrix("/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AEI/17/resInAluWithADDP_over1","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/17_barcodes.csv","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AEI/17/geneWithRegular_res.txt","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/regular/17/gene_mat_ad_over1.h5","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/regular/17/gene_mat_dp_over1.h5","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/regular/17/gene_mat_over1.h5")
    get_gene_matrix("/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/SPRINT_identified_all.res","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/17_barcodes.csv","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_17/geneWithRes_all.txt","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/all/17/gene_mat_ad.h5","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/all/17/gene_mat_dp.h5","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/all/17/gene_mat.h5")
    #get_gene_matrix("/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AEI/57/resInAluWithADDP_over1","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/57_barcodes.csv","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/AEI/57/geneWithRegular_res.txt","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/regular/57/gene_mat_ad_over1.h5","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/regular/57/gene_mat_dp_over1.h5","/disk1/wenqing/tmp_data/ASD/asd_male_pfc/gene_mat/regular/57/gene_mat_over1.h5")
    #ctr
    #get_gene_matrix("/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AEI/20/resInAluWithADDP_over1","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/20_barcodes.csv","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AEI/20/geneWithRegular_res.txt","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/20/gene_mat_ad_over1.h5","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/20/gene_mat_dp_over1.h5","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/20/gene_mat_over1.h5")
    #get_gene_matrix("/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AEI/23/resInAluWithADDP_over1","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/23_barcodes.csv","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AEI/23/geneWithRegular_res.txt","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/23/gene_mat_ad_over1.h5","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/23/gene_mat_dp_over1.h5","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/23/gene_mat_over1.h5")
    #get_gene_matrix("/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AEI/56/resInAluWithADDP_over1","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/56_barcodes.csv","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AEI/56/geneWithRegular_res.txt","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/56/gene_mat_ad_over1.h5","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/56/gene_mat_dp_over1.h5","/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_mat/regular/56/gene_mat_over1.h5")


    print('Done')

main()


        