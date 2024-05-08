#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import subprocess,os
import re
import pandas as pd
import numpy as np
import h5py  
import concurrent.futures
from scipy.sparse import csr_matrix 
from scipy import sparse




def main():

    #def get_AEI(bed_in_dir,bed_anno_dir,bed_out_dir,res_in_dir,AEI_out_dir)
    def get_AEI(element,region,cell_name):

        ###################获取gene区域的reads
        print 'Get reads in gene:'

        #筛选位于gene区域的reads        
        # 定义你的命令和参数  
        command = ['bedtools', 'intersect', '-wo', '-a', 'regular_combined.zz.sorted.modified.'+cell_name, '-b', 'geneWzRes_ad.'+cell_name+'.txt']  
        
        # 执行命令并将输出重定向到一个文件中  
        with open('regular_combined.zz.sorted.modified_InGene.'+cell_name, 'w') as output_file:  
            # 使用Popen执行命令，并捕获标准输出  
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)  
            # 等待命令执行完成，并获取输出  
            stdout_data, stderr_data = process.communicate()  
            
            # 检查命令是否执行成功  
            if process.returncode != 0:  
                print("Error executing command: ", stderr_data)  
            else:  
                # 将标准输出写入文件  
                output_file.write(stdout_data)

        ###################step4:根据gene范围，剪切reads长度
        print 'Trim reads according to Gene:'

        fi=open(element+'_combined.zz.sorted.modified_InGene.'+cell_name)
        fo=open(element+'_combined.zz.sorted.modified_InGene_trimmed_reads_OnlyAG.'+cell_name,'w')

        for line in fi:
            seq=line.rstrip().split("\t")
            read=seq[5]
            reads_start=int(seq[1])
            reads_end=int(seq[2])
            gene_start=int(seq[8])
            gene_end=int(seq[9])

            #剪切reads的几种可能的场景
            #1.该条read包含该gene
            if gene_start > reads_start and gene_end < reads_end:
                new_start=gene_start
                new_end=gene_end
                trimmed_read=read[(gene_start-reads_start-1):(gene_end-reads_start)]
            #2.该条read被包含于该gene
            elif gene_start < reads_start and gene_end > reads_end:
                new_start=reads_start
                new_end=reads_end
                trimmed_read=read[0:(reads_end-reads_start)]
            #3.该read的下游区域未覆盖到gene
            elif gene_start == reads_start and gene_end > reads_end:
                new_start=reads_start
                new_end=reads_end
                trimmed_read=read[0:(reads_end-reads_start)] 
            #4.该read的上游区域未覆盖到gene
            elif gene_start < reads_start and gene_end == reads_end:
                new_start=reads_start
                new_end=reads_end
                trimmed_read=read[0:(reads_end-reads_start)] 
            #5.该read下游过覆盖gene
            elif gene_start == reads_start and gene_end < reads_end:
                new_start=reads_start
                new_end=gene_end
                trimmed_read=read[0:(gene_end-gene_start)]   
            #6.该read上游过覆盖gene
            elif gene_start > reads_start and gene_end == reads_end:
                new_start=gene_start
                new_end=gene_end
                trimmed_read=read[(gene_start-reads_start):(reads_end-reads_start)]                 
            #7.该reads上游过覆盖、下游未覆盖
            elif gene_start > reads_start and gene_end > reads_end:
                new_start=gene_start
                new_end=reads_end
                trimmed_read=read[(gene_start-reads_start):(reads_end-reads_start)]                 
            #8.该reads上游未覆盖、下游过覆盖
            elif gene_start < reads_start and gene_end < reads_end:
                new_start=reads_start
                new_end=gene_end
                trimmed_read=read[0:(gene_end-reads_start)]                 

            snv=seq[3]
            AG_num=len(re.findall("AG:",snv))
            TC_num=len(re.findall("TC:",snv))
            As_num=AG_num+TC_num
            AA_num=len(re.findall("A",trimmed_read))
            #total_num=As_num+AA_num
            #AG mismatch > TC mismatch只算AG
            total_num=AG_num+AA_num
            gene_mark=seq[7]+'_'+seq[8]+'_'+seq[9]
            fo.write(seq[0]+'\t'+str(new_start)+'\t'+str(new_end)+'\t'+seq[3]+'\t'+seq[4]+'\t'+trimmed_read+'\t'+seq[10]+'\t'+seq[11]+'\t'+str(total_num)+'\t'+gene_mark+'\n')

        fi.close()
        fo.close()

        ###########step5:获取AEI
        print 'Get AEI:'

        gene_counts_in_dir=element+'_combined.zz.sorted.modified_InGene_trimmed_reads_OnlyAG.'+cell_name

        fi=open(gene_counts_in_dir)
        #######待改！！！！！！！！
        gene_dict1={}
        for line in fi:
            seq=line.rstrip().split('\t')
            gene_id=seq[6]  
            res_As_num=int(seq[7])
            reads_As_num=int(seq[8])
            try:
                gene_dict1[gene_id][0]=res_As_num
                gene_dict1[gene_id][1]=gene_dict1[gene_id][1]+reads_As_num
            except Exception, e:
                gene_dict1[gene_id]=[0,0]
                gene_dict1[gene_id][0]=res_As_num
                gene_dict1[gene_id][1]=reads_As_num
        fi.close()

        EI=[]
        for key in gene_dict1.keys():
            if gene_dict1[key][0]==0 or gene_dict1[key][1]==0:
                EI.append(float(0))
            else:
                EI.append(float(gene_dict1[key][0])/int(gene_dict1[key][1]))
        
        data=pd.DataFrame()
        #data["gene_name"]=gene_dict1.keys()
        data["EI"]=EI
        data['gene']=gene_dict1.keys()

        return(data)


    def process_cell(one_cell, fi_total_res, fi_total_reads):  
        # 初始化结果文件  
        with open('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_ei/20/tmp_one_cell/geneWzRes_ad.'+one_cell+'.txt', 'w') as fo:  
            # 提取每个细胞的 res  
            sub_res_in_one_cell = fi_total_res[fi_total_res.iloc[:, 6] == one_cell]  
            
            # 提取每个细胞的 reads  
            sub_reads_in_one_cell = fi_total_reads[fi_total_reads.iloc[:, 6] == one_cell]  
            
            # 生成唯一的基因区域  
            uniq_gene_region = np.unique(  
                sub_res_in_one_cell.iloc[:, 0].astype(str) + '_' +  
                sub_res_in_one_cell.iloc[:, 1].astype(str) + '_' +  
                sub_res_in_one_cell.iloc[:, 2].astype(str) + '_' +  
                sub_res_in_one_cell.iloc[:, 3].astype(str)  
            )  
            
            # 写入每个基因区域的结果  
            for line in uniq_gene_region:  
                seq = line.split('_')  
                chr = seq[0]  
                start = seq[1]  
                end = seq[2]  
                gene_name = seq[3]  
                res_in_gene = sum(sub_res_in_one_cell[sub_res_in_one_cell.iloc[:, 3] == gene_name].iloc[:, 7])  
                fo.write(chr + '\t' + start + '\t' + end + '\t' + gene_name + '\t' + str(res_in_gene) + '\n')  
            
        # 保存 reads 数据  
        sub_reads_in_one_cell.to_csv('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_ei/20/tmp_one_cell/regular_combined.zz.sorted.modified.'+one_cell, sep='\t', index=False, header=False)  
        
        # 获取 AEI 结果  
        tmp_EI_res = get_AEI('regular', 'gene', one_cell)  
        tmp_EI_res['cell'] = [one_cell] * tmp_EI_res.shape[0]  
        
        return tmp_EI_res

 
    # 假设 cells 是一个包含所有细胞名称的列表  
    fi_total_res=pd.read_table('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_ei/20/tmp_one_cell/geneWithRegular_res.txt',sep='\t',header=None)
    fi_total_reads=pd.read_table('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_ei/20/regular_combined.zz.sorted.modified',sep='\t',header=None)  
    cells=np.unique(fi_total_res.iloc[:,6])
        
    # 初始化空的结果列表  
    all_EI_res = []  
        
    # 使用 ThreadPoolExecutor 运行多线程  
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:  # 你可以调整 max_workers 的数量  
        futures = {executor.submit(process_cell, cell, fi_total_res, fi_total_reads): cell for cell in cells}  
            
        for future in concurrent.futures.as_completed(futures):  
            cell = futures[future]  
            try:  
                tmp_EI_res = future.result()  
                all_EI_res.append(tmp_EI_res)  
                print("Cell %s processing completed!" % cell)  
            except Exception as e:  
                print("Error processing cell %s: %s" % (cell, e))  
        
    # 合并所有结果  
    all_EI_res_df = pd.concat(all_EI_res, ignore_index=True) 
    all_EI_res_df.to_csv('/disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_ei/20/tmp_one_cell/all_EI_res.txt', sep='\t', index=False, header=False)  
            
    subprocess.Popen('rm /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_ei/20/tmp_one_cell/regular_combined.zz.sorted.modified*',shell=True)
    subprocess.Popen('rm /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/gene_ei/20/tmp_one_cell/geneWzRes_ad.*',shell=True)                               
        

    cells=np.unique(all_EI_res_df['cell'])
    genes=np.unique(all_EI_res_df['gene'])
    matrix=np.zeros((len(genes),len(cells)),dtype='float32')
    df=pd.DataFrame(matrix, columns=cells, index=genes) 

    for index, row in all_EI_res_df.iterrows():
        gene=row['gene']
        cell=row['cell']
        df.at[gene,cell]=row['EI']

    sparse_matrix=sparse.csr_matrix(df) 

    genes_str = [str(gene) for gene in genes]
    cells_str = [str(cell) for cell in cells]
    with h5py.File('GEI_sparse_matrix.h5','w') as f:
            sparse_matrix.sort_indices()  # Compatibility with R
            expr = f.create_group("exprs")
            expr.create_dataset("data", data=sparse_matrix.data)
            expr.create_dataset("indices", data=sparse_matrix.indices)
            expr.create_dataset("indptr", data=sparse_matrix.indptr)
            expr.create_dataset("shape", data=sparse_matrix.shape)
            
            gene_name = f.create_group('gene')
            gene_name.create_dataset('gene',data=genes_str)
            barcode_grp = f.create_group('barcode')
            barcode_grp.create_dataset('barcode',data=cells_str)

main() 