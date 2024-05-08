#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import pandas as pd

def main():
    def get_exon_cov(bed_in_dir,bed_out_dir):
        fi=open(bed_in_dir)
        fo=open(bed_out_dir,'w')


        final_cov=[]
        exon_region_cov=[]
        for line in fi:
            seq=line.rstrip().split('\t')
            cov_start=int(seq[3])
            cov_end=int(seq[4])
            exon_region_start=int(seq[1])
            exon_region_end=int(seq[2])
            cov_depth=int(seq[5])
            flag1=cov_start >= exon_region_start and cov_end > exon_region_end
            flag2=cov_start >= exon_region_start and cov_end <= exon_region_end
            flag3=cov_end <= exon_region_end and cov_end <= exon_region_end
            flag4=cov_end <= exon_region_end and cov_end > exon_region_end
            if cov_depth!=0:
                if flag1:
                    new_start=cov_start
                    new_end=exon_region_end
                elif flag2:
                    new_start=cov_start
                    new_end=cov_end
                elif flag3:
                    new_start=exon_region_start
                    new_end=cov_end
                elif flag4:
                    new_start=exon_region_start
                    new_end=exon_region_end
                final_cov.append(new_end-new_start)
                exon_region_cov.append(exon_region_end-exon_region_start)

                fo.write(seq[0]+'\t'+str(new_start)+'\t'+str(new_end)+'\n')

        fi.close()
        fo.close()
        cov_rate=float(sum(final_cov))/sum(exon_region_cov)

        return(cov_rate)

    cov_res_10X=[]
    cov_res_bulk=[]
    cov_res_smart=[]
    gene_names=pd.read_csv("/disk1/wenqing/tmp_data/gene_name.bed",sep='\t')
    for tmp_gene in gene_names:
        input_10X='/disk1/wenqing/tmp_data/coverage_compare/coverage_res/10X/'+tmp_gene+'_depth_clean_10X.exon.bed'
        input_bulk='/disk1/wenqing/tmp_data/coverage_compare/coverage_res/bulk/'+tmp_gene+'_depth_clean_bulk.exon.bed'
        input_smart='/disk1/wenqing/tmp_data/coverage_compare/coverage_res/smart_seq/'+tmp_gene+'_depth_clean_smart.exon.bed'
        cov_res_10X.append(get_exon_cov(input_10X))
        cov_res_bulk.append(get_exon_cov(input_bulk))
        cov_res_smart.append(get_exon_cov(input_smart))

    



        




