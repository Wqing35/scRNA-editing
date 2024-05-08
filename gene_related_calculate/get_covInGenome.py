#!/usr/bin/python2
# coding=utf-8


def main():

    def get_covInChr(bed_in_dir=0):
        fi=open(bed_in_dir)
        all_cov=0
        whole_len=0
        for line in fi:
            seq=line.rstrip().split('\t')
            if int(eval(seq[3]))!=0:
                cov=int(seq[2])-int(seq[1])
            else:
                cov=0
            tmp_len=int(seq[2])-int(seq[1])
            all_cov=all_cov+cov
            whole_len=whole_len+tmp_len 
        fi.close()       
        ave_cov=float(all_cov)/whole_len
        return(ave_cov)

    input_path='/disk1/wenqing/tmp_data/coverage_compare/'
    #output='/local/wenqing/data/geneElements/results/GW16/'
    tmp_aver=[]
    chrs='chr1'
    #chrs='chr16'
    #element=['exon']
    for chr in chrs:
        input=input_path+'/RBM8A_depth_clean_bulk.bed'
        tmp_aver.append(get_covInChr(input))
    aver_covInGenome=sum(set(tmp_aver))/len(chrs)
    print(aver_covInGenome)



main()
