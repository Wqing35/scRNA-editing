#!/

import re
from collections import OrderedDict
import numpy as np

# 指定输入的GTF文件和输出文件名
input_file = "/disk1/wenqing/tmp_data/hg19.refGene.gtf"
output_file = "/disk1/wenqing/tmp_data/hg19.refGene.processed.txt"

fi=open(input_file)

genes_dict=OrderedDict()
for line in fi:
    seq=line.rstrip().split('\t')
    gene_info=seq[8]
    gene_name=gene_info.rstrip().split(' ')[1]
    gene_name=re.sub(r'"', '', gene_name)
    gene_name=re.sub(r';', '', gene_name)
    try:
        genes_dict[gene_name].append([seq[0]+'\t'+seq[2]+'\t'+seq[3]+'\t'+seq[4]])
    except Exception as e:
        genes_dict[gene_name]=[seq[0]+'\t'+seq[2]+'\t'+seq[3]+'\t'+seq[4]]
fi.close()

tmp_key=''
subset_keys=[]
for one in genes_dict.keys():
    for tmp_ele in  genes_dict[one]:
        new_str=''.join(tmp_ele)
        if '3UTR' in new_str:
            tmp_key=one
        else:
            continue
    subset_keys.append(tmp_key)

subset_keys=np.unique(subset_keys)

subset_dict_3UTR = OrderedDict((key, genes_dict[key]) for key in subset_keys if key in genes_dict)

fo = open(output_file, 'w')
# 提取基因和对应的3UTR区域信息
for one in subset_dict_3UTR.keys():
    gene_info = []
    utr_info = []
    for tmp_ele in subset_dict_3UTR[one]:
        region_info = ''.join(tmp_ele)
        if 'transcript' in region_info:
            seq = region_info.rstrip().split('\t')
            gene_start = seq[2]
            gene_end = seq[3]
            gene_info.append((gene_start, gene_end, one))
        if '3UTR' in region_info:
            seq = region_info.rstrip().split('\t')
            utr_start = seq[2]
            utr_end = seq[3]
            utr_info.append((utr_start, utr_end))

    # 写入所有 transcript 和 UTR 区域信息
    for gene_start, gene_end, transcript in gene_info:
        for utr_start, utr_end in utr_info:
            fo.write(seq[0] + '\t' + gene_start + '\t' + gene_end + '\t' + transcript + '\t' + utr_start + '\t' + utr_end + '\n')
    
    fo.write('\n')  # 每个基因信息之间添加空行

fo.close()

#将上述文件按3’UTR整理成chr名称、（trancript的起点、3UTR的起点）/（3UTR的起点、转录本的终点）、转录本的名字

output_file_5 = '/disk1/wenqing/tmp_data/gtf_files/refseq/transcriptsTabByUTR_5end_sorted.txt'  
output_file_3 = '/disk1/wenqing/tmp_data/gtf_files/refseq/transcriptsTabByUTR_3end_sorted.txt'  
  
# 初始化一个列表来保存筛选后的行  
filtered_lines_5 = []  
filtered_lines_3 = []  
  
# 初始化一个字典来跟踪已经写入的基因名称和对应的数字后缀  
gene_count = {}  
  
# 打开输入文件  
with open(output_file, 'r') as fi_1:  
    for line in fi_1:  
        if line.strip() != '':  
            seq = line.rstrip().split('\t')  
            s = seq[0]  
            transcript_start = seq[1]  
            transcript_end = seq[2]  
            gene_name = seq[3]  
            utr_start = seq[4]  
            utr_end = seq[5]  
              
            if s.startswith('chr') and (len(s) <= 5 and (s[3] in ['X', 'Y'] or s[3:].isdigit())):  
                # 检查基因名称是否已经写入过  
                if gene_name in gene_count:  
                    gene_name_with_suffix = "{}_{}".format(gene_name, gene_count[gene_name])  
                    gene_count[gene_name] += 1  
                else:  
                    gene_name_with_suffix = gene_name  
                    gene_count[gene_name] = 1  
                  
                if transcript_start < utr_start and transcript_end > utr_end:  
                    # 分别添加到两个列表，用于不同的输出文件  
                    filtered_lines_5.append("{}\t{}\t{}\t{}\n".format(s, transcript_start, utr_start, gene_name_with_suffix))  
                    filtered_lines_3.append("{}\t{}\t{}\t{}\n".format(s, int(utr_start) + 1, transcript_end, gene_name_with_suffix))  
  
# 对筛选后的行进行排序  
# Python 2.7 需要手动指定比较函数的参数类型  
def sort_key(line):  
    seq = line.split('\t')  
    s = seq[0]  
    return (int(s[3:]) if s[3:].isdigit() else (ord(s[3]), s))  
  
# 使用 sorted 函数进行排序，并指定 key 函数  
sorted_lines_5 = sorted(filtered_lines_5, key=sort_key)  
sorted_lines_3 = sorted(filtered_lines_3, key=sort_key)  
  
# 打开输出文件并写入排序后的行  
with open(output_file_5, 'w') as fo_1:  
    fo_1.writelines(sorted_lines_5)  
  
with open(output_file_3, 'w') as fo_2:  
    fo_2.writelines(sorted_lines_3) 
##hg19_processed.txt有三个版本，0：保留了最后一个transcript和utr；1:保留第一个出现的transcipt和第一个出现的utr；2:如果两者存在多对一的关系，均保留（transcript区间为唯一标识符）
##注：是否下载UCSC的基因区域信息，将transcript的其实位置转化为gene的位置（统一化的操作），这样做上一步区分几个txt文件版本没啥意义




#######################Ver1:只保留第一个出现的transcript的起始位置和最后一个UTR的起始位置
fo = open(output_file, 'w')    
# 遍历subset_dict_3UTR字典中的每个基因  
for gene, regions in subset_dict_3UTR.items():  
    # 初始化变量来存储第一个transcript和最后一个3'UTR的信息  
    first_transcript_info = None  
    last_utr_info = None  
      
    # 遍历该基因的所有区域信息  
    for region_element in regions:  
        region_info = ''.join(region_element)  
        seq = region_info.rstrip().split('\t')  
        if 'transcript' in region_info:  
            gene_start = seq[2]  
            gene_end = seq[3]  
              
            # 如果first_transcript_info还没有被设置，则设置它  
            if first_transcript_info is None:  
                first_transcript_info = (gene_start, gene_end, seq[0])  # 假设seq[0]是transcript ID  
          
        elif '3UTR' in region_info:  
            utr_start = seq[2]  
            utr_end = seq[3]  
              
            # 更新last_utr_info  
            last_utr_info = (utr_start, utr_end)  
      
    # 检查是否找到了transcript和3'UTR信息  
    if first_transcript_info and last_utr_info:  
        gene_start, gene_end, transcript = first_transcript_info  
        utr_start, utr_end = last_utr_info  
          
        # 写入第一个transcript和最后一个3'UTR的信息  
        fo.write(transcript + '\t' + gene_start + '\t' + gene_end + '\t' + gene + '\t' + utr_start + '\t' + utr_end + '\n')  
      
    fo.write('\n')  # 每个基因信息之间添加空行  
  
fo.close()


with open(output_file, 'r') as fi_1:  
    for line in fi_1:  
        if line.strip() != '':  
            seq = line.rstrip().split('\t')  
            s = seq[0]  
            transcript_start = seq[1]  
            transcript_end = seq[2]  
            gene_name = seq[3]  
            utr_start = seq[4]  
            utr_end = seq[5]  
              
            if s.startswith('chr') and (len(s) <= 5 and (s[3] in ['X', 'Y'] or s[3:].isdigit())):                    
                if transcript_start < utr_start and transcript_end > utr_end:  
                    # 分别添加到两个列表，用于不同的输出文件  
                    filtered_lines_5.append("{}\t{}\t{}\t{}\n".format(s, transcript_start, utr_start, gene_name))  
                    filtered_lines_3.append("{}\t{}\t{}\t{}\n".format(s, int(utr_start) + 1, transcript_end, gene_name))  
##################################Ver1#####################################




