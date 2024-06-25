from collections import defaultdict
import re

def main():

    gene_records = {}
    current_gene_id = None

    fi=open('FANTOM_CAT.lv1_raw.only_mRNA.gtf')
    for line in fi:
        line=line.rstrip()
        fields = line.rstrip().split('\t')
        
        # 处理gene行
        if fields[2] == 'gene':
            match = re.search(r'gene_id "([^"]+)"', fields[8])
            if match:
                current_gene_id=match.group(1)# 提取gene_id
            gene_records[current_gene_id] = {
                'chrom': fields[0],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'exons': [],
            }
        # 处理transcript和exon行，收集属于当前gene的exon信息
        elif fields[2] in ['transcript', 'exon'] and 'gene_id "' + current_gene_id + '";' in fields[8]:
            if fields[2] == 'exon':
                exon_length = int(fields[4]) - int(fields[3])
                exon_start = int(fields[3]) - gene_records[current_gene_id]['start']
                gene_records[current_gene_id]['exons'].append((exon_length, exon_start))
    fi.close()


    def unique_elements_with_order(lst):
        seen = set()
        result = []
        for item in lst:
            if item not in seen:
                seen.add(item)
                result.append(item)
        return result





    with open('FANTOM.mRNA.psiVer.gtf', 'w') as output_file:
        for gene_id, transcripts in gene_records.iteritems():  # 使用iteritems()兼容Python 2.7
            all_exons = transcripts['exons']
            uniq_all_exons=unique_elements_with_order(all_exons)
            exons_length=[]
            exons_start=[]
            for one in transcripts['exons']:
                exons_length.append(one[0])
                exons_start.append(one[1])
            uniq_exons_length=unique_elements_with_order(exons_length)
            uniq_exons_start=unique_elements_with_order(exons_start)


            exon_lengths_str = ','.join(str(length) for length, _ in uniq_exons_length)
            exon_starts_str = ','.join(str(start) for _, start in uniq_exons_start)

            # 使用老式字符串格式化方法，兼容Python 2.7
            output_line = "%s %d %d %s 0 + %d %d 0 %d %s 0,%s," % (
                gene_records['chrom'], gene_records['start'], gene_records['end'], gene_id,
                gene_records['start'], gene_records['end'], len(unique_exons), exon_lengths_str, exon_starts_str)
            output_file.write(output_line + '\n')

    # 注意：以上代码假定gene_records是一个字典，键是基因ID，值是包含多个转录本记录的列表
    # 每个转录本记录又是一个字典，包含'chrom', 'start', 'end', 'exons'等键
    # 输出重组后的信息
    with open('FANTOM.mRNA.psiVer.gtf', 'w') as output_file:
        for gene_id, record in gene_records.items():
            exon_lengths_str = ','.join(str(length) for length, _ in sorted(record['exons']))
            exon_starts_str = ','.join(str(start) for _, start in sorted(record['exons']))
            output_line = "{0} {1} {2} {3} 0 + {1} {2} 0 {4} {5} 0,{6},".format(
                record['chrom'], record['start'], record['end'], gene_id,
                len(record['exons']), exon_lengths_str, exon_starts_str)
            output_file.write(output_line + '\n')  # 写入每行数据后换行



    for gene_id, record in gene_records.items():
        exon_lengths_str = ','.join(str(length) for length, _ in sorted(record['exons']))
        exon_starts_str = ','.join(str(start) for _, start in sorted(record['exons']))
        print "{0} {1} {2} {3} 0 + {1} {2} 0 {4} {5} 0,{6},".format(
            record['chrom'], record['start'], record['end'], gene_id,
            len(record['exons']), exon_lengths_str, exon_starts_str)


with open('FANTOM_CAT.lv1_raw.only_mRNA.gtf', 'r') as file: 
    content = file.read()
    parse_gtf(content)



chrX\tFANTOM\tgene\t99860347\t99923322\t.\t-\t.\tgene_id "ENSG00000000003.10"; geneSuperClass "all_mRNA";  geneClass "coding_mRNA";  geneSubClass "protein_coding"; gene_type "protein_coding"; gene_name "TSPAN6"; coding_status "coding"; cumulative_support "ENCODE:FANTOM:GENCODE:HUBDMAP:MITRANS"; geneCategory "coding_mRNA"; DHS_type "DHS_promoter";


chrX\tFANTOM\texon\t99880228\t99884983\t.\t-\t.\tgene_id "ENSG00000000003.10"; transcript_id "HBMT00001550123.1"; exon_number 1;



12358-12111