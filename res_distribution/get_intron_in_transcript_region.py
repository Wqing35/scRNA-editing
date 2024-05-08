#!/disk1/wenqing/anaconda3/envs/wq_py2/bin/python2
# coding=utf-8

import pandas as pd  
  
# 读取转录本的BED文件  
transcripts_bed = pd.read_csv('/disk1/wenqing/tmp_data/my_transcript_wzGeneName.bed', sep='\t', header=None, names=['chrom', 'start', 'end', 'transcript_name'])  
  
# 读取内含子的BED文件  
introns_bed = pd.read_csv('/disk1/wenqing/tmp_data/my_intron_wzGeneName.bed', sep='\t', header=None, names=['chrom', 'start', 'end', 'intron_name'])  
  
# 创建一个空DataFrame来存储结果  
results = pd.DataFrame(columns=['chrom', 'transcript_start', 'transcript_end', 'transcript_name', 'intron_starts', 'intron_ends'])  
  
# 遍历每个转录本：196520次迭代
for index, transcript in transcripts_bed.iterrows():  
    # 提取转录本的染色体、起始位置、终止位置和名称  
    chrom = transcript['chrom']  
    transcript_start = transcript['start']  
    transcript_end = transcript['end']  
    transcript_name = transcript['transcript_name']  
      
    # 查找该转录本内的所有内含子  
    introns_in_transcript = introns_bed[(introns_bed['chrom'] == chrom) &   
                             (introns_bed['start'] >= transcript_start) &   
                             (introns_bed['end'] <= transcript_end)]  
      
    # 如果存在内含子，则提取它们的起始和终止位置  
    if not introns_in_transcript.empty:  
        intron_starts = ','.join(map(str, introns_in_transcript['start'].tolist()))  
        intron_ends = ','.join(map(str, introns_in_transcript['end'].tolist()))  
    else:  
        intron_starts = 'None'  
        intron_ends = 'None'  
      
    # 将结果添加到DataFrame中  
    results = results.append({'chrom': chrom,   
                               'transcript_start': transcript_start,   
                               'transcript_end': transcript_end,   
                               'transcript_name': transcript_name,   
                               'intron_starts': intron_starts,   
                               'intron_ends': intron_ends},   
                              ignore_index=True)  
  
# 将结果保存到新的BED文件中  
results.to_csv('/disk1/wenqing/tmp_data/gencode_hg19_intron_in_transcript.txt', sep='\t', index=False, header=False)