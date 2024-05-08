#!/bin/bash

#phase=('17' '18' '57' '20' '56' '32')
phase='57'
celltype=('L2_3_ExN' 'L4' 'IN_VIP')
# 'ExN' 'L5_6' 'L5_6_CC' 'Olig')

for tmp_phase in ${phase[@]};
do
    for one_type in ${celltype[@]};
    do
        salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/celltype_fq/"$tmp_phase"_"$one_type"_2.fastq -o /disk1/wenqing/tmp_data/ASD/asd_male_pfc//"$tmp_phase"_"$one_type"_transcripts_quant_gencode -p 10
    done
done

#对每一个细胞进行isoform定量
#phase='18'
#barcode=$(cat /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_18/random_sample_25_barcodes.txt)

#for tmp_phase in ${phase[@]};
#do
#    for one_barcode in ${barcode[@]};
#    do
#        salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/"$tmp_phase"_"$one_barcode"_2.fastq -o /disk1/wenqing/tmp_data/ASD/onecell_isoform_quant/"$tmp_phase"_"$one_barcode"_transcripts_quant_gencode -p 2
#    done
#done

#phase='57'
#barcode=$(cat /disk1/wenqing/tmp_data/ASD/asd_male_pfc/asd_57/random_sample_25_barcodes.txt)

#for tmp_phase in ${phase[@]};
#do
#    for one_barcode in ${barcode[@]};
#    do
#        salmon quant -i /disk1/wenqing/tmp_data/hg19/gencode.v44lift37.transcripts.index -l A -r /disk1/wenqing/tmp_data/ASD/"$tmp_phase"_"$one_barcode"_2.fastq -o /disk1/wenqing/tmp_data/ASD/onecell_isoform_quant/"$tmp_phase"_"$one_barcode"_transcripts_quant_gencode -p 2
#    done
#done

