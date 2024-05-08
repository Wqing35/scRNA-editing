#!/bin/bash

directory=('20' '32' '56')
celltypes=('L2_3_ExN' 'L4' 'L5_6' 'L5_6_CC' 'Olig' 'IN_VIP')

for tmp_directory in ${directory[@]};
do 
    for tmp_celltype in ${celltypes[@]}; 
    do
        echo $tmp_directory
        echo $tmp_celltype
        cd /disk1/wenqing/tmp_data/ASD/ctr_male_pfc/AEI/"$tmp_directory"/celltype_AEI/"$tmp_celltype"
        python2 /disk1/wenqing/SPRINT/SPRINT_master/sprint/get_EI_inGene.py
    done
done


