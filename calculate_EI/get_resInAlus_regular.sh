cat regular_merged.res.anno|grep Alu |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > regular_merged.res.alu.modified
bedtools intersect -wb -a regular_merged.res.alu.modified -b ~/tmp_data/Alu_seq/AsInAlus_modified.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$12"_"$13}' > resInAluWithADDP

cat resInAluWithADDP |awk '{print $5}'|awk -F ':' '{print $2}' > dp.txt
paste -d '\t' resInAluWithADDP dp.txt > resInAluWithADDP_WizDepth
