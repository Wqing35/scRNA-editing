##########3组res对应的AS事件的psi如何变化
#直接在已鉴定的psi-change事件里定位编辑位点
#不同事件组成的区间信息不同，分开鉴定
paste <(awk '{if($1!="GW08") print $1}' psi_no_change_events.txt | sed 's/:/\t/g'| sed 's/;/\t/g' ) <(awk '{if($1!="GW08") print $1}' psi_no_change_events.txt) > no_change/psi_no_change_events_for_res_format.txt
cd no_change
#1.RI
#up组：up_res0, down_res0, stable_res0
#down组：up_res0, down_res0, stable_res0
#stable组：up_res0, down_res0, stable_res0
#no-change组：up_res0, down_res0, stable_res0

cat psi_no_change_events_for_res_format.txt| awk '{if($2=="RI") print $3"\t"$4"\t"$6}'|sort|uniq -c |awk '{print $2"\t"$3"\t"$4}'> RI_region_exd2exon.bed

bedtools sort -faidx ~/tmp_data/names.txt -i RI_region_exd2exon.bed > RI_region_exd2exon.srt.bed

bedtools intersect -wo -a RI_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res

bedtools intersect -wo -a RI_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res

bedtools intersect -wo -a RI_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res

#2.SE
#up组：up_res9, down_res3, stable_res0
#down组：up_res0, down_res0, stable_res0
#stable组：up_res62, down_res13, stable_res35
#no-change组：up_res115, down_res15, stable_res101
#——严格的按照skipped exon区间，没有res富集，仍然采用exd2exon的版本
cat psi_no_change_events_for_res_format.txt|awk '{if($2=="SE") print $0}'|sed 's/-/\t/g'|awk '{print $3"\t"$4"\t"$7}'|sort|uniq -c|awk '{print $2"\t"$3"\t"$4}' > SE_region_exd2exon.bed

bedtools sort -faidx ~/tmp_data/names.txt -i SE_region_exd2exon.bed > SE_region_exd2exon.srt.bed

bedtools intersect -wo -a SE_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res

bedtools intersect -wo -a SE_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res

bedtools intersect -wo -a SE_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res

#3.MX
#up组：up_res0, down_res0, stable_res0
#down组：up_res0, down_res0, stable_res0
#stable组：up_res2, down_res2, stable_res5
#no-change组：up_res4, down_res3, stable_res29
cat psi_no_change_events_for_res_format.txt|awk '{if($2=="MX") print $0}'|sed 's/-/\t/g'|awk '{print $3"\t"$4"\t"$11}'|sort|uniq -c|awk '{print $2"\t"$3"\t"$4}' > MX_region_exd2exon.bed

bedtools sort -faidx ~/tmp_data/names.txt -i MX_region_exd2exon.bed > MX_region_exd2exon.srt.bed

bedtools intersect -wo -a MX_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res

bedtools intersect -wo -a MX_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res

bedtools intersect -wo -a MX_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res

rm MX_region_exd2exon.bed

#4.A5——up_res0, down_res0, stable_res0
#down组：up_res0, down_res0, stable_res0
#stable组：up_res4, down_res2, stable_res2
#no-change组：up_res4, down_res4, stable_res2
cat psi_no_change_events_for_res_format.txt|awk '{if($2=="A5") print $0}'|sed 's/-/\t/g'|awk '{print $3"\t"$4"\t"$7}'|sort|uniq -c|awk '{print $2"\t"$3"\t"$4}' > A5_region_exd2exon.bed

bedtools sort -faidx ~/tmp_data/names.txt -i A5_region_exd2exon.bed > A5_region_exd2exon.srt.bed

bedtools intersect -wo -a A5_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res|wc -l

bedtools intersect -wo -a A5_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res|wc -l

bedtools intersect -wo -a A5_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res|wc -l

rm A5_region_exd2exon.bed

#5.A3——up_res0, down_res0, stable_res0
#down组：up_res0, down_res0, stable_res0
#stable组：up_res1, down_res0, stable_res2
#no-change组：up_res9, down_res4, stable_res2
cat psi_no_change_events_for_res_format.txt|awk '{if($2=="A3") print $0}'|sed 's/-/\t/g'|awk '{print $3"\t"$4"\t"$7}'|sort|uniq -c|awk '{print $2"\t"$3"\t"$4}' > A3_region_exd2exon.bed

bedtools sort -faidx ~/tmp_data/names.txt -i A3_region_exd2exon.bed > A3_region_exd2exon.srt.bed

bedtools intersect -wo -a A3_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res|wc -l

bedtools intersect -wo -a A3_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res|wc -l

bedtools intersect -wo -a A3_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res|wc -l

rm A3_region_exd2exon.bed

#6.AF——up_res1, down_res0, stable_res0
#down组：up_res0, down_res0, stable_res0
#stable组：up_res27, down_res13, stable_res14
#no-change组：up_res64, down_res26, stable_res39
cat psi_no_change_events_for_res_format.txt|awk '{if($2=="AF") print $0}'|sed 's/-/\t/g'|awk '{print $3"\t"$4"\t"$9}'|sort|uniq -c|awk '{print $2"\t"$3"\t"$4}' > AF_region_exd2exon.bed

bedtools sort -faidx ~/tmp_data/names.txt -i AF_region_exd2exon.bed > AF_region_exd2exon.srt.bed

bedtools intersect -wo -a AF_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res|wc -l

bedtools intersect -wo -a AF_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res|wc -l

bedtools intersect -wo -a AF_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res|wc -l

rm AF_region_exd2exon.bed

#7.AL——up_res2, down_res0, stable_res1
#down组：up_res0, down_res0, stable_res0
#stable组：up_res13, down_res3, stable_res0
#no-change组：up_res14, down_res2, stable_res29
cat psi_no_change_events_for_res_format.txt|awk '{if($2=="AL") print $0}'|sed 's/-/\t/g'|awk '{print $3"\t"$4"\t"$9}'|sort|uniq -c|awk '{print $2"\t"$3"\t"$4}' > AL_region_exd2exon.bed

bedtools sort -faidx ~/tmp_data/names.txt -i AL_region_exd2exon.bed > AL_region_exd2exon.srt.bed

bedtools intersect -wo -a AL_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res|wc -l

bedtools intersect -wo -a AL_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res|wc -l

bedtools intersect -wo -a AL_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res|wc -l

rm AL_region_exd2exon.bed



###################################先锚定有res在的region，再利用sub psi matrix鉴定差异的as event和是否成组
paste <(awk '{if($1!="GW08") print $1}' MX_psi.txt | sed 's/[:;+-]/\t/g') <(awk '{if($1!="GW08") print $1}' MX_psi.txt) > MX/MX_psi_for_res_format.txt

cd MX

bedtools sort -faidx ~/tmp_data/names.txt -i <(awk '{print $3"\t"$4"\t"$11"\t"$12}' MX_psi_for_res_format.txt|sort|uniq -c |awk '{print $2"\t"$3"\t"$4"\t"$5}') > MX_region_exd2exon.srt.bed

bedtools intersect -wo -a MX_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/up.srt.res |awk '{print $4}'|sort|uniq -c|awk '{print $2}' > MX_region_has_up_res.txt

bedtools intersect -wo -a MX_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/down.srt.res |awk '{print $4}'|sort|uniq -c|awk '{print $2}' > MX_region_has_down_res.txt


bedtools intersect -wo -a MX_region_exd2exon.srt.bed -b /disk1/wenqing/tmp_data/PFC_s2/all_analysis_result/editing_level_change/notwith1/cor_threshold_0.3/stable.srt.res |awk '{print $4}'|sort|uniq -c|awk '{print $2}' > MX_region_has_stable_res.txt


cat MX/MX_region_has_up_res.txt A5/A5_region_has_up_res.txt A3/A3_region_has_up_res.txt AF/AF_region_has_up_res.txt AL/AL_region_has_up_res.txt SE/SE_region_has_up_res.txt RI/RI_region_has_up_res.txt > allType_region_has_up.txt
cat MX/MX_region_has_down_res.txt A5/A5_region_has_down_res.txt A3/A3_region_has_down_res.txt AF/AF_region_has_down_res.txt AL/AL_region_has_down_res.txt SE/SE_region_has_down_res.txt RI/RI_region_has_down_res.txt > allType_region_has_down.txt
cat MX/MX_region_has_stable_res.txt A5/A5_region_has_stable_res.txt A3/A3_region_has_stable_res.txt AF/AF_region_has_stable_res.txt AL/AL_region_has_stable_res.txt SE/SE_region_has_stable_res.txt RI/RI_region_has_stable_res.txt > allType_region_has_stable.txt