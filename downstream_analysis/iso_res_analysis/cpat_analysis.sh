make_hexamer_tab.py -c /disk1/wenqing/tmp_data/ASD/IsoformSequence/isoformSwitchAnalyzeR_isoform_nt.fasta   -n ../gencode.v44lift37.lncRNA_transcripts.fa > Human_Hexamer.tsv

make_logitModel.py  -x Human_Hexamer.tsv -c /disk1/wenqing/tmp_data/ASD/IsoformSequence/isoformSwitchAnalyzeR_isoform_nt.fasta -n ../gencode.v44lift37.lncRNA_transcripts.fa -o Human

cpat.py -g /disk1/wenqing/tmp_data/ASD/IsoformSequence/isoformSwitchAnalyzeR_isoform_nt.fasta  -d ./Human_logitModel.RData -x ./Human_Hexamer.tsv -o output1



bedtools intersect -wo -a /disk1/wenqing/tmp_data/gene.bed -b resInAluWithADDP_over0| awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"_"$7"\t"$9}'|more

cat regular_combined.zz.sorted.modified|awk '{if($2 > 65192179){print $0}else if($3 < 65192438){print $0}else{print "pass this read"}}'|more



wget http://eddylab.org/software/hmmer/hmmer-3.2.tar.gz
tar -xzvf  hmmer-3.2.1.tar.gz
cd hmmer-3.2
./configure
make
make check
make install 
vim ~/.bashrc#添加环境变量
export PATH=/usr/local/bin:$PATH
source ~/.bashrc#激活环境变量


wget ftp://ftp.ebi.ac.uk:21/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk:21/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget ftp://ftp.ebi.ac.uk:21/pub/databases/Pfam/current_release/active_site.dat.gz
gunzip *.gz

hmmpress Pfam-A.hmm#建库

pfam_scan.pl -fasta /disk1/wenqing/tmp_data/ASD/IsoformSequence/isoformSwitchAnalyzeR_isoform_AA.fasta -dir /disk1/wenqing/tmp_data/ASD/prepForIso/pfam_hmm -outfile /disk1/wenqing/tmp_data/ASD/prepForIso/pfam_res/pfam_scan_result.fa -as