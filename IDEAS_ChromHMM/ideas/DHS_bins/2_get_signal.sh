#!/bin/bash

# -- Kaili
# This script is for getting MACS2 p-value signal for given DHS-center bins file.
# 1. v1_100_300bp
# 2. v2_1_300bp

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/"

cd ${workDir}

# 1. v1_100_300bp
sort -k1,1 mm10_OCR-center_bins_v1.bed > mm10_OCR-center_bins_v1.sorted.bed
python ${scriptDir}get_66samples_10marks_signal_pvalue.py /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v1.sorted.bed \
/data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/
#
cd ${workDir}v1_100_300bp/
mkdir signal
mkdir code
# divide code file
awk '{FS=OFS="\t"}{if(NR<133){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_1.sh
awk '{FS=OFS="\t"}{if(NR>(1*132) && NR<(2*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_2.sh
awk '{FS=OFS="\t"}{if(NR>(2*132) && NR<(3*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_3.sh
awk '{FS=OFS="\t"}{if(NR>(3*132) && NR<(4*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_4.sh
awk '{FS=OFS="\t"}{if(NR>(4*132) && NR<(5*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_5.sh
awk '{FS=OFS="\t"}{if(NR>(5*132) && NR<(6*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_6.sh
awk '{FS=OFS="\t"}{if(NR>(6*132) && NR<(7*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_7.sh
awk '{FS=OFS="\t"}{if(NR>(7*132) && NR<(8*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_8.sh
awk '{FS=OFS="\t"}{if(NR>(8*132) && NR<(9*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_9.sh
awk '{FS=OFS="\t"}{if(NR>(9*132) && NR<(10*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_10.sh
# run code to get signal file
nohup bash ./code/codes_for_getting_signal_pvalue_1.sh > ./code/nohup.codes_for_getting_signal_pvalue_1.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_2.sh > ./code/nohup.codes_for_getting_signal_pvalue_2.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_3.sh > ./code/nohup.codes_for_getting_signal_pvalue_3.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_4.sh > ./code/nohup.codes_for_getting_signal_pvalue_4.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_5.sh > ./code/nohup.codes_for_getting_signal_pvalue_5.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_6.sh > ./code/nohup.codes_for_getting_signal_pvalue_6.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_7.sh > ./code/nohup.codes_for_getting_signal_pvalue_7.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_8.sh > ./code/nohup.codes_for_getting_signal_pvalue_8.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_9.sh > ./code/nohup.codes_for_getting_signal_pvalue_9.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_10.sh > ./code/nohup.codes_for_getting_signal_pvalue_10.out 2>&1&

# .bed file
awk '{FS="\t";OFS=" "}{print $1,$2,$3,$4}' ${workDir}mm10_OCR-center_bins_v1.sorted.bed > mm10_OCR-center_bins_v1.sorted_space.bed
awk '{FS="\t";OFS=" "}{if(NR==FNR){a[$4]=$0}else{print a[$1]}}' /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v1.bed \
/data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal/heart_0_ATAC.tab > mm10_OCR-center_bins_v1_signal_based.bed
awk '{FS="\t";OFS=" "}{print $1,$2,$3,$4}' mm10_OCR-center_bins_v1_signal_based.bed > mm10_OCR-center_bins_v1_signal_based_space.bed
# .input file
cp signal_pvalue.input DHS_v1_100-300bp.input
# .sh file
cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.sh DHS_v1_100-300bp.sh
vim DHS_v1_100-300bp.sh
# .parafile file
cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.parafile DHS_v1_100-300bp.parafile
vim DHS_v1_100-300bp.parafile

# run IDEAS
nohup bash DHS_v1_100-300bp.sh > ./nohup.DHS_v1_100-300bp.out 2>&1&


# 2. v2_1_300bp
sort -k1,1 mm10_OCR-center_bins_v2.bed > mm10_OCR-center_bins_v2.sorted.bed
python ${scriptDir}get_66samples_10marks_signal_pvalue.py /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v2.sorted.bed \
/data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/
#
cd ${workDir}v2_1_300bp/
mkdir signal
mkdir code
# divide code file
awk '{FS=OFS="\t"}{if(NR<133){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_1.sh
awk '{FS=OFS="\t"}{if(NR>(1*132) && NR<(2*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_2.sh
awk '{FS=OFS="\t"}{if(NR>(2*132) && NR<(3*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_3.sh
awk '{FS=OFS="\t"}{if(NR>(3*132) && NR<(4*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_4.sh
awk '{FS=OFS="\t"}{if(NR>(4*132) && NR<(5*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_5.sh
awk '{FS=OFS="\t"}{if(NR>(5*132) && NR<(6*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_6.sh
awk '{FS=OFS="\t"}{if(NR>(6*132) && NR<(7*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_7.sh
awk '{FS=OFS="\t"}{if(NR>(7*132) && NR<(8*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_8.sh
awk '{FS=OFS="\t"}{if(NR>(8*132) && NR<(9*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_9.sh
awk '{FS=OFS="\t"}{if(NR>(9*132) && NR<(10*132+1)){print $0}}' codes_for_getting_signal_pvalue.sh > ./code/codes_for_getting_signal_pvalue_10.sh
# run code to get signal file
nohup bash ./code/codes_for_getting_signal_pvalue_1.sh > ./code/nohup.codes_for_getting_signal_pvalue_1.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_2.sh > ./code/nohup.codes_for_getting_signal_pvalue_2.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_3.sh > ./code/nohup.codes_for_getting_signal_pvalue_3.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_4.sh > ./code/nohup.codes_for_getting_signal_pvalue_4.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_5.sh > ./code/nohup.codes_for_getting_signal_pvalue_5.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_6.sh > ./code/nohup.codes_for_getting_signal_pvalue_6.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_7.sh > ./code/nohup.codes_for_getting_signal_pvalue_7.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_8.sh > ./code/nohup.codes_for_getting_signal_pvalue_8.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_9.sh > ./code/nohup.codes_for_getting_signal_pvalue_9.out 2>&1&
nohup bash ./code/codes_for_getting_signal_pvalue_10.sh > ./code/nohup.codes_for_getting_signal_pvalue_10.out 2>&1&

# .bed file
awk '{FS="\t";OFS=" "}{print $1,$2,$3,$4}' ${workDir}mm10_OCR-center_bins_v2.sorted.bed > mm10_OCR-center_bins_v2.sorted_space.bed
awk '{FS="\t";OFS=" "}{if(NR==FNR){a[$4]=$0}else{print a[$1]}}' /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v2.bed \
/data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal/heart_0_DNAme.tab > mm10_OCR-center_bins_v2_signal_based.bed
awk '{FS="\t";OFS=" "}{print $1,$2,$3,$4}' mm10_OCR-center_bins_v2_signal_based.bed > mm10_OCR-center_bins_v2_signal_based_space.bed
# .input file
cp signal_pvalue.input DHS_v2_1-300bp.input
# .sh file
cp /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp.sh DHS_v2_1-300bp.sh
vim DHS_v2_1-300bp.sh
# .parafile file
cp /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp.parafile DHS_v2_1-300bp.parafile
vim DHS_v2_1-300bp.parafile

# run IDEAS
nohup bash DHS_v2_1-300bp.sh > ./nohup.DHS_v2_1-300bp.out 2>&1&
