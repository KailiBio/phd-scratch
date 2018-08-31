#!/bin/bash

# -- Kaili
# This script is for repeating yu's IDEAS result.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/test_ideas/"
workDir="/data/zusers/fankaili/ideas/run_ideas_repeat/"

cd ${workDir}

# 1. get file from json and metadata
python ${scriptDir}get_ideas_input_hm_bam.py
#
awk '{FS=OFS="\t"}{if(NR<450){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_1.sh
awk '{FS=OFS="\t"}{if(NR>(1*450-1) && NR<(2*450+1-2)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_2.sh
awk '{FS=OFS="\t"}{if(NR>(2*450-2) && NR<(3*450+1+1)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_3.sh
awk '{FS=OFS="\t"}{if(NR>(3*450-1+2) && NR<(4*450+1)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_4.sh
awk '{FS=OFS="\t"}{if(NR>(4*450) && NR<(5*450+1-1)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_5.sh
awk '{FS=OFS="\t"}{if(NR>(5*450-1) && NR<(6*450+1-2)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_6.sh
awk '{FS=OFS="\t"}{if(NR>(6*450-2) && NR<(7*450+1+1)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_7.sh
awk '{FS=OFS="\t"}{if(NR>(7*450-1+2) && NR<(8*450+1)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_8.sh
awk '{FS=OFS="\t"}{if(NR>(8*450) && NR<(9*450+1)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_9.sh
awk '{FS=OFS="\t"}{if(NR>(9*450)){print $0}}' codes_for_getting_signal_bam.sh > ./code/codes_for_getting_signal_bam_10.sh
#
nohup bash ./code/codes_for_getting_signal_bam_1.sh > ./code/nohup.codes_for_getting_signal_bam_1.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_2.sh > ./code/nohup.codes_for_getting_signal_bam_2.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_3.sh > ./code/nohup.codes_for_getting_signal_bam_3.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_4.sh > ./code/nohup.codes_for_getting_signal_bam_4.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_5.sh > ./code/nohup.codes_for_getting_signal_bam_5.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_6.sh > ./code/nohup.codes_for_getting_signal_bam_6.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_7.sh > ./code/nohup.codes_for_getting_signal_bam_7.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_8.sh > ./code/nohup.codes_for_getting_signal_bam_8.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_9.sh > ./code/nohup.codes_for_getting_signal_bam_9.out 2>&1&
nohup bash ./code/codes_for_getting_signal_bam_10.sh > ./code/nohup.codes_for_getting_signal_bam_10.out 2>&1&


# 2. modify all the files
grep "ATAC" /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input >> run_IDEAS_8hm_atac_dname_repeat.input
grep "DNAme" /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input >> run_IDEAS_8hm_atac_dname_repeat.input


# 3. run IDEAS
nohup bash run_IDEAS_8hm_atac_dname_repeat.sh > nohup.run_IDEAS_8hm_atac_dname_repeat.out 2>&1&
