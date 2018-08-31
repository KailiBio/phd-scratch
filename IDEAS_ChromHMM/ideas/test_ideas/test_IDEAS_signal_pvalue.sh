#!/bin/bash

# -- Kaili
# This scripts is for testing IDEAS for 66 biosamples (use signal p-value for histone mark data).

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/test_ideas/"
workDir="/data/zusers/fankaili/ideas/run_ideas_p_value/"

cd ${workDir}

# 1. get file from json and metadata
python ${scriptDir}get_ideas_input_hm_pvalue.py
#
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
#
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

# 2. modify all the files

# 3. run IDEAS
nohup bash run_IDEAS_8hm_atac_dname_pvalue.sh > nohup.run_IDEAS_8hm_atac_dname_pvalue.out 2>&1&
