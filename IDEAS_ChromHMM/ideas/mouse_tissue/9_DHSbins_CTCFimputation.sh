#!/bin/bash

# -- Kaili
# This script is for using dhs-bins to do CTCF imputation.
# 1. run 9 CTCF sample using dhs-bins (p0)
# 2. impute to 66 ctcf_samples
# 3. use liver14.5 & lung14.5 to validate imputation.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 1. run 9 CTCF sample using dhs-bins
# get input file
awk '{print $1}' /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.input | sort -u > all_ctcf_sample.txt
cp all_ctcf_sample.txt subset_ctcf_sample.txt
vim subset_ctcf_sample.txt
awk '{FS=OFS=" "}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' \
subset_ctcf_sample.txt /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.input > dhs_ctcf.input
sed -i 's/_normal/_dhs/g' dhs_ctcf.input

# get .parafile & .sh file
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.sh dhs_ctcf.sh
vim dhs_ctcf.sh
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.parafile dhs_ctcf.parafile
vim dhs_ctcf.parafile

# other file
cp -r /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/bin ./
cp -r /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/data ./
cp /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_signal_based_space.bed ./

# run
nohup bash dhs_ctcf.sh > nohup.dhs_ctcf.out 2>&1&





# 2. impute to 66 ctcf_samples


# 3. use liver0 & lung0 to validate imputation.
