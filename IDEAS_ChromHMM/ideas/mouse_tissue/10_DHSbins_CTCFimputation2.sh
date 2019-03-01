#!/bin/bash

# -- Kaili
# This script is for using dhs-bins to do CTCF imputation.
# impute time point
# 1. use 9 sample p0 model to impute e14.5 CTCF
# 2. use 9 sample p0 model to run all 11 samples
#  impute celltype
# 3. get 8 sample model (without brain tissue)
# 4. use 8 sample without_brain model to impute brain CTCF
# 5. use 8 sample without_brain model to run all 11 samples
### all have two runs

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}


# 1. use 9 sample p0 model to impute e14.5 CTCF
## get input file
awk '{FS=OFS=" "}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' \
all_66_sample.txt /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.input > ctcf_9sample_to_11sample.input
sed -i 's/_normal/_dhs/g' ctcf_9sample_to_11sample.input
#
cp ctcf_9sample_to_11sample.input ctcf_9sample_impute_11sample.input
vim ctcf_9sample_impute_11sample.input
## get sh file
cp ctcf_9sample_impute.sh ctcf_9sample_impute_11sample.sh
vim ctcf_9sample_impute_11sample.sh
## get parafile
cp ctcf_9sample_impute.parafile ctcf_9sample_impute_11sample.parafile
vim ctcf_9sample_impute_11sample.parafile

nohup bash ctcf_9sample_impute_11sample.sh > ./nohup/nohup.ctcf_9sample_impute_11sample.out 2>&1&
# z008 53820

#----------
# second run
cp ctcf_9sample_impute_11sample.input ctcf_9sample_impute_11sample_2.input
cp ctcf_9sample_impute_11sample.sh ctcf_9sample_impute_11sample_2.sh
vim ctcf_9sample_impute_11sample_2.sh
cp ctcf_9sample_impute_11sample.parafile ctcf_9sample_impute_11sample_2.parafile
vim ctcf_9sample_impute_11sample_2.parafile
nohup bash ctcf_9sample_impute_11sample_2.sh > ./nohup/nohup.ctcf_9sample_impute_11sample_2.out 2>&1&
# z008 45675

# 2. use 9 sample p0 model to run all 11 samples
## get sh file
cp ctcf_9sample_impute.sh ctcf_9sample_to_11sample.sh
vim ctcf_9sample_to_11sample.sh
## get parafile
cp ctcf_9sample_impute.parafile ctcf_9sample_to_11sample.parafile
vim ctcf_9sample_to_11sample.parafile

nohup bash ctcf_9sample_to_11sample.sh > ./nohup/nohup.ctcf_9sample_to_11sample.out 2>&1&
# z010 9052

#----------
# second run
cp ctcf_9sample_to_11sample.input ctcf_9sample_to_11sample_2.input
cp ctcf_9sample_to_11sample.sh ctcf_9sample_to_11sample_2.sh
vim ctcf_9sample_to_11sample_2.sh
cp ctcf_9sample_to_11sample.parafile ctcf_9sample_to_11sample_2.parafile
vim ctcf_9sample_to_11sample_2.parafile
nohup bash ctcf_9sample_to_11sample_2.sh > ./nohup/nohup.ctcf_9sample_to_11sample_2.out 2>&1&
# z003 34734

# 3. get 8 sample model (without brain tissue)
## get input file
cp all_ctcf_sample.txt subset_ctcf_sample2.txt
vim subset_ctcf_sample2.txt
awk '{FS=OFS=" "}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' \
subset_ctcf_sample2.txt /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.input > ctcf_8sample_model.input
sed -i 's/_normal/_dhs/g' ctcf_8sample_model.input
## get sh file
cp dhs_ctcf.sh ctcf_8sample_model.sh
vim ctcf_8sample_model.sh
## get parafile
cp dhs_ctcf.parafile ctcf_8sample_model.parafile
vim ctcf_8sample_model.parafile

nohup bash ctcf_8sample_model.sh > ./nohup/nohup.ctcf_8sample_model.out 2>&1&
# z010 11707


# 4. use 8 sample without_brain model to impute brain CTCF
# get 100 fold para0 file
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model_result/ctcf_8sample_model.para0 > \
/data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model_result/ctcf_8sample_model_100fold.para0
# get input file
cp ctcf_9sample_to_11sample.input ctcf_8sample_impute_11sample.input
vim ctcf_8sample_impute_11sample.input
## get sh file
cp ctcf_9sample_impute.sh ctcf_8sample_impute_11sample.sh
vim ctcf_8sample_impute_11sample.sh
## get parafile
cp ctcf_9sample_impute.parafile ctcf_8sample_impute_11sample.parafile
vim ctcf_8sample_impute_11sample.parafile

nohup bash ctcf_8sample_impute_11sample.sh > ./nohup/nohup.ctcf_8sample_impute_11sample.out 2>&1&
# z010 7457

#----------
# second run
cp ctcf_8sample_impute_11sample.input ctcf_8sample_impute_11sample_2.input
cp ctcf_8sample_impute_11sample.sh ctcf_8sample_impute_11sample_2.sh
vim ctcf_8sample_impute_11sample_2.sh
cp ctcf_8sample_impute_11sample.parafile ctcf_8sample_impute_11sample_2.parafile
vim ctcf_8sample_impute_11sample_2.parafile
nohup bash ctcf_8sample_impute_11sample_2.sh > ./nohup/nohup.ctcf_8sample_impute_11sample_2.out 2>&1&


# 5. use 8 sample without_brain model to run all 11 samples
# get input file
cp ctcf_9sample_to_11sample.input ctcf_8sample_to_11sample.input
## get sh file
cp ctcf_9sample_to_11sample.sh ctcf_8sample_to_11sample.sh
vim ctcf_8sample_to_11sample.sh
## get parafile
cp ctcf_9sample_to_11sample.parafile ctcf_8sample_to_11sample.parafile
vim ctcf_8sample_to_11sample.parafile

nohup bash ctcf_8sample_to_11sample.sh > ./nohup/nohup.ctcf_8sample_to_11sample.out 2>&1&
# z010 7565

#----------
# second run
cp ctcf_8sample_to_11sample.input ctcf_8sample_to_11sample_2.input
cp ctcf_8sample_to_11sample.sh ctcf_8sample_to_11sample_2.sh
vim ctcf_8sample_to_11sample_2.sh
cp ctcf_8sample_to_11sample.parafile ctcf_8sample_to_11sample_2.parafile
vim ctcf_8sample_to_11sample_2.parafile
nohup bash ctcf_8sample_to_11sample_2.sh > ./nohup/nohup.ctcf_8sample_to_11sample_2.out 2>&1&
