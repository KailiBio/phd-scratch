#!/bin/bash

# -- Kaili
# This script is for testing imputation by adding rep2.
# 1. impute rep2 by rep1 (11sample, all have CTCF)
# 2. run rep2 only
# 3. rep1-imput-rep2 one-sample-only
# 4.  rep1 replicates
# 5. rep1-impute-rep1
# 6. rep1-impute-rep1 one sample

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 1. impute rep2 by rep1 (11sample, all have CTCF)
## 1) .sh
cp ctcf_samples.sh ctcf_samples_rep1-impute-rep2.sh
vim ctcf_samples_rep1-impute-rep2.sh
## 2) .parafile
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.para0 > /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples_100fold.para0
#
cp ctcf_samples.parafile ctcf_samples_rep1-impute-rep2.parafile
vim ctcf_samples_rep1-impute-rep2.parafile
## 3) .input
cp ctcf_samples.input ctcf_samples_rep1-impute-rep2.input
sed -i 's/rep1/rep2/g' ctcf_samples_rep1-impute-rep2.input
vim ctcf_samples_rep1-impute-rep2.input

nohup bash ctcf_samples_rep1-impute-rep2.sh > ./nohup/nohup.ctcf_samples_rep1-impute-rep2.out 2>&1&


# 2. run rep2 only
cp ctcf_samples_rep1-impute-rep2.sh ctcf_samples_rep2.sh
vim ctcf_samples_rep2.sh
cp ctcf_samples_rep1-impute-rep2.parafile ctcf_samples_rep2.parafile
vim ctcf_samples_rep2.parafile
cp ctcf_samples.input ctcf_samples_rep2.input
sed -i 's/rep1/rep2/g' ctcf_samples_rep2.input

nohup bash ctcf_samples_rep2.sh > ./nohup/nohup.ctcf_samples_rep2.out 2>&1&

# 3. rep1-imput-rep2 one-sample-only
cp ctcf_samples_rep1-impute-rep2.sh ctcf_samples_rep1-impute-rep2_oneSample.sh
vim ctcf_samples_rep1-impute-rep2_oneSample.sh
cp ctcf_samples_rep1-impute-rep2.parafile ctcf_samples_rep1-impute-rep2_oneSample.parafile
vim ctcf_samples_rep1-impute-rep2_oneSample.parafile
grep "lung_14.5" ctcf_samples_rep1-impute-rep2.input > ctcf_samples_rep1-impute-rep2_oneSample.input

nohup bash ctcf_samples_rep1-impute-rep2_oneSample.sh > ./nohup/nohup.ctcf_samples_rep1-impute-rep2_oneSample.out 2>&1&

# 4.  rep1 replicates
cp ctcf_samples.input ctcf_samples_rep1_2.input
cp ctcf_samples.sh ctcf_samples_rep1_2.sh
vim ctcf_samples_rep1_2.sh
cp ctcf_samples.parafile ctcf_samples_rep1_2.parafile
vim ctcf_samples_rep1_2.parafile

nohup bash ctcf_samples_rep1_2.sh > ./nohup/nohup.ctcf_samples_rep1_2.out 2>&1&

# 5. rep1-impute-rep1
cp ctcf_samples.input ctcf_samples_rep1-impute-rep1.input
vim ctcf_samples_rep1-impute-rep1.input
cp ctcf_samples_rep1-impute-rep2.sh ctcf_samples_rep1-impute-rep1.sh
vim ctcf_samples_rep1-impute-rep1.sh
cp ctcf_samples_rep1-impute-rep2.parafile ctcf_samples_rep1-impute-rep1.parafile
vim ctcf_samples_rep1-impute-rep1.parafile

nohup bash ctcf_samples_rep1-impute-rep1.sh > ./nohup/nohup.ctcf_samples_rep1-impute-rep1.out 2>&1&

# 6. rep1-impute-rep1 one sample
grep "lung_14.5" ctcf_samples_rep1-impute-rep1.input > ctcf_samples_rep1-impute-rep1_oneSample.input
cp ctcf_samples_rep1-impute-rep2_oneSample.sh ctcf_samples_rep1-impute-rep1_oneSample.sh
vim ctcf_samples_rep1-impute-rep1_oneSample.sh
cp ctcf_samples_rep1-impute-rep2_oneSample.parafile ctcf_samples_rep1-impute-rep1_oneSample.parafile
vim ctcf_samples_rep1-impute-rep1_oneSample.parafile

nohup bash ctcf_samples_rep1-impute-rep1_oneSample.sh > ./nohup/nohup.ctcf_samples_rep1-impute-rep1_oneSample.out 2>&1&
