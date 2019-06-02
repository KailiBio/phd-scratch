#!/bin/bash

# -- Kaili
# This script is for runing IDEAS for testing imputation between reps with given CTCF data.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 1. rep1-impute-rep1 given 9 CTCF
cp ctcf_samples_rep1-impute-rep1.sh ctcf_samples_rep1-impute-rep1_given9samples.sh
vim ctcf_samples_rep1-impute-rep1_given9samples.sh
cp ctcf_samples_rep1-impute-rep1.parafile ctcf_samples_rep1-impute-rep1_given9samples.parafile
vim ctcf_samples_rep1-impute-rep1_given9samples.parafile
cp ctcf_samples.input ctcf_samples_rep1-impute-rep1_given9samples.input
vim ctcf_samples_rep1-impute-rep1_given9samples.input

nohup bash ctcf_samples_rep1-impute-rep1_given9samples.sh > ./nohup/nohup.ctcf_samples_rep1-impute-rep1_given9samples.out 2>&1&

# 2. rep1-impute-rep2 given 9 CTCF
cp ctcf_samples_rep1-impute-rep2.sh ctcf_samples_rep1-impute-rep2_given9samples.sh
vim ctcf_samples_rep1-impute-rep2_given9samples.sh
cp ctcf_samples_rep1-impute-rep2.parafile ctcf_samples_rep1-impute-rep2_given9samples.parafile
vim ctcf_samples_rep1-impute-rep2_given9samples.parafile
cp ctcf_samples_rep2.input ctcf_samples_rep1-impute-rep2_given9samples.input
vim ctcf_samples_rep1-impute-rep2_given9samples.input

nohup bash ctcf_samples_rep1-impute-rep2_given9samples.sh > ./nohup/nohup.ctcf_samples_rep1-impute-rep2_given9samples.out 2>&1&

# 3. rep1-impute-rep2 given 9 CTCF , 10K
awk '{FS=OFS}{if(NR==1){print $0}else{printf 10000*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 10000*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.para0 > /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples_10k.para0
#
cp ctcf_samples_rep1-impute-rep2_given9samples.sh ctcf_samples_rep1-impute-rep2_given9samples_10k.sh
vim ctcf_samples_rep1-impute-rep2_given9samples_10k.sh
cp ctcf_samples_rep1-impute-rep2_given9samples.parafile ctcf_samples_rep1-impute-rep2_given9samples_10k.parafile
vim ctcf_samples_rep1-impute-rep2_given9samples_10k.parafile
cp ctcf_samples_rep1-impute-rep2_given9samples.input ctcf_samples_rep1-impute-rep2_given9samples_10k.input

nohup bash ctcf_samples_rep1-impute-rep2_given9samples_10k.sh > ./nohup/nohup.ctcf_samples_rep1-impute-rep2_given9samples_10k.out 2>&1&
