#!/bin/bash

# -- Kaili
# This script is for imputing using 11sample trained model.

mkdir /data/zusers/fankaili/ideas/CTCF_impute/imputation/

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/imputation/"

cd ${workDir}

# 1. first try
## set otherparap and imputation="CTCF"
### .input file
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.input imputation1.input
cat /data/zusers/fankaili/ideas/dhs_bins/normal_bins/66samples_10marks_normal_bins.input imputation1.input > tmp.input
sort -u  tmp.input > imputation1.input
rm tmp.input
### .parafile
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.parafile imputation1.parafile
vim imputation1.parafile
### .sh
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples.sh imputation1.sh
vim imputation1.sh
### other
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/mm10.bed ./
cp -r /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/bin ./
cp -r /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/data ./
### run
nohup bash imputation1.sh > nohup.imputation1.out 2>&1&
