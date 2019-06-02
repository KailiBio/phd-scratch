#!/bin/bash

# -- Kaili
# This script is for running IDEAS using all 66 samples, with or without CTCF.
# To see how many new states will generate.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 1. all 66 samples with CTCF
## get input file
cat ctcf_9sample_impute.input ctcf_9sample_to_11sample.input | sort -u > all_data_impuatation.input
## get sh file
cp ctcf_8sample_impute_11sample.sh all_data_impuatation.sh
vim all_data_impuatation.sh
## get parafile
cp ctcf_8sample_impute_11sample.parafile all_data_impuatation.parafile
vim all_data_impuatation.parafile

nohup bash all_data_impuatation.sh > ./nohup/nohup.all_data_impuatation.out 2>&1&
# z001 22638
# get 44 states in total


# 2. all 66 samples without CTCF, set state number the same as the last run
## get input file
awk '{FS=OFS}{if($2!="CTCF"){print $0}}' all_data_impuatation.input > all_data_44states.input
## get sh file
cp all_data_impuatation.sh all_data_44states.sh
vim all_data_44states.sh
##  get parafile
# set state num=43
cp all_data_impuatation.parafile all_data_44states.parafile
vim all_data_44states.parafile

nohup bash all_data_44states.sh > ./nohup/nohup.all_data_44states.out 2>&1&
# z001 28615
# z001 28965


# 3. make hierarchical clutering, look at state changes
# Rscript make_hclust_after_imputate_CTCF.R
