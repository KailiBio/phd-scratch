#!/bin/bash

# -- Kaili
# This script is for record ideas imputation runs.
# 1. using E14.5&P0, run IDEAS on 8HM+ATAC+DNAme+CTCF
## 1) with all tissues in e14.5&P0
## 2) Only use tissues with CTCF in e14.5&P0

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/"

cd ${workDir}

# 1. using E14.5&P0, run IDEAS on 8HM+ATAC+DNAme+CTCF
mkdir e14.5p0_8hm_ATAC_DNAme_CTCF
cd ${workDir}e14.5p0_8hm_ATAC_DNAme_CTCF


## 1) with all tissues in e14.5&P0
mkdir e14.5p0_8hm_ATAC_DNAme_CTCF_1
cp /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS/mm10.bed ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/
## .input file
# 8HM+ATAC+DNAme
grep "_14.5" /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input > ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input
grep "_0" /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input >> ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input
# CTCF
grep "CTCF" /data/zusers/fankaili/ideas/mm10_tissue_used_list.txt | grep "_0_day" | awk '{FS="\t";OFS=" "}{print $1,$2,$3}' >> ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input
grep "CTCF" /data/zusers/fankaili/ideas/mm10_tissue_used_list.txt | grep "_14.5_day" | awk '{FS="\t";OFS=" "}{print $1,$2,$3}' >> ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input
## 221 files in total
sed -i 's/_0_day /_0 /g' ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input
sed -i 's/_14.5_day /_14.5 /g' ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input

## .sh file
cp /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS/mouse_tissue_8hm_ATAC_WGBS.sh ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.sh
vim ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.sh

## .parafile file
cp /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS/mouse_tissue_8hm_ATAC_WGBS.parafile ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.parafile
vim ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.parafile

nohup bash ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.sh > ./e14.5p0_8hm_ATAC_DNAme_CTCF_1/nohup.e14.5p0_8hm_ATAC_DNAme_CTCF_1.out 2>&1&



## 2) Only use tissues with CTCF in e14.5&P0
mkdir e14.5p0_8hm_ATAC_DNAme_CTCF_2
cp /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS/mm10.bed ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/
## .input file
# 8HM+ATAC+DNAme
grep "lung_14.5" /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input > ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.input
grep "liver_14.5" /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input >> ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.input
grep "_0" /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input >> ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.input
# CTCF
grep "CTCF" /data/zusers/fankaili/ideas/mm10_tissue_used_list.txt | grep "_0_day" | awk '{FS="\t";OFS=" "}{print $1,$2,$3}' >> ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.input
grep "CTCF" /data/zusers/fankaili/ideas/mm10_tissue_used_list.txt | grep "_14.5_day" | awk '{FS="\t";OFS=" "}{print $1,$2,$3}' >> ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.input
## 121 files in total
sed -i 's/_0_day /_0 /g' ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.input
sed -i 's/_14.5_day /_14.5 /g' ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.input

## .sh file
cp /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS/mouse_tissue_8hm_ATAC_WGBS.sh ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.sh
vim ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.sh

## .parafile file
cp /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS/mouse_tissue_8hm_ATAC_WGBS.parafile ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.parafile
vim ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.parafile

nohup bash ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/e14.5p0_8hm_ATAC_DNAme_CTCF_2.sh > ./e14.5p0_8hm_ATAC_DNAme_CTCF_2/nohup.e14.5p0_8hm_ATAC_DNAme_CTCF_2.out 2>&1&


# datafile="DHS_v2_1-300bp.input"
# parafile="DHS_v2_1-300bp.parafile"
# tmpfolder="/data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v2_1-300bp_result/"
#
#
#
# datafile="DHS_v1_100-300bp.input"
# parafile="DHS_v1_100-300bp.parafile"
# tmpfolder="/data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_result/"
#
#
# datafile="run_IDEAS_8hm_atac_dname_pvalue.input"
# parafile="run_IDEAS_8hm_atac_dname_pvalue.parafile"
# tmpfolder="/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/"
