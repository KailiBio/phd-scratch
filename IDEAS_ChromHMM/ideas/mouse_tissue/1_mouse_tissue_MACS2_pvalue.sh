#!/bin/bash

# -- Kaili
# This script is for getting mouse tissue signal files.
##### Still can't access all the data.
# Some data seems not included in json (ENCSR371KFW)
# Some data (P54 WGBS) seems not download yet and also no bw file.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/"

# 1. get signal code from ENCODE json
python ${scriptDir}get_signal_file.py

cd /data/zusers/fankaili/ideas/code/get_macs2_pvalue_code/
for i in {1..10}
do
    awk -v i="$i" '{FS=OFS="\t"}{if(NR>150*(i-1) && NR<=150*i){print $0}}' get_macs2_pvalue_signal_code.sh \
    > get_macs2_pvalue_signal_code_${i}.sh
    nohup bash get_macs2_pvalue_signal_code_${i}.sh > nohup.get_macs2_pvalue_signal_code_${i}.out 2>&1&
done

# 2. run IDEAS on MACS2 p-value

## 1) mouse_tissue_8hm_ATAC_WGBS
mkdir ${workDir}mouse_tissue_8hm_ATAC_WGBS
cd ${workDir}mouse_tissue_8hm_ATAC_WGBS
cp /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed ./

# make the marks array
declare -A array
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3 ATAC-seq WGBS
do
    array[$mark]=1
done

# clear the file
if [ -f mouse_tissue_8hm_ATAC_WGBS.input ]; then rm mouse_tissue_8hm_ATAC_WGBS.input; else touch mouse_tissue_8hm_ATAC_WGBS.input; fi

# get .input file
while read line
do
    mark=`awk '{print $2}' <<< $line`
    if [[ ${array[$mark]} ]];then
        awk '{FS="\t";OFS=" "}{print $1,$2,$3}' <<< $line >> mouse_tissue_8hm_ATAC_WGBS.input
    fi
done < /data/zusers/fankaili/ideas/mm10_tissue_used_list.txt

cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.parafile mouse_tissue_8hm_ATAC_WGBS.parafile
vim mouse_tissue_8hm_ATAC_WGBS.parafile

cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.sh mouse_tissue_8hm_ATAC_WGBS.sh
vim mouse_tissue_8hm_ATAC_WGBS.sh

nohup bash mouse_tissue_8hm_ATAC_WGBS.sh > nohup.mouse_tissue_8hm_ATAC_WGBS.out 2>&1&



## 2) mouse_tissue_8hm_ATAC_WGBS_DNase
mkdir ${workDir}mouse_tissue_8hm_ATAC_WGBS_DNase
cd ${workDir}mouse_tissue_8hm_ATAC_WGBS_DNase
cp /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed ./

# make the marks array
declare -A array
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3 ATAC-seq WGBS DNase-seq
do
    array[$mark]=1
done

# clear the file
if [ -f mouse_tissue_8hm_ATAC_WGBS_DNase.input ]; then rm mouse_tissue_8hm_ATAC_WGBS_DNase.input; else touch mouse_tissue_8hm_ATAC_WGBS_DNase.input; fi

# get .input file
while read line
do
    mark=`awk '{print $2}' <<< $line`
    if [[ ${array[$mark]} ]];then
        awk '{FS="\t";OFS=" "}{print $1,$2,$3}' <<< $line >> mouse_tissue_8hm_ATAC_WGBS_DNase.input
    fi
done < /data/zusers/fankaili/ideas/mm10_tissue_used_list.txt

cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.parafile mouse_tissue_8hm_ATAC_WGBS_DNase.parafile
vim mouse_tissue_8hm_ATAC_WGBS_DNase.parafile

cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.sh mouse_tissue_8hm_ATAC_WGBS_DNase.sh
vim mouse_tissue_8hm_ATAC_WGBS_DNase.sh

nohup bash mouse_tissue_8hm_ATAC_WGBS_DNase.sh > nohup.mouse_tissue_8hm_ATAC_WGBS_DNase.out 2>&1&




## 3) mouse_tissue_8hm_ATAC_WGBS_DNase_TF
mkdir ${workDir}mouse_tissue_8hm_ATAC_WGBS_DNase_TF
cd ${workDir}mouse_tissue_8hm_ATAC_WGBS_DNase_TF
cp /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed ./

# make the marks array
declare -A array
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3 ATAC-seq WGBS DNase-seq CTCF EP300 GATA4 POLR2A
do
    array[$mark]=1
done

# clear the file
if [ -f mouse_tissue_8hm_ATAC_WGBS_DNase_TF.input ]; then rm mouse_tissue_8hm_ATAC_WGBS_DNase_TF.input; else touch mouse_tissue_8hm_ATAC_WGBS_DNase_TF.input; fi

# get .input file
while read line
do
    mark=`awk '{print $2}' <<< $line`
    if [[ ${array[$mark]} ]];then
        awk '{FS="\t";OFS=" "}{print $1,$2,$3}' <<< $line >> mouse_tissue_8hm_ATAC_WGBS_DNase_TF.input
    fi
done < /data/zusers/fankaili/ideas/mm10_tissue_used_list.txt

cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.parafile mouse_tissue_8hm_ATAC_WGBS_DNase_TF.parafile
vim mouse_tissue_8hm_ATAC_WGBS_DNase_TF.parafile

cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.sh mouse_tissue_8hm_ATAC_WGBS_DNase_TF.sh
vim mouse_tissue_8hm_ATAC_WGBS_DNase_TF.sh

nohup bash mouse_tissue_8hm_ATAC_WGBS_DNase_TF.sh > nohup.mouse_tissue_8hm_ATAC_WGBS_DNase_TF.out 2>&1&



#####################
# Oct 02
cd /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS/
mv mouse_tissue_8hm_ATAC_WGBS.input mouse_tissue_8hm_ATAC_WGBS.input0
sort -u mouse_tissue_8hm_ATAC_WGBS.input0 > mouse_tissue_8hm_ATAC_WGBS.input
vim mouse_tissue_8hm_ATAC_WGBS.sh
# change result folder
nohup bash mouse_tissue_8hm_ATAC_WGBS.sh > nohup.mouse_tissue_8hm_ATAC_WGBS_2.out 2>&1&



cd /data/zusers/fankaili/ideas/mouse_tissue_8hm_ATAC_WGBS_DNase/
mv mouse_tissue_8hm_ATAC_WGBS_DNase.input mouse_tissue_8hm_ATAC_WGBS_DNase.input0
cp mouse_tissue_8hm_ATAC_WGBS_DNase.input0 ss.txt
vim ss.txt
sed -i 's/.50_day/.5_day/g' ss.txt
sort -u ss.txt > mouse_tissue_8hm_ATAC_WGBS_DNase.input
vim mouse_tissue_8hm_ATAC_WGBS_DNase.sh
# change result folder
nohup bash mouse_tissue_8hm_ATAC_WGBS_DNase.sh > nohup.mouse_tissue_8hm_ATAC_WGBS_DNase_2.out 2>&1&



#####################
# heart only
mkdir /data/zusers/fankaili/ideas/heart/
cd /data/zusers/fankaili/ideas/heart/
