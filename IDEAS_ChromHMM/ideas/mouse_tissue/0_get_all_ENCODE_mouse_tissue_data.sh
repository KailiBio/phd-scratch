#!/bin/bash

# -- Kaili
# This is script for getting all the 870 mouse tissue file ID from ENCODE.
# 1. try to get the expID of all mouse data
# 2. get all the ENCODE mouse rep1 data, for rerunning IDEAS.
## 1) ATAC-seq
## 2) make DNAme bigWig in 1bp resolution
## 3) get 8 histone mark ChIP-seq
# 3. get rep1 signal value for all bins

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/"

cd ${workDir}

# 1. try to get the expID of all mouse data
# !. need double check if you use this file.
python ${scriptDir}get_mouse_data_from_ENCODE_json.py


# 2. get all the ENCODE mouse rep1 data, for rerunning IDEAS.
## 1) ATAC-seq
python ${scriptDir}get_ENCODE_mouse_rep1_ATAC_data.py
## 2) make DNAme bigWig in 1bp resolution
python ${scriptDir}get_ENCODE_mouse_rep1_DNAme_data.py
## 3) get 8 histone mark ChIP-seq
python ${scriptDir}get_ENCODE_mouse_rep1_8HM_data.py
## 4) get CTCF data
python ${scriptDir}get_ENCODE_mouse_rep1_CTCF_data.py
##
cat ENCODE_mouse_rep1_8HM_filelist.txt ENCODE_mouse_rep1_ATAC_filelist.txt ENCODE_mouse_rep1_CTCF_filelist.txt ENCODE_mouse_rep1_DNAme_filelist.txt \
| sort -k1,1 -k2,2 > ENCODE_mouse_rep1_filelist.txt

# 3. get rep1 signal value for all bins
normal_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/"
dhs_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
mkdir ${normal_bins_signal_path}
mkdir ${dhs_bins_signal_path}
#
encode_data_path="/data/projects/encode/data/"
#
cat ENCODE_mouse_rep1_8HM_filelist.txt ENCODE_mouse_rep1_ATAC_filelist.txt \
ENCODE_mouse_rep1_CTCF_filelist.txt > tmp_ENCODE_mouse_rep1_filelist.txt
#
if [ -f ENCODE_rep1_signal_code.sh ]
then
    rm ENCODE_rep1_signal_code.sh
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "cp "${encode_data_path}${expID}"/"${fileID}".bigWig /tmp/" >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}".tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""5}' /tmp/"${sample}"_"${assay}".tab > /tmp/"${sample}"_"${assay}".txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".tab" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".txt" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}".tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""5}' /tmp/"${sample}"_"${assay}".tab > /tmp/"${sample}"_"${assay}".txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".tab" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".txt" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "rm /tmp/"${fileID}".bigWig" >> ENCODE_rep1_signal_code.sh
done < tmp_ENCODE_mouse_rep1_filelist.txt
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "cp "${encode_data_path}${expID}"/"${fileID}".bigWig /tmp/" >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}".tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""6}' /tmp/"${sample}"_"${assay}".tab > /tmp/"${sample}"_"${assay}".txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".tab" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".txt" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}".tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""6}' /tmp/"${sample}"_"${assay}".tab > /tmp/"${sample}"_"${assay}".txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".tab" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}".txt" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "rm /tmp/"${fileID}".bigWig" >> ENCODE_rep1_signal_code.sh
done < ENCODE_mouse_rep1_DNAme_filelist.txt
##### cut code into 10 samples each
mkdir /data/zusers/fankaili/ideas/code/get_rep1_signal/
#
for i in {1..68}
do
    awk -v i="$i" '{if(NR>((i-1)*100) && NR<=(100*i)){print $0}}' ENCODE_rep1_signal_code.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_code_${i}.sh
done
