#!/bin/bash

# -- Kaili
# This script is for comparing with ChromInpute & Avocado -- run whole genome.

# 0. get whole genome 25bp bins signal
# 1. ChromImpute
#2. Avocado

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/imputation_comparison/"

cd ${workDir}

# 0. get whole genome 25bp bins signal
## 1) get bed file
sort -k1,1 -k2,2n /data/zusers/fankaili/ideas/dhs_ctcf/mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed > tmp.bed
bedtools merge -i tmp.bed > mm10_dhs_bins_merged.bed
bedtools makewindows -b mm10_dhs_bins_merged.bed -w 25 | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $0,"bin_"NR}' > mm10_25bp.bed
rm tmp.bed
## 2) calculate signal for 25bp bins
### get chrom.size files. (original one and modified one)
cp /home/fankaili/genome/mm10.chrom.sizes.clean ./chromimpute_wholeGenome/
awk '{FS=OFS="\t"}{a=int($2/25);b=$2%25;if(b!=0){a+=1};print $1,a*25}' ./chromimpute_wholeGenome/mm10.chrom.sizes.clean > ./chromimpute_wholeGenome/mm10.chrom.sizes.clean2
#
encode_data_path="/data/projects/encode/data/"
signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_25bp_bins/"
bedfile="/data/zusers/fankaili/ideas/imputation_comparison/mm10_25bp.bed"
mkdir ${signal_path}
#
if [ -f ENCODE_rep1_signal_all_25bp_code.sh ]; then rm ENCODE_rep1_signal_all_25bp_code.sh; fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed "${encode_data_path}${expID}"/"${fileID}".bigWig "${bedfile}" "${signal_path}${sample}"_"${assay}"_25bp.tab" >> ENCODE_rep1_signal_all_25bp_code.sh
    echo "awk '{print ""\$""5}' "${signal_path}${sample}"_"${assay}"_25bp.tab > "${signal_path}${sample}"_"${assay}"_25bp.txt" >> ENCODE_rep1_signal_all_25bp_code.sh
done < ENCODE_mouse_rep1_CTCFsample_filelist1.txt
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig "${bedfile}" "${signal_path}${sample}"_"${assay}"_25bp.tab" >> ENCODE_rep1_signal_all_25bp_code.sh
    echo "awk '{print ""\$""6}' "${signal_path}${sample}"_"${assay}"_25bp.tab > "${signal_path}${sample}"_"${assay}"_25bp.txt" >> ENCODE_rep1_signal_all_25bp_code.sh
done < ENCODE_mouse_rep1_CTCFsample_filelist2.txt
###
for i in {1..11}
do
    awk -v i="$i" '{if(NR>((i-1)*22) && NR<=(22*i)){print $0}}' ENCODE_rep1_signal_all_25bp_code.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_all_25bp_code_${i}.sh
    nohup bash /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_all_25bp_code_${i}.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/nohup.ENCODE_rep1_signal_all_25bp_code_${i}.out 2>&1&
done
# z001
## 3) get bed in signal order
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$0}else{print a[$1]}}' mm10_25bp.bed /data/zusers/fankaili/ideas/signal/rep1_signal_25bp_bins/forebrain_0_ATAC_25bp.tab > mm10_25bp_inOrder.bed
awk '{OFS=" "}{print $1,$2,$3,$4}' mm10_25bp_inOrder.bed > mm10_25bp_inOrder_space.bed

mkdir chromimpute_wholeGenome
mkdir avocado_wholeGenome

# 1. ChromImpute
## 1) convert singal into wig
nohup bash ${scriptDir}Convert_signal_to_zippedWig_ChromImpute.sh > ./nohup/nohup.Convert_signal_to_zippedWig_ChromImpute.out 2>&1&
# z001 42395
## 2) get master table
if [ -f ./chromimpute_wholeGenome/master_table.txt ];then rm ./chromimpute_wholeGenome/master_table.txt; fi
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF DNAme
    do
        echo -e ${sample}"\t"${mark}"\t"${sample}"_"${mark}"_25bp" >> ./chromimpute_wholeGenome/master_table.txt
    done
done
vim ./chromimpute_wholeGenome/master_table.txt
## 3) run ChromImpute
nohup bash ${scriptDir}run_ChromImpute_pipeline.sh /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_wholeGenome/ /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_wholeGenome/mm10.chrom.sizes.clean > ./nohup/nohup.run_ChromImpute_pipeline.out 2>&1&
# z001

#2. Avocado
## 1) convert signal file into compressed numpy
nohup python ${scriptDir}convert_signal_to_zippedNumpy_Avocado.py /data/zusers/fankaili/ideas/signal/rep1_signal_25bp_bins/ > ./nohup/nohup.convert_signal_to_zippedNumpy_Avocado.out 2>&1&
# z001 50798
## 2) run avocado using default
nohup python ${scriptDir}run_Avocado_pipeline.py /data/zusers/fankaili/ideas/signal/rep1_signal_25bp_bins/ /data/zusers/fankaili/ideas/imputation_comparison/avocado_wholeGenome/ > ./nohup/nohup.run_Avocado_pipeline.out 2>&1&
# z001
