#!/bin/bash

# -- Kaili
# This script is for calculating RAMPAGE signal in merged-TSSs for all files, then get RAMPAGE signal for merged-TSSs based on strand.
# 1. calculate rampage signal in all files for merged-TSS
# 2. get rampage signal for merged-TSS by strand

rampage_signal_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_mergedTSS/"
rampage_signal_tissue_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/"
tss_plus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/GRCh38_merged-TSS_gene_length_plus.txt"
tss_minus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/GRCh38_merged-TSS_gene_length_minus.txt"
mergedTSS_bed="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/GRCh38_merged-TSS_gene_50bp_bed6.bed"

## function for getting RAMPAGE signal
getRampageSignal(){
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$4}else{print $4,a[$4]}}' ${rampage_signal_dir}$1.tab  ${tss_plus} > ${rampage_signal_tissue_dir}$3_rampage.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$4}else{print $4,a[$4]}}' ${rampage_signal_dir}$2.tab  ${tss_minus} >> ${rampage_signal_tissue_dir}$3_rampage.txt
}
#

## calculate signal
while read line
do
    echo ${line}
    id=`awk '{print $1}' <<< ${line}`
    plus=`awk '{print $2}' <<< ${line}`
    minus=`awk '{print $3}' <<< ${line}`
    biosample=`awk '{print $4}' <<< ${line}`
    # 1. calculate rampage signal in all files for merged-TSS
    bigWigAverageOverBed /data/projects/encode/data/${id}/${plus}.bigWig ${mergedTSS_bed} ${rampage_signal_dir}${plus}.tab
    bigWigAverageOverBed /data/projects/encode/data/${id}/${minus}.bigWig ${mergedTSS_bed} ${rampage_signal_dir}${minus}.tab
    # 2. get rampage signal for merged-TSS by strand
    getRampageSignal ${plus} ${minus} ${biosample}
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt
