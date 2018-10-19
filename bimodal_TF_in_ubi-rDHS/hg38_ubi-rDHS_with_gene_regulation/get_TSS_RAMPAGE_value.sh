#!/bin/bash

# -- Kaili
# This script is for getting TSS RAMPAGE signal after bigWigAverageOverBed result.
# INPUT:
# OUTPUT:
# EXP: bash get_TSS_RAMPAGE_value.sh /data/zusers/fankaili/ccre/hg38_ubi-rDHS/ /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/ \
#      /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/ TSS.Filtered.uniq.bed hg38_RAMPAGE_list.txt

work_dir=$1
rampage_file_path=$2
out_fir=$3
tss_bed_file=$4
file_list=$5

# work_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"
# rampage_file_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/"
# out_fir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/"
# tss_bed_file="TSS.Filtered.uniq.bed"
# file_list="hg38_RAMPAGE_list.txt"


cd ${work_dir}

# 1. divide tss_bed_file based on strand
file_name=${tss_bed_file%.bed}
awk '{FS=OFS="\t"}{if($6=="+"){print $0}}' ${tss_bed_file} > ${filename}_plus.bed
awk '{FS=OFS="\t"}{if($6=="-"){print $0}}' ${tss_bed_file} > ${filename}_minus.bed

# 2. get RAMPAGE signal

## function for getting RAMPAGE signal
getRampageSignal(){
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$4}else{print $8,a[$4]}}' ${rampage_file_path}$1.tab ${work_dir}${filename}_plus.bed > \
    ${out_fir}$3_rampage.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$4}else{print $8,a[$4]}}' ${rampage_file_path}$2.tab ${work_dir}${filename}_minus.bed > \
    ${out_fir}$3_rampage.txt
}


while read line
do
    plus=`awk '{FS=OFS="\t"}{print $2}' <<< $line`
    minus=`awk '{FS=OFS="\t"}{print $3}' <<< $line`
    biosample=`awk '{FS=OFS="\t"}{print $4}' <<< $line`
    echo ${plus}
    #getRampageSignal ${plus} ${minus} ${biosample}
done < ${file_list}
