#!/bin/bash

# -- Kaili
# This script is for getting DNase signal for ubi-rOCRs and non-ubi active-rOCRs.

# INPUT: path for DNase signal of all rOCRs
#              path for storing not-ubi active-rOCRs in each biosample. (can be empty, will make those files in this script.)
#              path for output file
# OUTPUT: sample_DNase_signal_comparison.txt (matrix with all DNase signal with label)
# EXP: bash get_DNase_signal_two_groups.sh /data/projects/psychencode/Registry/V1/GRCh38/Signal-Files/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/

DNase_signal_path=$1
non_ubi_outPath=$2
outPath=$3

hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
hg38_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed"
matched_file_list="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_matched_DNase_RNA_RAMPAGE_list.txt"

######################

## this function for getting non-ubi cell-type active rOCRs.
get_non_ubi_active_rOCRs(){
    matched_file_list=$1
    DNase_signal_path=$2
    hg38_rOCRs=$3
    hg38_ubi_rOCRs=$4
    non_ubi_outPath=$5
    # get non-ubi-rOCRs
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' ${hg38_ubi_rOCRs} ${hg38_rOCRs} \
    > tmp.non_ubi_rOCRs.bed
    #
    while read line
    do
        dnase_exp_id=`awk '{print $1}' <<< ${line}`
        dnase_file_id=`awk '{print $2}' <<< ${line}`
        sample=`awk '{print $7}' <<< ${line}`
        #
        awk '{FS=OFS="\t"}{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if(a[$4]){print $0}}}' \
        ${DNase_signal_path}${dnase_exp_id}-${dnase_file_id}.txt tmp.non_ubi_rOCRs.bed > ${non_ubi_outPath}${dnase_exp_id}_OCR.bed
    done < ${matched_file_list}
}

## This function for getting DNase signal of ubi-rOCRs and non-ubi active rOCRs.
get_DNase_signal(){
    matched_file_list=$1
    DNase_signal_path=$2
    hg38_ubi_rOCRs=$3
    non_ubi_outPath=$4
    outPath=$5
    #
    if [ -f ${outPath}sample_DNase_signal_comparison.txt ]; then rm ${outPath}sample_DNase_signal_comparison.txt; fi
    #
    while read line
    do
        dnase_exp_id=`awk '{print $1}' <<< ${line}`
        dnase_file_id=`awk '{print $2}' <<< ${line}`
        sample=`awk '{print $7}' <<< ${line}`
        # get ubi-rOCRs signal
        num_ubi=`wc -l ${hg38_ubi_rOCRs} | awk '{print $1}'`
        awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$2,"ubi-rOCR",sample,num}}}' \
        ${hg38_ubi_rOCRs} ${DNase_signal_path}${dnase_exp_id}-${dnase_file_id}.txt >> ${outPath}sample_DNase_signal_comparison.txt
        # get non-ubi active-rOCRs signal
        num=`wc -l ${non_ubi_outPath}${dnase_exp_id}_OCR.bed | awk '{print $1}'`
        awk -v sample="$sample" -v num="$num" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$2,"non-ubi_active-rOCR",sample,num}}}' \
        ${non_ubi_outPath}${dnase_exp_id}_OCR.bed ${DNase_signal_path}${dnase_exp_id}-${dnase_file_id}.txt >> \
        ${outPath}sample_DNase_signal_comparison.txt
    done < ${matched_file_list}
}

######################

get_non_ubi_active_rOCRs ${matched_file_list} ${DNase_signal_path} ${rOCRs_file} ${hg38_ubi_rOCRs} ${non_ubi_outPath}
get_DNase_signal ${matched_file_list} ${DNase_signal_path} ${hg38_ubi_rOCRs} ${non_ubi_outPath} ${outPath}
