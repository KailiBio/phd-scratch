#!/bin/bash

# -- Kaili
# This script is for getting RAMPAGE signal for TSSs that overlap ubi-rOCRs and non-ubi cell-type active-rOCRs.

# INPUT:
# OUTPUT: sample_RAMPAGE_signal_comparison.txt
# EXP: bash get_RAMPAGE_signal_two_groups.sh
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

tss_labeled_file=$1
tss_signal_path=$2
non_ubi_outPath=$3
outPath=$4

matched_file_list="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_matched_DNase_RNA_RAMPAGE_list.txt"
hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
######################

get_RAMPAGE_signa(){
    matched_file_list=$1
    tss_labeled_file=$2
    tss_signal_path=$3
    non_ubi_outPath=$4
    hg38_ubi_rOCRs=$5
    outPath=$6
    #
    if [ -f ${outPath}sample_RAMPAGE_signal_comparison.txt ]; then rm ${outPath}sample_RAMPAGE_signal_comparison.txt; fi
    # get non-ubi TSSs
    intersectBed -a ${tss_labeled_file} -b ${hg38_ubi_rOCRs} -v -wa | sort -u > ${outPath}tmp.bed
    #
    while read line
    do
        dnase_exp_id=`awk '{print $1}' <<< ${line}`
        sample=`awk '{print $7}' <<< ${line}`
        #
        signal_file=`ls ${tss_signal_path} | grep  ${sample}`
        # get signal of TSSs that overlap ubi-rOCRs
        num_ubi=`grep "overlap_with_ubi-rOCRs" ${tss_labeled_file} | wc -l`
        awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){if($8=="overlap_with_ubi-rOCRs"){split($4,b,"_");a[b[1]b[2]]=1}}else{if(a[$1]){print $1,$2,"TSS_overlap_ubi-rOCRs",sample, num}}}' ${tss_labeled_file} ${tss_signal_path}${signal_file} >> ${outPath}sample_RAMPAGE_signal_comparison.txt
        # get signal of TSSs that overlap non-ubi cell-type active-rOCRs
        intersectBed -a ${outPath}tmp.bed -b ${non_ubi_outPath}${dnase_exp_id}_OCR.bed -wa -u > ${outPath}tmp.txt
        num=`wc -l ${outPath}tmp.txt | awk '{print $1}'`
        awk -v sample="$sample" -v num="$num" '{FS=OFS="\t"}{if(NR==FNR){split($4,b,"_");a[b[1]b[2]]=1}else{if(a[$1]){print $1,$2,"TSS_overlap_non-ubi_active-rOCRs",sample, num}}}' ${outPath}tmp.txt ${tss_signal_path}${signal_file} >> ${outPath}sample_RAMPAGE_signal_comparison.txt
    done < ${matched_file_list}
    rm ${outPath}tmp.bed ${outPath}tmp.txt
}


######################
get_RAMPAGE_signa ${matched_file_list} ${tss_labeled_file} ${tss_signal_path} ${non_ubi_outPath} ${hg38_ubi_rOCRs} ${outPath}
