#!/bin/bash

# -- Kaili
# This script is for getting RAMPAGE signal for TSSs that overlap ubi-rOCRs and non-ubi cell-type active-rOCRs.

# INPUT:
# OUTPUT: sample_RAMPAGE_signal_comparison.txt
# EXP: bash get_RAMPAGE_signal_two_groups.sh
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/hg38_v28_comprehensive_TSS_filtered_with_uniqID.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

tss_labeled_file=$1
tss_with_id_file=$2
tss_signal_path=$3
non_ubi_outPath=$4
outPath=$5


matched_file_list="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_matched_DNase_RNA_RAMPAGE_list.txt"
hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
######################

get_RAMPAGE_signa(){
    matched_file_list=$1
    tss_labeled_file=$2
    tss_with_id_file=$3
    tss_signal_path=$4
    non_ubi_outPath=$5
    hg38_ubi_rOCRs=$6
    outPath=$7
    #
    if [ -f ${outPath}sample_RAMPAGE_signal_comparison.txt ]; then rm ${outPath}sample_RAMPAGE_signal_comparison.txt; fi
    # get non-ubi TSSs
    intersectBed -a ${tss_labeled_file} -b ${hg38_ubi_rOCRs} -v -wa | sort -u > ${outPath}tmp.bed
    #
    while read line
    do
        dnase_exp_id=`awk '{print $1}' <<< ${line}`
        rampage_exp_id=`awk '{print $5}' <<< ${line}`
        sample=`awk '{print $7}' <<< ${line}`
        #
        # signal_file=`ls ${tss_signal_path} | grep  ${sample}`
        # get signal of TSSs that overlap ubi-rOCRs
        num_ubi=`grep "overlap_with_ubi-rOCRs" ${tss_labeled_file} | cut -f 4 | sort -u | wc -l`
        # awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){if($8=="overlap_with_ubi-rOCRs"){split($4,b,"_");a[b[1]b[2]]=1}}else{if(a[$1]){print $1,$2,"TSS_overlap_ubi-rOCRs",sample, num}}}' ${tss_labeled_file} ${tss_signal_path}${rampage_exp_id}.txt >> ${outPath}sample_RAMPAGE_signal_comparison.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){if($8=="overlap_with_ubi-rOCRs"){a[$4]=1}}else{if(a[$8]){print $4,$8}}}' ${tss_labeled_file} ${tss_with_id_file} > ${outPath}tmp.ubi.txt
        awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$1]){print a[$1],$2,"TSS_overlap_ubi-rOCRs",sample, num}}}' ${outPath}tmp.ubi.txt ${tss_signal_path}${rampage_exp_id}.txt | sort -u >> ${outPath}sample_RAMPAGE_signal_comparison.txt
        # get signal of TSSs that overlap non-ubi cell-type active-rOCRs
        intersectBed -a ${outPath}tmp.bed -b ${non_ubi_outPath}${dnase_exp_id}_OCR.bed -wa -u > ${outPath}tmp.txt
        num=`wc -l ${outPath}tmp.txt | cut -f 4 | sort -u | awk '{print $1}'`
        # awk -v sample="$sample" -v num="$num" '{FS=OFS="\t"}{if(NR==FNR){split($4,b,"_");a[b[1]b[2]]=1}else{if(a[$1]){print $1,$2,"TSS_overlap_non-ubi_active-rOCRs",sample, num}}}' ${outPath}tmp.txt ${tss_signal_path}${rampage_exp_id}.txt >> ${outPath}sample_RAMPAGE_signal_comparison.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$8]){print $4,$8}}}' ${outPath}tmp.txt ${tss_with_id_file} > ${outPath}tmp.active.txt
        awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$1]){print a[$1],$2,"TSS_overlap_non-ubi_active-rOCRs",sample, num}}}' ${outPath}tmp.active.txt ${tss_signal_path}${rampage_exp_id}.txt | sort -u >> ${outPath}sample_RAMPAGE_signal_comparison.txt
    done < ${matched_file_list}
    rm ${outPath}tmp.bed ${outPath}tmp.txt ${outPath}tmp.ubi.txt ${outPath}tmp.active.txt
}

######################
get_RAMPAGE_signa ${matched_file_list} ${tss_labeled_file} ${tss_with_id_file} ${tss_signal_path} ${non_ubi_outPath} ${hg38_ubi_rOCRs} ${outPath}
