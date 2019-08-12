#!/bin/bash

# -- Kaili
# This script is for enrichment analysis for new promoter peak shape result.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/"
promoterShapeResult="/data/zusers/zhangx/projects/rampage/0_rampage_peak/peak_uniq/"

cd ${workDir}

#####

non_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed"
ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
rampageDir="/data/zusers.ds/zhangx/projects/rampage/0_rampage_peak/kaili/20190624/filtered/motif/"
#
if [ -f peak_shape_count_155.txt ];then rm peak_shape_count_155.txt; fi
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    sample=`awk '{print $4}' <<< ${line}`
    echo ${expID}
    #
    non_ubi_np=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_np.txt -u | wc -l`
    non_ubi_bp=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_bp.txt -wa -u | wc -l`
    non_ubi_wp=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_wp.txt -wa -u | wc -l`
    #
    ubi_np=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_np.txt -u | wc -l`
    ubi_bp=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_bp.txt -wa -u | wc -l`
    ubi_wp=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_wp.txt -wa -u | wc -l`
    #
    echo -e ${expID}"\t"${sample}"\t"${ubi_np}"\t"${ubi_bp}"\t"${ubi_wp}"\t"${non_ubi_np}"\t"${non_ubi_bp}"\t"${non_ubi_wp} >> peak_shape_count_155.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt



#########
# TSS filter
non_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/non-ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed"
ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed"
rampageDir="/data/zusers.ds/zhangx/projects/rampage/0_rampage_peak/kaili/20190624/filtered/motif/"
#
if [ -f peak_shape_count_155_TSSfiltered.txt ];then rm peak_shape_count_155_TSSfiltered.txt; fi
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    sample=`awk '{print $4}' <<< ${line}`
    echo ${expID}
    #
    non_ubi_np=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_np.txt -u | wc -l`
    non_ubi_bp=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_bp.txt -wa -u | wc -l`
    non_ubi_wp=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_wp.txt -wa -u | wc -l`
    #
    ubi_np=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_np.txt -u | wc -l`
    ubi_bp=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_bp.txt -wa -u | wc -l`
    ubi_wp=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_wp.txt -wa -u | wc -l`
    #
    num_np=`wc -l ${rampageDir}${expID}_rampage_peaks_np.txt | awk '{print $1}'`
    num_bp=`wc -l ${rampageDir}${expID}_rampage_peaks_bp.txt | awk '{print $1}'`
    num_wp=`wc -l ${rampageDir}${expID}_rampage_peaks_wp.txt | awk '{print $1}'`
    #
    echo -e ${expID}"\t"${sample}"\t"${ubi_np}"\t"${ubi_bp}"\t"${ubi_wp}"\t"${non_ubi_np}"\t"${non_ubi_bp}"\t"${non_ubi_wp}"\t"${num_np}"\t"${num_bp}"\t"${num_wp} >> peak_shape_count_155_TSSfiltered.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt


intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed -b /data/zusers.ds/zhangx/projects/rampage/0_rampage_peak/kaili/20190624/filtered/motif/ENCSR497BYB_rampage_peaks_np.txt -wa -wb | head



intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -b /data/zusers.ds/zhangx/projects/rampage/0_rampage_peak/kaili/20190624/filtered/motif/ENCSR201ARN_rampage_peaks_np.txt | head

############

awk '{FS=OFS="\t"}{if($2>50){print $1,$2-50,$3+50,$4,$5,$6,$7}else{print $1,0,$3+50,$4,$5,$6,$7}}' /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed > hg38_v28_basic_TSS_filtered_50bp.bed
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -b hg38_v28_basic_TSS_filtered_50bp.bed -u | sort -k1,1 -k2,2n > non-ubi-rOCRs_overlap_hg38_v28_basic_TSS_50bp.bed
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -b hg38_v28_basic_TSS_filtered_50bp.bed -u | sort -k1,1 -k2,2n > ubi-rOCRs_overlap_hg38_v28_basic_TSS_50bp.bed


non_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/non-ubi-rOCRs_overlap_hg38_v28_basic_TSS_50bp.bed"
ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS_50bp.bed"
rampageDir="/data/zusers.ds/zhangx/projects/rampage/0_rampage_peak/kaili/20190624/filtered/motif/"
#
if [ -f peak_shape_count_155_TSSfiltered_50bp.txt ];then rm peak_shape_count_155_TSSfiltered_50bp.txt; fi
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    sample=`awk '{print $4}' <<< ${line}`
    echo ${expID}
    #
    non_ubi_np=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_np.txt -u | wc -l`
    non_ubi_bp=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_bp.txt -wa -u | wc -l`
    non_ubi_wp=`intersectBed -a ${non_ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_wp.txt -wa -u | wc -l`
    #
    ubi_np=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_np.txt -u | wc -l`
    ubi_bp=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_bp.txt -wa -u | wc -l`
    ubi_wp=`intersectBed -a ${ubi_rOCRs} -b ${rampageDir}${expID}_rampage_peaks_wp.txt -wa -u | wc -l`
    #
    echo -e ${expID}"\t"${sample}"\t"${ubi_np}"\t"${ubi_bp}"\t"${ubi_wp}"\t"${non_ubi_np}"\t"${non_ubi_bp}"\t"${non_ubi_wp} >> peak_shape_count_155_TSSfiltered_50bp.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt
