#!/bin/bash

# -- Kaili
# This script is for getting RNA-seq signal for gene whose TSSs overlap ubi-rOCRs and non-ubi active-rOCRs.

# INPUT:
# OUTPUT: sample_RNA_signal_comparison.txt
# EXP: bash get_RNAseq_signal_two_groups.sh
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

gene_labeled_file=$1
tss_labeled_file=$2
exp_signal_path=$3
non_ubi_outPath=$4
outPath=$5

matched_file_list="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_matched_DNase_RNA_RAMPAGE_list.txt"

######################

get_RNAseq_signal(){
    matched_file_list=$1
    gene_labeled_file=$2
    tss_labeled_file=$3
    exp_signal_path=$4
    non_ubi_outPath=$5
    outPath=$6
    #
    if [ -f ${outPath}sample_RNA_signal_comparison.txt ]; then rm ${outPath}sample_RNA_signal_comparison.txt; fi
    #
    while read line
    do
        dnase_exp_id=`awk '{print $1}' <<< ${line}`
        rna_exp_id=`awk '{print $3}' <<< ${line}`
        sample=`awk '{print $7}' <<< ${line}`
        # get signal of gene that overlapping ubi-rOCRs
        num_ubi=`grep "genes_whose_TSSs_overlap_ubi-rOCRs" ${gene_labeled_file} | wc -l`
        awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){if($2=="genes_whose_TSSs_overlap_ubi-rOCRs"){a[$1]=1}}else{if(a[$1]){print $1,$2,"genes_whose_TSSs_overlap_ubi-rOCRs",sample,num}}}' ${gene_labeled_file} ${exp_signal_path}${rna_exp_id}.txt >> ${outPath}sample_RNA_signal_comparison.txt
        # get signal of gene that overlapping not-ubi cell-type active-rOCRs
        intersectBed -a ${tss_labeled_file} -b ${non_ubi_outPath}${dnase_exp_id}_OCR.bed -wa -u | cut -f 7 | sort -u > ${outPath}tmp.gene.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){if($2=="genes_whose_TSSs_overlap_ubi-rOCRs"){a[$1]=1}}else{if(a[$1]!=1){print $0}}}' ${gene_labeled_file} ${outPath}tmp.gene.txt > ${outPath}tmp.gene_2.txt
        num=`wc -l ${outPath}tmp.gene_2.txt | awk '{print $1}'`
        awk -v sample="$sample" -v num="$num" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2,"genes_whose_TSSs_overlap_other_active_rOCRs",sample,num}}}' ${outPath}tmp.gene_2.txt ${exp_signal_path}${rna_exp_id}.txt >> ${outPath}sample_RNA_signal_comparison.txt
    done < ${matched_file_list}
    rm ${outPath}tmp.gene.txt ${outPath}tmp.gene_2.txt
}


######################
get_RNAseq_signal ${matched_file_list} ${gene_labeled_file} ${tss_labeled_file} ${exp_signal_path} ${non_ubi_outPath} ${outPath}
