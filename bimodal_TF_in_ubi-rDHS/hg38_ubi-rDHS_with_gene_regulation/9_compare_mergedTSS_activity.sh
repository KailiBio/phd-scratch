#!/bin/bash

# -- Kaili
# This script is for comparing RAMPAGE signal between ubi-rOCR overlapped and cell type-active OCR overlapped merged-TSSs.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/"

cd ${workDir}

# 1. get non ubi-rOCRs overlapped merged-TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' GRCh38_ubi-rOCR_overlapped_merged-TSS.bed \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed > GRCh38_non_ubi-rOCR_overlapped_mergedTSS.bed

# 2. gene expression comparison
cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/
#
mkdir mergedTSS_exp_comparison_file
mkdir mergedTSS_exp_comparison_pdf
#
while read line
do
    echo ${line};
    rampage_exp_id=`awk '{FS=OFS="\t"}{print $7}' <<< ${line}` ;
    sample=`awk '{FS=OFS="\t"}{print $9}' <<< ${line}` ;
    dnase_exp_id=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    dnase_file_id=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    #
    bedtools intersect -a ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed -b ./merged-TSS/GRCh38_non_ubi-rOCR_overlapped_mergedTSS.bed -wa -wb \
    | cut -f 8 | sort -u > temp.txt ;
    # ubi-rOCR overlapped TSS
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]==1){print $1,$2,"ubi-rOCR_overlapped"}}}' \
    ./merged-TSS/GRCh38_ubi-rOCR_overlapped_merged-TSS.bed ./rampage_tissue_mergedTSS/${sample}_rampage.txt > ./mergedTSS_exp_comparison_file/${sample}.txt ;
    # rOCR overlapped TSS
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"active-rOCR_overlapped"}}}' \
    temp.txt ./rampage_tissue_mergedTSS/${sample}_rampage.txt >> ./mergedTSS_exp_comparison_file/${sample}.txt ;
    #
    Rscript ${scriptDir}make_comparison_barplot.R ${sample} "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/mergedTSS_exp_comparison_file/" \
    "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/mergedTSS_exp_comparison_pdf/";
done < hg38_RAMPAGE_DNase_matched_list.txt
rm temp.txt
