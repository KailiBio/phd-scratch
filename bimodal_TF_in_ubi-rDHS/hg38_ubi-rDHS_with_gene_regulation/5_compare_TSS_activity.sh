#!/bin/bash

# -- Kaili
# This script is for comparing TSS activity. (ubi-rOCRs overlapped with cell-type active OCRs overlapped TSSs)
# 1. get RAMPAGE matched DNase data. (the same donor ID)
# 2. gene expression comparison

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. get RAMPAGE matched DNase data. (the same donor ID)

## 1) get RNA-RAMPAGE donor ID list
python ${scriptDir}get_donorID.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list_duplicate.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list_donor0.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$4}else{print $1,$2,a[$1],$4,$5,$6}}' hg38_RAMPAGE_list.txt hg38_RAMPAGE_list_donor0.txt \
> hg38_RAMPAGE_list_donor.txt

## 2) matching
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=1;b[$6]=1;c[$5]=$0;d[$6]=$0}else{if(a[$5] && ($5!="---")){print $0,c[$5]}else if(a[$6] && ($6!="---")){print $0,c[$6]}else if(b[$5] && ($5!="---")){print $0,d[$5]}else if(b[$6] && ($6!="---")){print $0,d[$6]}}}' \
hg38_RAMPAGE_list_donor.txt hg38_DNase_file_list_donor.txt > hg38_RAMPAGE_DNase_matched_list.txt
# get 16 matched biosamples


# 2. gene expression comparison
#
cut -f 2 ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID.txt > ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID_list.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$8,$5,$6,$7}}}' hg38_ubi-rOCR_non_overlapped_TSS_list.txt TSS.Filtered.uniqID.bed | \
sort -u > hg38_ubi-rOCR_non_overlapped_TSS_uniqID.bed
#
mkdir TSS_exp_comparison_file
mkdir TSS_exp_comparison_pdf
#
while read line
do
    echo ${line};
    rampage_exp_id=`awk '{FS=OFS="\t"}{print $7}' <<< ${line}` ;
    sample=`awk '{FS=OFS="\t"}{print $9}' <<< ${line}` ;
    dnase_exp_id=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    dnase_file_id=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if(a[$4]){print $0}}}' \
    /data/projects/psychencode/Registry/V1/GRCh38/Signal-Files/${dnase_exp_id}"-"${dnase_file_id}.txt hg38_non_ubi-rOCRs.bed \
    > ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed ;
    #
    bedtools intersect -a ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed -b hg38_ubi-rOCR_non_overlapped_TSS_uniqID.bed -wa -wb \
    | cut -f 8 | sort -u > temp.txt ;
    # ubi-rOCR overlapped TSS
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"ubi-rOCR_overlapped"}}}' \
    ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID_list.txt ./rampage_tissue/${sample}_rampage.txt > ./TSS_exp_comparison_file/${sample}.txt ;
    # rOCR overlapped TSS
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"active-rOCR_overlapped"}}}' \
    temp.txt ./rampage_tissue/${sample}_rampage.txt >> ./TSS_exp_comparison_file/${sample}.txt ;
    #
    Rscript ${scriptDir}make_comparison_barplot.R ${sample} "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS_exp_comparison_file/" \
    "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS_exp_comparison_pdf/";
done < hg38_RAMPAGE_DNase_matched_list.txt
rm temp.txt
