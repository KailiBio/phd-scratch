#!/bin/bash

# -- Kaili
# This script is for making DNase&RNA&RAMPAGE signal comparison.
# fig1

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. get match 16 files
cut -f 1-2,7-9 hg38_RNA_DNase_matched_list.txt > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$7;b[$1]=$8}else{if(a[$1]){print $1,$2,$3,$4,a[$1],b[$1],$5}}}' \
hg38_RAMPAGE_DNase_matched_list.txt tmp.txt > hg38_matched_DNase_RNA_RAMPAGE_list.txt


# 2. DNase-seq: ubi-rOCRs vs. non-ubi active-rOCRs
if [ -f GRCh38_OCRs_DNase_signal_comparison.txt ]
then
    rm GRCh38_OCRs_DNase_signal_comparison.txt
fi
#
signal_path="/data/projects/psychencode/Registry/V1/GRCh38/Signal-Files/"
while read line
do
    dnase_exp_id=`awk '{print $1}' <<< ${line}`
    dnase_file_id=`awk '{print $2}' <<< ${line}`
    sample=`awk '{print $7}' <<< ${line}`
    #
    signal_file=`ls /data/projects/psychencode/Registry/V1/GRCh38/Signal-Files/ | grep ${dnase_file_id}`
    # get ubi-rOCRs signal
    awk -v sample="$sample" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$2,"ubi-rOCR",sample,11888}}}' \
    GRCh38_ubi-rOCRs_EDGEid.bed ${signal_path}${dnase_exp_id}-${dnase_file_id}.txt >> GRCh38_OCRs_DNase_signal_comparison.txt
    # get non-ubi active-rOCRs signal
    num=`wc -l ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed | awk '{print $1}'`
    awk -v sample="$sample" -v num="$num" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$2,"non-ubi-active-rOCR",sample,num}}}' \
    ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed ${signal_path}${signal_file} >> \
    GRCh38_OCRs_DNase_signal_comparison.txt
done < hg38_matched_DNase_RNA_RAMPAGE_list.txt


# 3. RNA-seq: ubi-rOCRs vs. non-ubi active-rOCRs overlapped genes
if [ -f GRCh38_OCRs_RNA_signal_comparison.txt ]
then
    rm GRCh38_OCRs_RNA_signal_comparison.txt
fi
#
while read line
do
    dnase_exp_id=`awk '{print $1}' <<< ${line}`
    dnase_file_id=`awk '{print $2}' <<< ${line}`
    rna_exp_id=`awk '{print $3}' <<< ${line}`
    rna_file_id=`awk '{print $4}' <<< ${line}`
    sample=`awk '{print $7}' <<< ${line}`
    # get non-ubi active-rOCRs overlapped genes
    intersectBed -a ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed -b \
    hg38_TSS_in_non_ubi-rOCR_overlapped_gene.bed -wa -wb | cut -f 11 | sort -u > temp.txt ;
    # num of genes
    num_rOCR=`wc -l temp.txt | awk '{print $1}'`
    num_ubi=`wc -l GRCh38_ubi-rOCR_overlapped_gene_id.txt | awk '{print $1}'`
    # get non-ubi active-rOCRs overlapped genes expression
    awk -v sample="$sample" -v num="$num_rOCR" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"non-ubi-active-rOCR_overlapped_genes",sample,num}}}' \
    temp.txt ./all_gene_exp/${rna_exp_id}.txt >> GRCh38_OCRs_RNA_signal_comparison.txt
    #
    awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"ubi-rOCR_overlapped_genes",sample,num}}}' \
    GRCh38_ubi-rOCR_overlapped_gene_id.txt ./all_gene_exp/${rna_exp_id}.txt >> GRCh38_OCRs_RNA_signal_comparison.txt
done < hg38_matched_DNase_RNA_RAMPAGE_list.txt


# 4. RAMPAGE: ubi-rOCRs vs. non-ubi active-rOCRs overlapped TSSs
if [ -f GRCh38_OCRs_RAMPAGE_signal_comparison.txt ]
then
    rm GRCh38_OCRs_RAMPAGE_signal_comparison.txt
fi
#
while read line
do
    echo ${line};
    dnase_exp_id=`awk '{print $1}' <<< ${line}`
    dnase_file_id=`awk '{print $2}' <<< ${line}`
    rampage_exp_id=`awk '{print $5}' <<< ${line}`
    rampage_file_id=`awk '{print $6}' <<< ${line}`
    sample=`awk '{print $7}' <<< ${line}`
    # get non_ubi_active_OCRs overlapped TSSs
    bedtools intersect -a ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed -b \
    hg38_ubi-rOCR_non_overlapped_TSS_uniqID.bed -wa -wb | cut -f 8 | sort -u > temp.txt ;
    # num of TSSs
    num_rOCR=`wc -l temp.txt | awk '{print $1}'`
    num_ubi=`wc -l ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID_list.txt | awk '{print $1}'`
    # get non_ubi_active_OCRs overlapped TSSs signal
    if [ ${sample} == "GM23338" ]
    then
        sample2="GM23338_53_year"
    elif [ ${sample} == "myotube" ]
    then
        sample2="myotube_NA"
    elif [ ${sample} == "skeletal_muscle_myoblast" ]
    then
        sample2="skeletal_muscle_myoblast_NA"
    else
        sample2=${sample}
    fi
    awk -v sample="$sample" -v num="$num_rOCR" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"non_ubi-active-rOCR_overlapped_TSSs", sample, num}}}' \
    temp.txt ./rampage_tissue/${sample2}_rampage.txt >> GRCh38_OCRs_RAMPAGE_signal_comparison.txt ;
    # ubi-rOCR overlapped TSS
    awk -v sample="$sample" -v num="$num_ubi" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"ubi-rOCR_overlapped_TSSs", sample, num}}}' \
    ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID_list.txt ./rampage_tissue/${sample2}_rampage.txt >> GRCh38_OCRs_RAMPAGE_signal_comparison.txt ;
done < hg38_matched_DNase_RNA_RAMPAGE_list.txt

rm temp.txt


# 5. make boxplot
# Rscipt fig1.signal_comparison_boxplot.R
