#!/bin/bash

# --Kaili
# This script is for comparing gene expression. (ubi-rOCRs overlapped with cell-type active OCRs overlapped genes)
# 1. get RNA-seq matched DNase data. (the same donor ID)
# 2. gene expression comparison

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. get RNA-seq matched DNase data. (the same donor ID)

## 1) get RNA-seq donor ID list
python ${scriptDir}get_donorID.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_closest_gene_exp_list.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_closest_gene_exp_list_donor.txt

## 2) get DNase-seq donor ID list
python ${scriptDir}get_donorID.py /data/projects/psychencode/Registry/V1/GRCh38/Biosample-Lists/DNase-List.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_DNase_file_list_donor.txt

## 3) matching
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=1;b[$6]=1;c[$5]=$0;d[$6]=$0}else{if(a[$5] && ($5!="---")){print $0,c[$5]}else if(a[$6] && ($6!="---")){print $0,c[$6]}else if(b[$5] && ($5!="---")){print $0,d[$5]}else if(b[$6] && ($6!="---")){print $0,d[$6]}}}' \
hg38_closest_gene_exp_list_donor.txt hg38_DNase_file_list_donor.txt > hg38_RNA_DNase_matched_list.txt

#### get 24 biosamples with matched RNA-seq and DNase-seq.




# 2. gene expression comparison

## 1) get non ubi-rOCRs overlapped TSS and non ubi rOCRs
### get non ubi-rOCRs overlapped TSS
awk '{FS=OFS="\t"}{if($15==0){print $8,$9,$10,$11,$12,$13,$14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | cut -f 4 \
> hg38_ubi-rOCR_overlapped_TSS_list.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]!=1){print $0}}}' hg38_ubi-rOCR_overlapped_TSS_list.txt TSS.Filtered.bed \
> hg38_ubi-rOCR_non_overlapped_TSS_list.txt

### get non ubi rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' GRCh38_ubi-rOCRs.bed GRCh38-rOCRs.bed > \
hg38_non_ubi-rOCRs.bed


## 2) get expression for ubi-rOCR overlapped gene and cell-type active OCR overlapped gene, plot
cut -f 14 ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u > ./closest_gene/GRCh38_ubi-rOCR_closest_gene_list.txt
#
mkdir non_ubi_active_OCR
mkdir gene_exp_comparison_file
mkdir gene_exp_comparison_pdf
#
while read line
do
    echo ${line};
    rna_exp_id=`awk '{FS=OFS="\t"}{print $7}' <<< ${line}` ;
    sample=`awk '{FS=OFS="\t"}{print $9}' <<< ${line}` ;
    dnase_exp_id=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    dnase_file_id=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if(a[$4]){print $0}}}' \
    /data/projects/psychencode/Registry/V1/GRCh38/Signal-Files/${dnase_exp_id}"-"${dnase_file_id}.txt hg38_non_ubi-rOCRs.bed \
    > ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed ;
    #
    bedtools intersect -a ./non_ubi_active_OCR/${dnase_exp_id}_OCR.bed -b hg38_ubi-rOCR_non_overlapped_TSS_list.txt -wa -wb \
    | cut -f 11 | sort -u > temp.txt ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"ubi-rOCR_overlapped"}}}' \
    ./closest_gene/GRCh38_ubi-rOCR_closest_gene_list.txt ./all_gene_exp/${rna_exp_id}.txt > ./gene_exp_comparison_file/${sample}.txt ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"active-rOCR_overlapped"}}}' \
    temp.txt ./all_gene_exp/${rna_exp_id}.txt >> ./gene_exp_comparison_file/${sample}.txt ;
    #
    Rscript ${scriptDir}make_comparison_barplot.R ${sample} ;
done < hg38_RNA_DNase_matched_list.txt
rm temp.txt
