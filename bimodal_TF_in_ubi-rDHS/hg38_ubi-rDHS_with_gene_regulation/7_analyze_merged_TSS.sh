#!/bin/bash

# -- Kaili
# This script is for analysing merged TSS.
# 0. get merged TSS
# 1. RAMPAGE signal in each tissue
# 2. boxplot of rampage signal
# 3. histogram: TSS overlapped with ubi-rOCRs/rOCRs
# 4. histogram: TS index of merged-TSS
# 5.


scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/"

# 0. get merged TSS
sort -k7,7 -k1,1 -k2,2n TSS.Filtered.bed > TSS.Filtered_sorted_gene.bed
python ${scriptDir}merge_TSS_50bp.py
sort -k1,1 -k2,2n hg38_merged_TSS_gene.bed > hg38_merged_TSS_gene_sorted.bed
## used to be 194,897 TSSs in 56,516 genes
## now we have 133,609 merged TSSs.


############################
# Nov 01
# assign merged-TSS to ubi-rOCRs overlapped or not based on percentage of overlapped ubi-rOCRs/rOCRs.
# 1） overlapped with ubi-rOCRs
intersectBed -a hg38_merged_TSS_gene_sorted.bed -b GRCh38_ubi-rOCRs_EDGEid.bed -wa -c | awk '{if($8>0){print $0}}' > \
./merged-TSS/GRC38_merged-TSS_num_overlapped_ubi-rOCR.txt
# 2) overlapped with rOCRs
intersectBed -a hg38_merged_TSS_gene_sorted.bed -b GRCh38-rOCRs.bed -wa -c | awk '{if($8>0){print $0}}' > \
./merged-TSS/GRC38_merged-TSS_num_overlapped_rOCR.txt
# 3) assign more ubi-rOCR overlapped merged-TSS as ubi-rOCR overlapped merged-TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$8}else{if(a[$4]){if((b[$4]/$8)>=0.5){print $4,b[$4],$8,(b[$4]/$8)}}}}' \
./merged-TSS/GRC38_merged-TSS_num_overlapped_ubi-rOCR.txt ./merged-TSS/GRC38_merged-TSS_num_overlapped_rOCR.txt > \
./merged-TSS/GRCh38_merged-TSS_overlapped_count.txt
#
cut -f 1 ./merged-TSS/GRCh38_merged-TSS_overlapped_count.txt | sort -k1,1 > ./merged-TSS/GRCh38_ubi-rOCR_overlapped_merged-TSS_list.txt

############################


cd ${workDir}
## 1) length histogram of merged-TSSs.
awk '{FS=OFS="\t"}{print $0,$3-$2}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene_sorted.bed > GRCh38_merged-TSS_gene_length.txt
# only overlapped ubi-rOCR
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' GRCh38_ubi-rOCR_overlapped_merged-TSS_list.txt \
GRCh38_merged-TSS_gene_length.txt > GRCh38_ubi-rOCR_overlapped_merged-TSS_gene_length.txt
# make histogram
# locally
# Rscript make_merged-TSS_length_histogram.R


## 3) tissue-specificity index of merged-TSSs
### get merged-TSS tissue RAMPAGE signal matrix
echo "mergedTSS_id" > hg38_tissue_mergedTSS_exp_matrix.txt
cut -f 1 /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/A172_rampage.txt >> hg38_tissue_mergedTSS_exp_matrix.txt
#
while read line
do
    id=`awk '{print $1}' <<< ${line}`
    biosample=`awk '{print $4}' <<< ${line}`
    echo ${id} > temp.txt
    awk '{print $2}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/${biosample}_rampage.txt >> temp.txt
    paste hg38_tissue_mergedTSS_exp_matrix.txt temp.txt > temp2.txt
    mv temp2.txt hg38_tissue_mergedTSS_exp_matrix.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_TSS_exp_list2.txt
rm temp.txt

### calculate tissue-specificity index
awk '{if(NR>1){print $0}}' hg38_tissue_mergedTSS_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_mergedTSS_exp_TSscore.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_mergedTSS_exp_TSscore.txt
rm tmp.txt



#-----------------------------------------------------------
# 1. scatterplot of RAMPAGE signal in each tissue

## 1) calculate RAMPAGE signal
awk '{FS=OFS="\t"}{if($6=="+"){print $0}}' GRCh38_merged-TSS_gene_length.txt > GRCh38_merged-TSS_gene_length_plus.txt
awk '{FS=OFS="\t"}{if($6=="-"){print $0}}' GRCh38_merged-TSS_gene_length.txt > GRCh38_merged-TSS_gene_length_minus.txt
#
awk '{FS=OFS="\t"}{print $1,$2-50,$3+50,$4,1,$6}' GRCh38_merged-TSS_gene_length.txt > GRCh38_merged-TSS_gene_50bp_bed6.bed
#
bash ${scriptDir}calculate_RAMPAGE_signal_merged-TSS.sh

## 2) make figures
Rscript ${scriptDir}make_rampage_signal_scatter_between_tissue_separate.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "K562" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_scatter_between_tissue_separate.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "liver_32_year" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_scatter_between_tissue_separate.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "GM12878" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_scatter_between_tissue_separate.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "H7-hESC_NA" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_scatter_between_tissue_separate.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "stomach_40_week" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_scatter_between_tissue_separate.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "lung_24_week" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_scatter_between_tissue_separate.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "spleen_53_year" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"



#-----------------------------------------------------------
# 2. boxplot of rampage signal
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"A172" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"K562" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"liver_32_year" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"GM12878" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"H7-hESC_NA" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"stomach_40_week" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"lung_24_week" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"
#
Rscript ${scriptDir}make_rampage_signal_boxplot_overlapped_merged-TSS.R "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" \
"spleen_53_year" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"



#-----------------------------------------------------------
# 3. histogram: TSS overlapped with ubi-rOCRs/rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' GRCh38_ubi-rOCR_overlapped_merged-TSS_list.txt \
../hg38_merged_TSS_gene.bed > GRCh38_ubi-rOCR_overlapped_merged-TSS.bed
#
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed \
-b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene_sorted.bed -wa -wb | awk '{FS=OFS="\t"}{print $4,$11}' | sort -u | \
cut -f 1 | uniq -c | awk '{FS=" ";OFS="\t"}{print $2,$1}' > GRCh38_ubi-rOCR_overlapped_merged-TSS_count.txt
# rOCRs
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed \
-b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene_sorted.bed -wa -wb | awk '{FS=OFS="\t"}{print $4,$8}' | sort -u | \
cut -f 1 | uniq -c | awk '{FS=" ";OFS="\t"}{print $2,$1}' > GRCh38_rOCR_overlapped_merged-TSS_count.txt
#
# locally
# Rscript make_overlapped_mergedTSS_count_histogram.R



#-----------------------------------------------------------
# 4. histogram: TS index of merged-TSS
## get merged-TSS of housekeeping genes
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]){print $4}}}' /home/fankaili/genome/hg38_housekeeping_geneID_geneType_geneSymbol.txt \
GRCh38_ubi-rOCR_overlapped_merged-TSS.bed | sort -u > GRCh38_HK_mergedTSS_list.txt


## overlapped
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{print a[$1],$1,$2}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed \
hg38_tissue_mergedTSS_exp_TSscore.txt | sort -u > hg38_tissue_mergedTSS_exp_TSscore_withGene.txt
#
cut -f 7 GRCh38_ubi-rOCR_overlapped_merged-TSS.bed | sort -u > GRCh38_ubi-rOCR_overlapped_gene_list.txt
cut -f 4 GRCh38_ubi-rOCR_overlapped_merged-TSS.bed | sort -u > GRCh38_ubi-rOCR_overlapped_mergedTSS_list.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' GRCh38_ubi-rOCR_overlapped_gene_list.txt \
hg38_tissue_mergedTSS_exp_TSscore_withGene.txt | sort -u > hg38_tissue_ubi-rOCR_overlapped_mergedTSS_exp_TSscore_withGene.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$2]){print $0,"overlapped"}else{print $0,"not_overlapped"}}}' \
GRCh38_ubi-rOCR_overlapped_mergedTSS_list.txt hg38_tissue_ubi-rOCR_overlapped_mergedTSS_exp_TSscore_withGene.txt > \
hg38_tissue_ubi-rOCR_overlapped_gene_mergedTSS_TSscore.txt

# Rscript analyze_mergedTSS_TSindex.R

#-----------------------------------------------------------
# 5. gene TS index calculation
## 1) mean of overlapped merged-TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7;b[$4]=1}else{if(b[$1]){print a[$1],$1,$2}}}' GRCh38_ubi-rOCR_overlapped_merged-TSS.bed \
hg38_tissue_mergedTSS_exp_TSscore.txt | sort -k1,1 > GRCh38_ubi-rOCR_overlapped_mergedTSS_TSindex.txt
# calculate mean for each gene
awk '{FS=OFS="\t"}{if(NR==1){id=$1;sum=$3;n=1}else{if(id!=$1){print id,sum/n;id=$1;sum=$3;n=1}else{sum+=$3;n+=1}}}END{print id,sum/n}' \
GRCh38_ubi-rOCR_overlapped_mergedTSS_TSindex.txt > GRCh38_ubi-rOCR_overlapped_mergedTSS_TSindex_average.txt

## 2) mean of all merged-TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{if(FNR!=1){print a[$1],$1,$2}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed \
hg38_tissue_mergedTSS_exp_TSscore.txt | sort -k1,1 > GRCh38_mergedTSS_TSindex.txt
# calculate mean for each gene
awk '{FS=OFS="\t"}{if(NR==1){id=$1;sum=$3;n=1}else{if(id!=$1){print id,sum/n;id=$1;sum=$3;n=1}else{sum+=$3;n+=1}}}END{print id,sum/n}' \
GRCh38_mergedTSS_TSindex.txt > GRCh38_mergedTSS_TSindex_average.txt


## 3) mean of rest merged-TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7;b[$4]=1}else{if(b[$1]!=1){print a[$1],$1,$2}}}' GRCh38_ubi-rOCR_overlapped_merged-TSS.bed \
hg38_tissue_mergedTSS_exp_TSscore.txt | sort -k1,1 > GRCh38_non_ubi-rOCR_overlapped_mergedTSS_TSindex.txt
# calculate mean for each gene
awk '{FS=OFS="\t"}{if(NR==1){id=$1;sum=$3;n=1}else{if(id!=$1){print id,sum/n;id=$1;sum=$3;n=1}else{sum+=$3;n+=1}}}END{print id,sum/n}' \
GRCh38_non_ubi-rOCR_overlapped_mergedTSS_TSindex.txt > GRCh38_non_ubi-rOCR_overlapped_mergedTSS_TSindex_average.txt

## 4) make figures
# Rscript compare_mergedTSS_TSindex.R
