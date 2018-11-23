#!/bin/bash

# -- Kaili
# This script is for getting RAMPAGE signal weighted TSS TS-index as gene TS-index.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/"

cd ${workDir}

# 1. get average RAMPAGE signal for each mergedTSS
awk '{FS=OFS="\t"}{if(NR>1){sum=0;for(i=2;i<=NF;i++){sum+=$i};print $1,sum/(NF-1)}}' hg38_tissue_mergedTSS_exp_matrix.txt \
> hg38_mergedTSS_ave_RAMPAGE.txt


# 2. for each gene, get TSS TSindex weight.

# get aveRAMPAGE with matched TSindex
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' hg38_mergedTSS_ave_RAMPAGE.txt GRCh38_mergedTSS_TSindex.txt \
> hg38_mergedTSS_TSindex_aveRAMPAGE.txt

# get weighted for each mergedTSS
awk '{FS=OFS="\t"}{if(NR==1){gene=$1;sum=$4}else{if(gene==$1){sum+=$4}else{print gene,sum;gene=$1;sum=$4}}}' \
hg38_mergedTSS_TSindex_aveRAMPAGE.txt > hg38_mergedTSS_geneSumRAMPAGE.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(a[$1]==0){print $0,0}else{print $0,$4/a[$1]}}}' hg38_mergedTSS_geneSumRAMPAGE.txt \
hg38_mergedTSS_TSindex_aveRAMPAGE.txt > hg38_mergedTSS_TSindex_aveRAMPAGE_weight.txt


# 3. get weigted TSindex
awk '{FS=OFS="\t"}{if(NR==1){gene=$1;sum=($3*$5)}else{if(gene==$1){sum+=($3*$5)}else{print gene,sum;gene=$1;sum=($3*$5)}}}' \
hg38_mergedTSS_TSindex_aveRAMPAGE_weight.txt > hg38_mergedTSS_gene_weigedt_TSindex.txt



# 4. comparison
# Rscript compare_geneTSindex_weightedTSindex.R
