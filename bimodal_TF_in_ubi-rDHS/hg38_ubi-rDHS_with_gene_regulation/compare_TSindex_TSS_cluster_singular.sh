#!/bin/bash

# -- Kaili
# This script is for comparing TSindex between TSS-clusters and singular TSSs.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/"

cd ${workDir}

# 1. divide TSSs into TSS-cluster and sigular-TSS
awk '{FS=OFS="\t"}{if($2==$3){print $0}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed > hg38_mergedTSS_singular.bed
awk '{FS=OFS="\t"}{if($2!=$3){print $0}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed > hg38_mergedTSS_cluster.bed

# 2. get corresponding TSindex
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$2]){print $0,"cluster"}else{print $0,"singular"}}}' hg38_mergedTSS_cluster.bed \
GRCh38_mergedTSS_TSindex.txt > GRCh38_mergedTSS_TSindex_cluster_singular.txt

# 3. comparison
# Rscript compare_TSS_TSindex_cluster_singular.R
