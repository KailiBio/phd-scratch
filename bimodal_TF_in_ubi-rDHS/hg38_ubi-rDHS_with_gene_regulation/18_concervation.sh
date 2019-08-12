#!/bin/bash

# -- Kaili
# This script is for comparing sequence conservation between ubi-rOCRs & rOCRs.

# 1. get sequence conservation score
# 2. between ubi-rOCRs & rOCRs
# 3. between ubi-rOCRs & rOCRs all overlap TSSs

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/conservation/"

cd ${workDir}

# 1. get sequence conservation score
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons7way/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.bw
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP7way/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP7way/hg38.phyloP7way.bw

# 2. between ubi-rOCRs & rOCRs
# phastCons
bigWigAverageOverBed hg38.phastCons7way.bw /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed hg38_rOCRs_phastCons7.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1, "ubi-rOCRs",$5}else{print $1,"rOCRs",$5}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed hg38_rOCRs_phastCons7.tab > hg38_rOCRs_phastCons7.txt
# phyloP
bigWigAverageOverBed hg38.phyloP7way.bw /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed hg38_rOCRs_phyloP7.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1, "ubi-rOCRs",$5}else{print $1,"rOCRs",$5}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed hg38_rOCRs_phyloP7.tab > hg38_rOCRs_phyloP7.txt

# 3. between ubi-rOCRs & rOCRs all overlap TSSs
# phastCons
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.bed hg38_rOCRs_phastCons7.txt > hg38_rOCRs_overlap_TSS_phastCons7.txt
# phyloP
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.bed hg38_rOCRs_phyloP7.txt > hg38_rOCRs_overlap_TSS_phyloP7.txt

# Rscript make_conservation_fig.R
