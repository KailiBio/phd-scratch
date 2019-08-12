#!/bin/bash

# -- Kaili
# This script is for comparing conservation between ubi-rOCRs and non-ubi-rOCRs.

workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/conservation/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"

cd ${workDir}

ubi="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
non_ubi="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed"

ubi_TSS="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed"
non_ubi_TSS="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/non-ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed"

# 1. phastCons
cut -f 1-4 ${ubi} > tmp.bed
bigWigAverageOverBed /home/fankaili/genome/hg38.phastCons100way.bw tmp.bed hg38_ubi-rOCRs_phastCons100way.tab
bigWigAverageOverBed /home/fankaili/genome/hg38.phastCons100way.bw ${non_ubi} hg38_non-ubi-rOCRs_phastCons100way.tab
awk '{FS=OFS="\t"}{print $1,$5,"ubi-rOCRs"}' hg38_ubi-rOCRs_phastCons100way.tab > hg38_rOCRs_phastCons100way.txt
awk '{FS=OFS="\t"}{print $1,$5,"non_ubi-rOCRs"}' hg38_non-ubi-rOCRs_phastCons100way.tab >> hg38_rOCRs_phastCons100way.txt

cut -f 1-4 ${ubi_TSS} > tmp_tss.bed
bigWigAverageOverBed /home/fankaili/genome/hg38.phastCons100way.bw tmp_tss.bed hg38_ubi-rOCRs_TSS_phastCons100way.tab
bigWigAverageOverBed /home/fankaili/genome/hg38.phastCons100way.bw ${non_ubi_TSS} hg38_non-ubi-rOCRs_TSS_phastCons100way.tab
awk '{FS=OFS="\t"}{print $1,$5,"ubi-rOCRs"}' hg38_ubi-rOCRs_TSS_phastCons100way.tab > hg38_rOCRs_TSS_phastCons100way.txt
awk '{FS=OFS="\t"}{print $1,$5,"non_ubi-rOCRs"}' hg38_non-ubi-rOCRs_TSS_phastCons100way.tab >> hg38_rOCRs_TSS_phastCons100way.txt

# 2. phyloP
bigWigAverageOverBed /home/fankaili/genome/hg38.phyloP100way.bw tmp.bed hg38_ubi-rOCRs_phyloP100way.tab
bigWigAverageOverBed /home/fankaili/genome/hg38.phyloP100way.bw ${non_ubi} hg38_non-ubi-rOCRs_phyloP100way.tab
awk '{FS=OFS="\t"}{print $1,$5,"ubi-rOCRs"}' hg38_ubi-rOCRs_phyloP100way.tab > hg38_rOCRs_phyloP100way.txt
awk '{FS=OFS="\t"}{print $1,$5,"non_ubi-rOCRs"}' hg38_non-ubi-rOCRs_phyloP100way.tab >> hg38_rOCRs_phyloP100way.txt

bigWigAverageOverBed /home/fankaili/genome/hg38.phyloP100way.bw tmp_tss.bed hg38_ubi-rOCRs_TSS_phyloP100way.tab
bigWigAverageOverBed /home/fankaili/genome/hg38.phyloP100way.bw ${non_ubi_TSS} hg38_non-ubi-rOCRs_TSS_phyloP100way.tab
awk '{FS=OFS="\t"}{print $1,$5,"ubi-rOCRs"}' hg38_ubi-rOCRs_TSS_phyloP100way.tab > hg38_rOCRs_TSS_phyloP100way.txt
awk '{FS=OFS="\t"}{print $1,$5,"non_ubi-rOCRs"}' hg38_non-ubi-rOCRs_TSS_phyloP100way.tab >> hg38_rOCRs_TSS_phyloP100way.txt


# 3. merge and make figures
# Rscript make_conservation_boxplot.R
