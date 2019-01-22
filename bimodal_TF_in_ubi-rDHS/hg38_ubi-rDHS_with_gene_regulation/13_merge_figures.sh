#!/bin/bash

# -- Kaili
# This script is for merging all the meaningful figures.
# fig1. boxplot for DNase-seq, RNA-seq, RAMPAGE in 16 identical samples.
# fig2. tissue-specificity index density plot
# fig3.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# fig1. boxplot for DNase-seq, RNA-seq, RAMPAGE in 16 identical samples.
bash ${scriptDir}fig1_DNase_RNA_RAMPAGE_signal_comparison_boxplot.sh

# fig2. tissue-specificity index density plot
