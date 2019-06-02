#!/bin/bash

# -- Kaili
# This script is for comparing TF preference between ubi-rOCRs and non-ubi-rOCRs.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/"

cd ${workDir}

# 1. get encode TF ChIP-seq peak datalist.
