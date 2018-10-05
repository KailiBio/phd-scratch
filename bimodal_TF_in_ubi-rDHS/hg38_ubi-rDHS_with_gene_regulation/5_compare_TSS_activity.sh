#!/bin/bash

# -- Kaili
# This script is for comparing TSS activity. (ubi-rOCRs overlapped with cell-type active OCRs overlapped TSSs)
# 1. get RAMPAGE matched DNase data. (the same donor ID)
# 2. gene expression comparison

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}
