#!/bin/bash

# -- Kaili
# This script is using vM4 result as example for sanity check.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

# 1. get TPM matrix (vM4)
# from 3_call_DE.sh - step one
# mm10_M4_TPM_matrix.txt

 # 2. do clustering
#Rscript do_cluster_M4.R
