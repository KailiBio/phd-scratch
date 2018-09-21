#!/bin/bash

# -- Kaili
# This script is for doing ubi-rDHS analysis on hg38.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scripts/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. get ubi-rDHS in hg38

# 2. get closest gene list

cp /data/zusers/moorej3/moorej.ghpcc.project/Reference/Human/hg38/GENCODE24/TSS.Filtered.bed ./
#


## basic info


# 3. get TSS expression data for hg38 TSS
