#!/bin/bash

# -- Kaili
# This script is for assigning uniq ID to TSS.

# INPUT: 7 column TSS.bed ($4 is transcriptID, $7 is geneID)
# OUTPUT: 7 column TSS.bed with $4 replaced by uniq ID
#                   8 column TSS.bed, adding uniqID to each transcript in $8.
# EXP: bash get_uniq_TSS.sh hg38_v28_basic_TSS_filtered.bed /home/fankaili/genome/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

tss_bed=$1
inPath=$2
workDir=$3

#####################
cd ${workDir}
name=${tss_bed%.bed}

# get uniq TSSs
awk '{FS=OFS="\t"}{print $1,$2,$3,$6,$1"_"$2"_"$3"_"$6}' ${inPath}${tss_bed} | sort -u | sort -k1,1 -k2,2n | \
awk '{FS=OFS="\t"}{print $0,"TSS_"NR}' > tmp.TSS_uniq.bed

# 7 column TSS.bed with $4 replaced by uniq ID
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$6}else{if(a[$1"_"$2"_"$3"_"$6]){print $1,$2,$3,a[$1"_"$2"_"$3"_"$6],$5,$6,$7}}}' \
tmp.TSS_uniq.bed ${inPath}${tss_bed} | sort -u | sort -k1,1 -k2,2n > ${name}_uniq.bed

#  8 column TSS.bed, adding uniqID to each transcript in $8.
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$6}else{print $0,a[$1"_"$2"_"$3"_"$6]}}' \
tmp.TSS_uniq.bed ${inPath}${tss_bed} | sort -u | sort -k1,1 -k2,2n > ${name}_with_uniqID.bed

rm tmp.TSS_uniq.bed
