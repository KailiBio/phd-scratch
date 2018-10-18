#!/bin/bash

# -- Kaili
# This script is for making DHS & DHS-centered bins track hub.
# 1. get bigBed file
# 2. make trackhub

workDir="/data/public_html_users/fankaili/IDEAS/bedFile/"

cd ${workDir}
# 1. get bigBed file
awk '{FS=OFS="\t";color="140,140,140"}{if(NR==FNR){a[$4]=$6;b[$4]=1}else{if(b[$4]){if(a[$4]=="CTCF-only"){color="0,176,240"}else if(a[$4]=="Enhancer-like"){color="255,205,0"}else if(a[$4]=="Promoter-like"){color="255,0,0"}};print $1,$2,$3,$4,$3-$2,".",$2,$3,color}}' \
/data/projects/psychencode/Registry/V1/mm10/mm10-EDGEs.bed /data/projects/psychencode/Registry/V1/mm10/mm10-rOCRs.bed > mm10_DHS_bed9.bed
#
bedToBigBed mm10_DHS_bed9.bed /home/fankaili/genome/mm10.chrom.sizes.clean mm10_DHS_bed9.bb

# 2. make trackhub
