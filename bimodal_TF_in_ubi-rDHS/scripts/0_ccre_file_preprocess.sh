#!/bin/bash

# -- Kaili
# This is script to pre-process the ccREs files.

mkdir -p /data/zusers/fankaili/ccre

# 1. get ccREs file
cp /data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-cREs.bed /data/zusers/fankaili/ccre/
cp /data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V4/Registry/V4-mm10/mm10-cREs.bed /data/zusers/fankaili/ccre/
# cp hg38-ccREs.bed

# 2. ID transform between rDHS and ccREs (hg19)
cd /data/zusers/fankaili/ccre/
bedtools intersect -a /data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-cREs.bed -b /data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-rDHS-FDR3.bed -wa -wb > ccREs_ID_transfer.bed
awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$13}' ccREs_ID_transfer.bed > ccREs_ID_transfer_clean.bed

# 3. get ubiquitously rDHS list (hg19)
## ubiquitously means 450 out of 462
## and all the ubi-rDHSs are in hg19-ccRE
cd /data/zusers/fankaili/ccre/
awk '{FS=OFS="\t"}{if($2>=450){print $1}}' /data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V4/Registry/V4-hg19/ccRE-DNase-Biosample-Counts.txt > ubi_ccREs_hg19_list.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){print a[$1]}}}' ccREs_ID_transfer_clean.bed ubi_ccREs_hg19_list.txt > ubi_ccREs_hg19_list_ccreid.txt
