#!/bin/bash/

# -- Kaili
# This script is for getting DHS center 200-bp bins.
# 1. check mm10 DHS sizes and DHS-gap sizes
# 2. version 1: 100-300bp bins
# 3. version 2: all regions

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/"

cd ${workDir}


# 1. check mm10 DHS sizes and DHS-gap sizes
## get genomic region Yu used.
awk '{FS=" ";OFS="\t"}{print $1,$2,$3}' /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed | sort -k1,1 -k2,2n > mm10.bed
bedtools merge -i mm10.bed > mm10_chrom.bed
grep -v "chrM" mm10_chrom.bed > mm10_chrom_clean.bed
# get gap region
sort -k1,1 -k2,2n /data/projects/psychencode/Registry/V1/mm10/mm10-rOCRs.bed | cut -f 1-3 > mm10-rOCRs.bed
bedtools subtract -a mm10_chrom_clean.bed -b mm10-rOCRs.bed > mm10_rOCR_gap.bed

# Rscript mm10_DHS_length_basic.R


# 2. version 1: 100-300bp bins
python ${scriptDir}DHS-center_bins_v1.py /data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed /data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap.bed \
/data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v1.bed


# 3. version 2: all regions, 1-300bp bins
python ${scriptDir}DHS-center_bins_v2.py /data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed /data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap.bed \
/data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v2.bed


# 4. version 3: all regions, 150-400bp bins
awk 'BEGIN{FS=OFS="\t";n=1;chr="chr1";end=0}{if($1==chr){gap=($2-end)}else{gap=($2-0)};print $0,"O"n,n,$3-$2,gap;n+=1;chr=$1;end=$3}' \
/data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed > mm10-rOCRs.bed_withLength.txt
#
python ${scriptDir}merging_DHS_gap.py
sort -k5,5n mm10-rOCRs_v3.txt | awk '{FS=OFS="\t"}{print $0,$3-$2,$3-$2-$6}' > mm10-rOCRs_v3_sorted.txt
#
sort -k1,1 -k2,2n mm10-rOCRs_v3_sorted.txt | cut -f 1-3 > temp.txt
bedtools subtract -a mm10_chrom_clean.bed -b temp.txt > mm10_rOCR_gap_v3.bed
## manage small gaps
awk '{FS=OFS="\t"}{if(($3-$2)<100){print $0,$3-$2}}' mm10_rOCR_gap_v3.bed > temp2.txt
cp mm10-rOCRs_v3_sorted.txt mm10-rOCRs_v3_sorted_2.txt
while read line
do
    a=`awk '{print $2}' <<< ${line}`
    b=`awk '{print $3}' <<< ${line}`
    c=`awk '{print $1}' <<< ${line}`
    num_a=`grep ${c}$'\t' mm10-rOCRs_v3_sorted_2.txt | grep $'\t'${a}$'\t' | wc -l`
    num_b=`grep ${c}$'\t' mm10-rOCRs_v3_sorted_2.txt | grep $'\t'${b}$'\t' | wc -l`
    if [ ${num_b} -ne 0 ]
    then
        sed -i "s/\t$b\t/\t$a\t/" mm10-rOCRs_v3_sorted_2.txt
    elif [ ${num_a} -ne 0 ]
    then
        sed -i "s/\t$a\t/\t$b\t/" mm10-rOCRs_v3_sorted_2.txt
    fi
done < temp2.txt
#
sort -k5,5n mm10-rOCRs_v3_sorted_2.txt | awk '{FS=OFS="\t"}{print $0,$3-$2,$3-$2-$8}' > mm10-rOCRs_v3_sorted_3.txt
sort -k1,1 -k2,2n mm10-rOCRs_v3_sorted_2.txt | cut -f 1-3 > temp3.txt
bedtools subtract -a mm10_chrom_clean.bed -b temp3.txt > mm10_rOCR_gap_v3_2.bed
# cut gaps into bins
python ${scriptDir}DHS-center_bins_v3.py
#
cut -f 1-4 mm10-rOCRs_v3_sorted_2.txt > mm10_OCR-center_bins_v3.bed
cat mm10_OCR-center_bins_v3_gap.bed >> mm10_OCR-center_bins_v3.bed


# Rscript make_length_distribution_dhs-bins.R
# Rscript make_DHS_bins_v3_length_distribution.R
