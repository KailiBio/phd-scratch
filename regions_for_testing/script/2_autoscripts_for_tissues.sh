#!/bin/bash

# -- Kaili
# This script is for making and running automatic pipeline to get high density enhancer tables for all 8 tissues in mouse e11.5

workDir="/data/zusers/fankaili/regions/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/regions_for_testing/script/"

ccreDir="/data/projects/encode/Registry/V2/mm10/"
signalDir="/data/projects/encode/Registry/V2/mm10/Signal-Files/"
rna_filelist="/data/zusers/fankaili/regions/mouse_e11.5_RNA-seq_list.txt"

cd ${workDir}

# 0. preparation
cat ./e11.5_dELS/e11.5_forebrain_dELS.bed ./e11.5_dELS/e11.5_midbrain_dELS.bed ./e11.5_dELS/e11.5_hindbrain_dELS.bed ./e11.5_dELS/e11.5_neural-tube_dELS.bed ./e11.5_dELS/e11.5_heart_dELS.bed ./e11.5_dELS/e11.5_facial_dELS.bed ./e11.5_dELS/e11.5_limb_dELS.bed ./e11.5_dELS/e11.5_liver_dELS.bed | sort -u | sort -k1,1 -k2,2n > ./e11.5_dELS/e11.5_dELS.bed

# 1. get split table
while read line
do
    sample=`awk '{print $1}' <<< $line`
    echo $sample
    #
    nohup bash ${scriptDir}generate_high_density_enhancer_given_tissue_e11.5.sh ${sample} > ./nohup/nohup.generate_high_density_enhancer_given_${sample}_e11.5.out 2>&1&
done < ${rna_filelist}

# 2. get merged bins
while read line
do
    sample=`awk '{print $1}' <<< $line`
    echo $sample
    #
    nohup bash ${scriptDir}generate_table_merged_region.sh ${sample} > ./nohup/nohup.generate_table_merged_region_${sample}_e11.5.out 2>&1&
done < ${rna_filelist}
####
# Rscript make_merge_length_histogram_merged.R

# 3. merge 8 tissues into one
bash ${scriptDir}get_high_density_table_for_merged_all_tissue.sh
