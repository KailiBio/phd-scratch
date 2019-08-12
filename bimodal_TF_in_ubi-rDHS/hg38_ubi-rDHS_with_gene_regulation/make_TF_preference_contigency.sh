#!/bin/bash

# -- Kaili
# This script is for making contigency table using different groups of rOCRs.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/"

cd ${workDir}

# 1. all rOCRs
for cellline in K562 GM12878 HepG2 H1
do
    echo -e "motif\tubi_TF\tnon_ubi_TF\tubi\tnon_ubi" > TFpeak_contigency_table_${cellline}.txt
    while read line
    do
        tf=`awk '{print $1}' <<< $line`
        expID=`awk '{print $2}' <<< $line`
        fileID=`awk '{print $3}' <<< $line`
        echo ${tf}
        #
        cp /data/projects/encode/data/${expID}/${fileID}.bed.gz ./
        gzip -d ${fileID}.bed.gz
        #
        ubi_TFpeak=`intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -b ${fileID}.bed -F 0.5 -wa -u | wc -l`
        non_ubi_TFpeak=`intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -b ${fileID}.bed -F 0.5 -wa -u | wc -l`
        ubi=`wc -l /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed | awk '{print $1}'`
        non_ubi=`wc -l /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed | awk '{print $1}'`
        #
        echo -e ${tf}"\t"${ubi_TFpeak}"\t"${non_ubi_TFpeak}"\t"${ubi}"\t"${non_ubi} >> TFpeak_contigency_table_${cellline}.txt
        #
        rm ${fileID}.bed*
    done < ${cellline}_TF_peak_filelist.txt
done

# 2. active rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR && $2>1.64){a[$1]=1}else{if(a[$4]){print $0}}}' /data/zusers/moorej3/moorej.ghpcc.project/PsychENCODE/Registry/V0-hg38/signal-output/ENCSR921NMD-ENCFF971AHO.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed > K562_active_non-ubi-rOCRs.bed
awk '{FS=OFS="\t"}{if(NR==FNR && $2>1.64){a[$1]=1}else{if(a[$4]){print $0}}}' /data/zusers/moorej3/moorej.ghpcc.project/PsychENCODE/Registry/V0-hg38/signal-output/ENCSR000ENP-ENCFF205TKQ.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed > HepG2_active_non-ubi-rOCRs.bed
awk '{FS=OFS="\t"}{if(NR==FNR && $2>1.64){a[$1]=1}else{if(a[$4]){print $0}}}' /data/zusers/moorej3/moorej.ghpcc.project/PsychENCODE/Registry/V0-hg38/signal-output/ENCSR000EMU-ENCFF131HMO.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed > H1_active_non-ubi-rOCRs.bed
#
for cellline in K562 HepG2 H1
do
    echo -e "motif\tubi_TF\tnon_ubi_TF\tubi\tnon_ubi" > TFpeak_contigency_table_active_${cellline}.txt
    while read line
    do
        tf=`awk '{print $1}' <<< $line`
        expID=`awk '{print $2}' <<< $line`
        fileID=`awk '{print $3}' <<< $line`
        echo ${tf}
        #
        cp /data/projects/encode/data/${expID}/${fileID}.bed.gz ./
        gzip -d ${fileID}.bed.gz
        #
        ubi_TFpeak=`intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -b ${fileID}.bed -F 0.5 -wa -u | wc -l`
        non_ubi_TFpeak=`intersectBed -a ${cellline}_active_non-ubi-rOCRs.bed -b ${fileID}.bed -F 0.5 -wa -u | wc -l`
        ubi=`wc -l /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed | awk '{print $1}'`
        non_ubi=`wc -l ${cellline}_active_non-ubi-rOCRs.bed | awk '{print $1}'`
        #
        echo -e ${tf}"\t"${ubi_TFpeak}"\t"${non_ubi_TFpeak}"\t"${ubi}"\t"${non_ubi} >> TFpeak_contigency_table_active_${cellline}.txt
        #
        rm ${fileID}.bed*
    done < ${cellline}_TF_peak_filelist.txt
done


# 3. rOCRs overlapping TSSs
for cellline in K562 HepG2 H1
do
    echo -e "motif\tubi_TF\tnon_ubi_TF\tubi\tnon_ubi" > TFpeak_contigency_table_overlapTSS_${cellline}.txt
    while read line
    do
        tf=`awk '{print $1}' <<< $line`
        expID=`awk '{print $2}' <<< $line`
        fileID=`awk '{print $3}' <<< $line`
        echo ${tf}
        #
        cp /data/projects/encode/data/${expID}/${fileID}.bed.gz ./
        gzip -d ${fileID}.bed.gz
        #
        ubi_TFpeak=`intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed -b ${fileID}.bed -F 0.5 -wa -u | wc -l`
        intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/non-ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed -b ${cellline}_active_non-ubi-rOCRs.bed -wa -u > tmp.bed
        non_ubi_TFpeak=`intersectBed -a tmp.bed -b ${fileID}.bed -F 0.5 -wa -u | wc -l`
        ubi=`wc -l /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed | awk '{print $1}'`
        non_ubi=`wc -l tmp.bed | awk '{print $1}'`
        #
        echo -e ${tf}"\t"${ubi_TFpeak}"\t"${non_ubi_TFpeak}"\t"${ubi}"\t"${non_ubi} >> TFpeak_contigency_table_overlapTSS_${cellline}.txt
        #
        rm ${fileID}.bed*
    done < ${cellline}_TF_peak_filelist.txt
done

Rscript ${scriptDir}make_volcano_cellline.R
